#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <zlib.h>

#ifdef USE_SAMTOOLS_LIBS
#include "samtools/khash.h"
#include "samtools/kseq.h"
#else
#include "khash.h"
#include "kseq.h"
#endif

#ifndef _LIB_ONLY
#define _SEQQS_MAIN
#endif

KSEQ_INIT(gzFile, gzread)
KHASH_MAP_INIT_STR(str, uint64_t)

#define INIT_SEQLEN 10
#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef enum {
  PHRED,
  SANGER,
  SOLEXA,
  ILLUMINA,
  NONE
} qual_type;

#define Q_OFFSET 0
#define Q_MIN 1
#define Q_MAX 2

static const int quality_contants[4][3] = {
  /* offset, min, max */
  {0, 4, 60}, /* PHRED */
  {33, 0, 93}, /* SANGER */
  {64, -5, 62}, /* SOLEXA */
  {64, 0, 62} /* ILLUMINA */
};

#define qrng(type) (quality_contants[(type)][Q_MAX] - quality_contants[(type)][Q_MIN] + 1)
#define qmin(type) (quality_contants[(type)][Q_MIN])
#define qmax(type) (quality_contants[(type)][Q_MAX])
#define qoffset(type) (quality_contants[(type)][Q_OFFSET])

#define CHAR_WIDTH 256

/* nucleotide table; all non-IUPAC codes go to X */
unsigned char seq_nt17_table[256] = {
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,16,16,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,1,14,2, 13,0,0,4, 11,0,0,12, 0,3,15,0,
    0,0,5,6, 8,0,7,9, 0,10,0,0, 0,0,0,0,
    0,1,14,2, 13,0,0,4, 11,0,0,12, 0,3,15,0,
    0,0,5,6, 8,0,7,9, 0,10,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

char *seq_nt17_rev_table = "XACMGRSVTWYHKDBN-";

typedef struct _qs_set_t {
  size_t l, m;
  unsigned k;
  uint64_t **ntm;
  uint64_t **qm;
  uint64_t *lm;
  qual_type qt;
  uint64_t n_uniq_kmer_pos;
  khash_t(str) *h;
} qs_set_t;

static void qs_printstr(const kstring_t *s, unsigned line_len) {
  if (line_len != UINT_MAX) {
    int i, rest = s->l;
    for (i = 0; i < s->l; i += line_len, rest -= line_len) {
      putchar('\n');
      if (rest > line_len) fwrite(s->s + i, 1, line_len, stdout);
      else fwrite(s->s + i, 1, rest, stdout);
    }
    putchar('\n');
  } else {
    putchar('\n');
    puts(s->s);
  }
}

void qs_printseq(const kseq_t *s, int line_len) {
  putchar(s->qual.l? '@' : '>');
  fputs(s->name.s, stdout);
  if (s->comment.l) {
    putchar(' '); fputs(s->comment.s, stdout);
  }
  qs_printstr(&s->seq, line_len);
  if (s->qual.l) {
    putchar('+');
    qs_printstr(&s->qual, line_len);
  }
}


qs_set_t *qs_init(qual_type qt, unsigned k) {
  /* 
     Allocate matrices for quality and nucleotides. Rows correspond to
     position in sequence, so growing them is simpler.
  */
  unsigned i;
  qs_set_t *qs = malloc(sizeof(qs_set_t));
  qs->qt = qt;
  qs->m = (size_t) INIT_SEQLEN;
  qs->l = 0;
  qs->n_uniq_kmer_pos = 0;
  qs->k = k;
  qs->h = k > 0 ? kh_init(str) : NULL;
  qs->ntm = malloc(qs->m*sizeof(uint64_t*));
  qs->qm = malloc(qs->m*sizeof(uint64_t*));
  qs->lm = calloc(qs->m, sizeof(uint64_t));
  for (i = 0; i < qs->m; i++) {
    qs->ntm[i] = calloc(17, sizeof(uint64_t));
    qs->qm[i] = calloc(qrng(qs->qt), sizeof(uint64_t));
  }
  return qs;
}

void qs_update(qs_set_t *qs, kseq_t *seq, int strict) {
  unsigned i, last_m, nt, non_iupac=0;
  char bq;
  if (seq->seq.l > qs->m) {
    /* grow all matrices */
    last_m = qs->m;
    qs->m = seq->seq.l + 1;
    kroundup32(qs->m);
    /* fprintf(stderr, "[%s] adding rows to matrix (old size: %d; new size: %d)\n", __func__, last_m, qs->m); */
    qs->ntm = realloc(qs->ntm, sizeof(uint64_t*)*qs->m);
    qs->qm = realloc(qs->qm, sizeof(uint64_t*)*qs->m);
    qs->lm = realloc(qs->lm, sizeof(uint64_t)*qs->m);
    for (i = last_m; i < qs->m; i++) {
      qs->ntm[i] = calloc(17, sizeof(uint64_t));
      qs->qm[i] = calloc(qrng(qs->qt), sizeof(uint64_t));
    }
  }

  /* update largest sequence encountered */
  if (seq->seq.l > qs->l) qs->l = seq->seq.l;
  
  /* update length (0-indexed) */
  qs->lm[seq->seq.l-1]++;

  char *kmer;
  khiter_t key;
  int is_missing, ret;
  if (qs->k) kmer = malloc(sizeof(char)*(qs->k + 2 + log10(UINT_MAX)));

  if (seq->qual.l && seq->seq.l != seq->qual.l) {
    fprintf(stderr, "[%s] warning: quality and sequence lengths differ in sequence '%s'\n", __func__, seq->name.s);
    if (strict) exit(1);
  }

  for (i = 0; i < seq->seq.l; i++) {
    /* update nucleotide composition */
    nt = seq_nt17_table[(int) seq->seq.s[i]];
    if (!nt) non_iupac++;
    qs->ntm[i][nt]++;
    
    /* update quality composition */
    if (seq->qual.l) {
      bq = (char) seq->qual.s[i];
      if (bq - qoffset(qs->qt) < qmin(qs->qt) || 
	  bq - qoffset(qs->qt) > qmax(qs->qt)) {
	fprintf(stderr, "[%s] warning: base quality '%d' out of range (%d <= b <= %d) in sequence '%s'\n", __func__, bq, qmin(qs->qt), qmax(qs->qt), seq->name.s);
	if (strict) exit(1);
      }
      qs->qm[i][bq - qoffset(qs->qt) - qmin(qs->qt)]++;
    }
    
    /* hash positional k-mers */
    if (qs->k) {
      if (qs->k > seq->seq.l) {
	fprintf(stderr, "[%s] warning: k-mer length longer than sequence '%s'\n", __func__, seq->name.s);
      } else {
	if (i <= seq->seq.l-qs->k) {
	  strncpy(kmer, seq->seq.s, (size_t) qs->k);
	  sprintf(kmer + qs->k, "-%u", i+1);
	  
	  /* hash kmer */
	  key = kh_get(str, qs->h, kmer);
	  is_missing = (key == kh_end(qs->h));
	  if (is_missing) {
	    key = kh_put(str, qs->h, strdup(kmer), &ret);
	    kh_value(qs->h, key) = 1;
	    qs->n_uniq_kmer_pos++;
	  } else {
	    kh_value(qs->h, key) = kh_value(qs->h, key) + 1;
	  }
	}
      }
    }
  }
  if (qs->k) free(kmer);
  
  if (non_iupac) {
    fprintf(stderr, "[%s] warning: %d non-IUPAC characters found in sequence '%s'.\n", __func__, non_iupac, seq->name.s);
    if (strict) exit(1);
  }
}

void qs_qm_fprint(FILE *file, qs_set_t *qs) {
  unsigned i, j;
  uint64_t cnt;
  char bq;

  for (j = 0; j < qrng(qs->qt); j++) {
    bq = j + qoffset(qs->qt) + qmin(qs->qt);
    fprintf(file, "Q%d", bq);
    if (j < qrng(qs->qt)-1) fputc('\t', file);
  }
  fputc('\n', file);

  for (i = 0; i < qs->l; i++) {
    for (j = 0; j < qrng(qs->qt); j++) {
      cnt = qs->qm[i][j];
      fprintf(file, "%llu", cnt);
      if (j < qrng(qs->qt)-1) fputc('\t', file);
    }
    fputc('\n', file);
  }
  fputc('\n', file);
}

void qs_ntm_fprint(FILE *file, qs_set_t *qs) {
  unsigned i, j;
  uint64_t cnt;

  for (j = 0; j < 17; ++j) {
    fputc(seq_nt17_rev_table[j], file);
    if (j < 16) fputc('\t', file);
  }
  fputc('\n', file);

  for (i = 0; i < qs->l; i++) {
    for (j = 0; j < 17; j++) {
      cnt = qs->ntm[i][j];
      fprintf(file, "%llu", cnt);
      if (j < 16) fputc('\t', file);
    }
    fputc('\n', file);
  }
  fputc('\n', file);
}

void qs_lm_fprint(FILE *file, qs_set_t *qs) {
  unsigned i;
  for (i = 0; i < qs->l; i++) {
    fprintf(file, "%d\t%llu\n", i+1, qs->lm[i]);
  }
  fputc('\n', file);
}

void qs_kmer_fprint(FILE *file, qs_set_t *qs) {
  if (!qs->k) return;
  khiter_t k;
  char *key;
  int i = 0;
  fprintf(file, "kmer\tpos\tcount\n");
  for (k = kh_begin(qs->h); k != kh_end(qs->h); ++k) {
    if (!kh_exist(qs->h, k)) continue;
    key = (char*) kh_key(qs->h, k);
    fprintf(file, "%.*s\t%s\t%llu\n", qs->k, key, key+(qs->k+1), kh_value(qs->h, k));
    i++;
    free((char *) kh_key(qs->h, k));
  }
  fputc('\n', file);
}


void qs_destroy(qs_set_t *qs) {
  unsigned i;
  for (i = 0; i < qs->m; i++) {
    free(qs->ntm[i]);
    free(qs->qm[i]);
  }
  kh_destroy(str, qs->h);
  free(qs->ntm);
  free(qs->qm);
  free(qs->lm);
  free(qs);
}

int is_interleaved_pair(const char *s1, const char *s2) {
  while (*s1 && *s2) {
    /* strings should be identical apart from trailing 1 or 2
       (i.e. seq-a/1 and seq-a/2) */
    if (*s1 != *s2 && (*s1 != '1' && *s2 != '2'))
      return 0;
    s1++; s2++;
  }
  return 1;
}

#ifdef _SEQQS_MAIN

int usage() {
  fputs("\
Usage: seqqs [options] <in.fq>\n\n\
Options: -q    quality type, either illumina, solexa, or sanger (default: sanger)\n\
         -p    prefix for output files (default: none)\n\
         -k    hash k-mers of length k (default: off)\n\
         -e    emit reads to stdout, for pipelining (default: off)\n\
         -i    input is interleaved, output statistics per each file (default: off)\n\
         -s    strict; some warnings become errors (default: off)\n\
Arguments:  <in.fq> or '-' for stdin.\n\n\
Output:\n\
<prefix> is output prefix name. The following output files will be created:\n\
<prefix>_qual.txt:  quality distribution by position matrix\n\
<prefix>_nucl.txt:  nucleotide distribution by position matrix\n\
<prefix>_len.txt:   length distribution by position matrix\n\
<prefix>_kmer.txt:  k-mer distribution by position matrix\n\
\
If -i is used, these will have \"_1.txt\" and \"_2.txt\" suffixes.", stderr);
  return 1;
}

int main(int argc, char *argv[]) {
  int c, pr=0, k=0, emit=0, strict=0, interleaved=0;
  char *qual_fn="qual.txt", *nucl_fn="nucl.txt", *len_fn="len.txt", *kmer_fn="kmer.txt";
  char *prefix="", *rname;
  FILE *qual_fp[2], *nucl_fp[2], *len_fp[2], *kmer_fp[2];
  qual_type qtype=SANGER;
  qs_set_t *qs[2];
  gzFile fp;
  kseq_t *seq;

  if (argc == 1) return usage();

  while ((c = getopt(argc, argv, "q:k:p:esi")) >= 0) {
    switch (c) {
    case 'q':
      if (strcmp(optarg, "illumina") == 0)
	qtype = ILLUMINA;
      else if (strcmp(optarg, "solexa") == 0)
	qtype = SOLEXA;
      else if (strcmp(optarg, "sanger") == 0)
	qtype = SANGER;
      else {
	fprintf(stderr, "Unknown quality type '%s'.\n", optarg);
	return(1);
      }
      break;
    case 'k':
      k = atoi(optarg);
      break;
    case 's':
      strict = 1;
      break;
    case 'i':
      interleaved = 1;
      break;
    case 'e':
      emit = 1;
      break;
    case 'p':
      prefix = calloc(strlen(optarg)+2, sizeof(char));
      sprintf(prefix, "%s_", optarg);
      break;
    case 'h':
    default:
      return usage();
    }
  }

  if (argc == optind) return usage();
  fp = strcmp(argv[optind], "-") ? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");

  int n = strlen(prefix) + 2*interleaved + 9;
  qual_fn = calloc(n, sizeof(char));
  nucl_fn = calloc(n, sizeof(char));
  len_fn = calloc(n, sizeof(char));
  kmer_fn = calloc(n, sizeof(char));
  
  char *suffix = calloc(7, sizeof(char));
  for (pr = 0; pr < interleaved+1; pr++) {
    if (interleaved)
      sprintf(suffix, "_%d.txt", pr+1);
    else
      sprintf(suffix, ".txt");
    
    sprintf(qual_fn, "%squal%s", prefix, suffix);
    sprintf(nucl_fn, "%snucl%s", prefix, suffix);
    sprintf(len_fn, "%slen%s", prefix, suffix);
    sprintf(kmer_fn, "%skmer%s", prefix, suffix);
    
    qual_fp[pr] = fopen(qual_fn, "w");
    nucl_fp[pr] = fopen(nucl_fn, "w");
    len_fp[pr] = fopen(len_fn, "w");
    kmer_fp[pr] = fopen(kmer_fn, "w");
    
    if (!(qual_fp[pr] && nucl_fp[pr] && len_fp[pr] && kmer_fp[pr])) {
      fprintf(stderr, "[%s] error: cannot open a file for output.\n", __func__);
      return(1);
    }
  }
  free(suffix); 
  if (strlen(prefix)) free(prefix);

  if (strlen(prefix)) {
    free(qual_fn); free(nucl_fn); free(len_fn); free(kmer_fn);
  }
  
  qs[0] = qs_init(qtype, k);
  if (interleaved) {
    qs[1] = qs_init(qtype, k);
  }
  seq = kseq_init(fp);

  while (kseq_read(seq) >= 0) {
    if (strict) rname = strdup(seq->seq.s);
    qs_update(qs[0], seq, strict);
    if (emit) qs_printseq(seq, seq->seq.l);
    rname = strdup(seq->name.s);

    /* for interleaved files, grab and process another entry */
    if (interleaved) {
      if (kseq_read(seq) >= 0) {
	qs_update(qs[1], seq, strict);
	if (emit) qs_printseq(seq, seq->seq.l);
	if (!is_interleaved_pair(seq->name.s, rname)) {
	  fprintf(stderr, "[%s] warning: interleaved reads names differ '%s' != '%s'\n", __func__, seq->name.s, rname);
	  if (strict) return 1;
	}
      } else {
	fprintf(stderr, "[%s] error: interleaved file length not multiple of two.", __func__);
	return 1;
      }
    }
    free(rname);
  }
  
  for (pr = 0; pr < interleaved+1; pr++) {
    qs_ntm_fprint(nucl_fp[pr], qs[pr]);
    qs_qm_fprint(qual_fp[pr], qs[pr]);
    qs_lm_fprint(len_fp[pr], qs[pr]);
    if (k) qs_kmer_fprint(kmer_fp[pr], qs[pr]);
    qs_destroy(qs[pr]);
    fclose(qual_fp[pr]); fclose(nucl_fp[pr]); fclose(len_fp[pr]); fclose(kmer_fp[pr]);
  }

  kseq_destroy(seq);
  gzclose(fp);
}
#endif /* _SEQQS_MAIN */
