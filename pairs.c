#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <limits.h>
#include <unistd.h>
#include <string.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

static int is_interleaved_pair(const char *s1, const char *s2) {
  while (*s1 && *s2) {
    /* strings should be identical apart from trailing 1 or 2
       (i.e. seq-a/1 and seq-a/2) 

       TODO: this checking is fairly simplistic, could be made better.
    */
    if (*s1 != *s2 && (*s1 != '1' && *s2 != '2'))
      return 0;
    s1++; s2++;
  }
  return 1;
}

static int usage() {
  fprintf(stderr, "\nInterleaves (pairs) and un-interleaves paired-end files");
  fprintf(stderr, "Usage <command> <arguments>\n\n");
  fprintf(stderr, "Command:  join       common transformation of FASTA/Q\n");
  fprintf(stderr, "          split      get the nucleotide composition of FASTA/Q\n");
  return 1;
}

static void printstr(FILE *stream, const kstring_t *s, unsigned line_len)
{
  /* from Heng's stk_printstr */
  if (line_len != UINT_MAX) {
    int i, rest = s->l;
    for (i = 0; i < s->l; i += line_len, rest -= line_len) {
      fputc('\n', stream);
      if (rest > line_len) fwrite(s->s + i, 1, line_len, stream);
      else fwrite(s->s + i, 1, rest, stream);
    }
    fputc('\n', stream);
  } else {
    fputc('\n', stream);
    fputs(s->s, stream);
  }
}	

void printseq(FILE *stream, const kseq_t *s, int line_len, int tag) {
  /* from Heng's seqtk */
  fputc(s->qual.l? '@' : '>', stream);
  fputs(s->name.s, stream);
  if (tag) {
    fputc('/', stream); putc((char)(((int)'0')+tag), stream);
  }

  if (s->comment.l) {
    fputc(' ', stream); fputs(s->comment.s, stream);
  }
  printstr(stream, &s->seq, line_len);
  if (s->qual.l) {
    fputc('+', stream);
    printstr(stream, &s->qual, line_len);
  }
}

int join_usage() {
  fputs("\
Usage:    pairs join [options] <in1.fq> <in2.fq>\n\n\
Options:  -t   tag interleaved pairs with '/1' and '/2' (before comment)\n\
          -s   error out when read names are different\n\
Interleaves two paired-end files.\n\n", stderr);
  return 1;
}

int pairs_join(int argc, char *argv[]) {
  gzFile fp[2];
  kseq_t *ks[2];
  int c, i, tag=0, strict=0, l[] = {0, 0};
  while ((c = getopt(argc, argv, "ts")) >= 0) {
    switch (c) {
    case 't': tag = 1; break;
    case 's': strict = 1; break;
    default: return 1;
    }
  }

  if (optind == argc) return join_usage();

  for (i = 0; i < 2; ++i) {
    fp[i] = gzopen(argv[optind + i], "r");
    ks[i] = kseq_init(fp[i]);
  }
  for (;;) {
    for (i = 0; i < 2; ++i) l[i] = kseq_read(ks[i]);
    if (l[0] < 0 || l[1] < 0)
      break;
    if (!is_interleaved_pair(ks[0]->name.s, ks[1]->name.s)) {
      fprintf(stderr, "[%s] warning: different sequence names: %s != %s\n", __func__, ks[0]->name.s, ks[1]->name.s);
      if (strict) return 1;
    }
   
    for (i = 0; i < 2; ++i) printseq(stdout, ks[i], ks[i]->seq.l, tag ? i+1 : 0);
  }
  
  if (l[0] > 0 || l[1] > 0) {
    fprintf(stderr, "[%s] error: paired end files have differing numbers of reads.\n", __func__);
    exit(1);
  }

  for (i = 0; i < 2; ++i) {
    kseq_destroy(ks[i]);
    gzclose(fp[i]);
  }
  return 0;
}

/* from seqtk.c */
static void cpy_kstr(kstring_t *dst, const kstring_t *src) {
  if (src->l == 0) return;
  if (src->l + 1 > dst->m) {
    dst->m = src->l + 1;
    kroundup32(dst->m);
    dst->s = realloc(dst->s, dst->m);
  }
  dst->l = src->l;
  memcpy(dst->s, src->s, src->l + 1);
}

static void cpy_kseq(kseq_t *dst, const kseq_t *src) {
  cpy_kstr(&dst->name, &src->name);
  cpy_kstr(&dst->comment, &src->comment);
  cpy_kstr(&dst->seq,  &src->seq);
  cpy_kstr(&dst->qual, &src->qual);
}

int split_usage() {
  fputs("\
Usage:    pairs split [options] -1 <out_1.fq> -2 <out_2.fq> -u <out_unpaired.fq> <in.fq>\n\n\
Options:  -n       not strict; won't error if names are not consistent (not recommended)\n\
          -1 FILE  output file name for 1 reads\n\
          -2 FILE  output file name for 2 reads\n\
          -u FILE  output file name for unpaired reads\n\n\
          -m INT   minimum length of reads, equal or shorter will go to unpaired\n\n\
Split interleaved reads into three files: paired-end 1, paired-end 2, and a file for \n\
orphaned reads.\n\n\
Orphaned/unpaired reads are determined by either an empty FASTQ sequence entry, \n\
or a single N. Consistent names mean that the FASTQ/FASTA ID is identical, \n\
except for optional trailing '/1' and '/2', and ignoring comments.\n", stderr);
    return 1;
}

int pairs_split(int argc, char *argv[]) {
  kseq_t **seq, *tmp;
  FILE *fpout[] = {NULL, NULL, NULL};
  gzFile fp;
  int c, l, i, strict=1, min_length=0;
  unsigned total[]={0, 0}, removed[]={0, 0}, both_removed = 0;
  unsigned is_empty[]={0, 0};
  char *p, *tags[] = {"/1", "/2"};
  while ((c = getopt(argc, argv, "1:2:u:n")) >= 0) {
    switch (c) {
    case '1': 
      fpout[0] = fopen(optarg, "w");
      break;
    case '2':
      fpout[1] = fopen(optarg, "w");
      break;
    case 'u':
      fpout[2] = fopen(optarg, "w");
      break;
    case 'm': 
      min_length = atoi(optarg); 
      if (min_length < 1) fprintf(stderr, "[%s] error: minimum length must be >= 1\n", __func__);
      return 1;
      break;
    case 'n': strict = 0; break;
    default: return 1;
    }
  }

  if (optind == argc) return split_usage();

  for (i = 0; i < 3; ++i) {
    if (!fpout[i]) {
      fprintf(stderr, "[%s] error: arguments -1, -2, and -u are required.", __func__);
      return 1;
    }
  }

  fp = (strcmp(argv[optind], "-") == 0) ? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
  
  seq = calloc(2, sizeof(kseq_t*));
  for (i = 0; i < 2; ++i) seq[i] = malloc(sizeof(kseq_t));
  tmp = kseq_init(fp);
  while ((l=kseq_read(tmp)) >= 0) {
    /* always read in chunks of two FASTX entries */
    cpy_kseq(seq[0], tmp);
    l = kseq_read(tmp);
    if (l < 0) break;
    cpy_kseq(seq[1], tmp);

    for (i = 0; i < 2; ++i) {
      /* remove /1 and /2 tags */
      p = strstr(seq[i]->name.s, tags[i]);
      if (p) {
	strncpy(p, "\0", 1);
      }
    }

    if (strcmp(seq[0]->name.s, seq[1]->name.s) != 0) {
      fprintf(stderr, "[%s] warning: interleaved reads names differ '%s' != '%s'\n", __func__, seq[0]->name.s, seq[1]->name.s);
      if (strict) return 1;
    }

    /* deal with unpaired cases (either no seq or single 'N') */
    for (i = 0; i < 2; ++i) {
      is_empty[i] = seq[i]->seq.l <= min_length || strcmp(seq[i]->seq.s, "N") == 0;
      removed[i] += is_empty[i];
      total[i] += 1;
    }
    
    if (!is_empty[0] && !is_empty[1]) {
      for (i = 0; i < 2; i++)      
	printseq(fpout[i], seq[i], seq[i]->seq.l, 0);
    } else if (is_empty[0] && is_empty[1]) {
      both_removed += 1;
      continue;
    } else {
      i = is_empty[0] ? 1 : 0;
      printseq(fpout[2], seq[i], seq[i]->seq.l, 0);
    }
  }
  if (total[0] != total[1]) {
    fprintf(stderr, "[%s] error: mismatched totals of interleaved pairs! %u != %u\n", __func__, total[0], total[1]);
    return 1;
  }
  fprintf(stderr, "totals: %u %u\nremoved: %u %u\n", total[0], total[1], removed[0], removed[1]);
  return 0;
}


int main(int argc, char *argv[]) {
  if (argc == 1) return usage();
  if (strcmp(argv[1], "join") == 0) return pairs_join(argc-1, argv+1);
  if (strcmp(argv[1], "split") == 0) return pairs_split(argc-1, argv+1);
  else {
    fprintf(stderr, "error: unrecognized subcommand '%s'.\n", argv[1]);
    return usage();
  }
  return 0;
}
