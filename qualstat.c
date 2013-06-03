#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#define _QS_MAIN

#ifdef USE_SAMTOOLS_LIBS
#include "samtools/khash.h"
#include "samtools/kseq.h"
#else
#include "khash.h"
#include "kseq.h"
#endif

KSEQ_INIT(gzFile, gzread)

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
  unsigned **ntm;
  unsigned **qm;
  unsigned *lm;
  qual_type qt;
} qs_set_t;

qs_set_t *qs_init(qual_type qt) {
  /* 
     Allocate matrices for quality and nucleotides. Rows correspond to
     position in sequence, so growing them is simpler.
  */
  unsigned i;

  qs_set_t *qs = malloc(sizeof(qs_set_t));
  qs->qt = qt;
  qs->m = (size_t) INIT_SEQLEN;
  qs->l = 0;
  qs->ntm = malloc(qs->m*sizeof(unsigned*));
  qs->qm = malloc(qs->m*sizeof(unsigned*));
  qs->lm = calloc(qs->m, sizeof(unsigned));
  for (i = 0; i < qs->m; i++) {
    qs->ntm[i] = calloc(16, sizeof(unsigned));
    qs->qm[i] = calloc(qrng(qs->qt), sizeof(unsigned));
  }
  return qs;
}

void qs_update(qs_set_t *qs, kseq_t *seq) {
  unsigned i, last_m, nt, non_iupac=0;
  char bq;
  if (seq->seq.l > qs->m) {
    last_m = qs->m;
    qs->m = seq->seq.l + 1;
    kroundup32(qs->m);
    /* fprintf(stderr, "[%s] adding rows to matrix (old size: %d; new size: %d)\n", __func__, last_m, qs->m); */
    qs->ntm = realloc(qs->ntm, sizeof(unsigned*)*qs->m);
    qs->qm = realloc(qs->qm, sizeof(unsigned*)*qs->m);
    qs->lm = realloc(qs->lm, sizeof(unsigned)*qs->m);
    for (i = last_m; i < qs->m; i++) {
      qs->ntm[i] = calloc(16, sizeof(unsigned));
      qs->qm[i] = calloc(qrng(qs->qt), sizeof(unsigned));
    }
  }
  if (seq->seq.l > qs->l) qs->l = seq->seq.l;

  if (seq->qual.l && seq->seq.l != seq->qual.l)
    fprintf(stderr, "[%s] quality and sequence lengths differ in sequence '%s'\n", __func__, seq->name.s);

  for (i = 0; i < seq->seq.l; i++) {
    nt = seq_nt17_table[(int) seq->seq.s[i]];
    if (!nt) non_iupac++;
    qs->ntm[i][nt]++;
    if (seq->qual.l) {
      bq = (char) seq->qual.s[i];
      if (bq - qoffset(qs->qt) < qmin(qs->qt) || 
	  bq - qoffset(qs->qt) > qmax(qs->qt))
	fprintf(stderr, "[%s] base quality '%d' out of range (%d <= b <= %d) in sequence '%s'\n", __func__, bq, qmin(qs->qt), qmax(qs->qt), seq->name.s);

      qs->qm[i][bq - qoffset(qs->qt) - qmin(qs->qt)]++;
    }
  }

  if (non_iupac)
    fprintf(stderr, "[%s] warning: %d non-IUPAC characters found in sequence '%s'.\n", __func__, non_iupac, seq->name.s);

}

void qs_qm_fprint(FILE *stream, qs_set_t *qs) {
  unsigned i, j, cnt;
  char bq;

  for (j = 0; j < qrng(qs->qt); j++) {
    bq = j + qoffset(qs->qt) + qmin(qs->qt);
    fprintf(stream, "Q%d", bq);
    if (j < qrng(qs->qt)-1) fputc(' ', stream);
  }
  fputc('\n', stream);

  for (i = 0; i < qs->l; i++) {
    for (j = 0; j < qrng(qs->qt); j++) {
      cnt = qs->qm[i][j];
      fprintf(stream, "%d", cnt);
      if (j < qrng(qs->qt)-1) fputc(' ', stream);
    }
    fputc('\n', stream);
  }
  fputc('\n', stream);
}

void qs_ntm_fprint(FILE *stream, qs_set_t *qs) {
  unsigned i, j, cnt;

  for (j = 0; j < 16; ++j) {
    fputc(seq_nt17_rev_table[j], stream);
    if (j < 15) fputc('\t', stream);
  }
  fputc('\n', stream);

  for (i = 0; i < qs->l; i++) {
    for (j = 0; j < 16; j++) {
      cnt = qs->ntm[i][j];
      fprintf(stream, "%d", cnt);
      if (j < 15) fputc('\t', stream);
    }
    fputc('\n', stream);
  }
  fputc('\n', stream);
}

void qs_lm_print(qs_set_t *qs) {
  
}


void qs_destory(qs_set_t *qs) {
  unsigned i;
  for (i = 0; i < qs->m; i++) {
    free(qs->ntm[i]);
    free(qs->qm[i]);
  }
  free(qs->ntm);
  free(qs->qm);
  free(qs->lm);
  free(qs);
}

#ifdef _QS_MAIN
int main(int argc, char *argv[]) {
  qs_set_t *qs = qs_init(SANGER);
  gzFile fp;
  kseq_t *seq;
  //fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
  fp = gzopen(argv[1], "r");

  seq = kseq_init(fp);
  while (kseq_read(seq) >= 0) {
    qs_update(qs, seq);
  }
  //  qs_ntm_fprint(stdout, qs);
  qs_qm_fprint(stdout, qs);
  kseq_destroy(seq);
  qs_destory(qs);
  gzclose(fp);
}
#endif
