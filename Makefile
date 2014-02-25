PROGRAM_NAME = seqqs
VERSION = 0.01
CC = gcc
DEBUG ?= 0
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION) -std=gnu99
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0
else 
	CFLAGS += -O3
endif
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
LDFLAGS = -lz
OBJS = seqqs.o
LOBJS = seqqs.o

.PHONY: clean all test

all: seqqs pairs

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

seqqs.o: kseq.h kseq.h

clean: 
	rm -f $(OBJS)
	rm -f $(PROGRAM_NAME)
	rm -f pairs pairs.o

seqqs: $(OBJS)
	$(CC) $(CFLAGS) $? -o $(PROGRAM_NAME) $(LDFLAGS) 

pairs: pairs.o
	$(CC) $(CFLAGS) $? -o pairs $(LDFLAGS) 

lib: libseqqs.so

test:
	(cd tests && python test_pairs.py in-1.fq in-2.fq)

libseqqs.so: CLFAGS += -fpic -D_LIB_ONLY
libseqqs.so: $(LOBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o $@ $^
