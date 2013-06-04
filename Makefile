PROGRAM_NAME = seqqs
VERSION = 0.01
CC = clang
DEBUG ?= 0
CFLAGS = -Wall -pedantic -DVERSION=$(VERSION) -std=c11
ifeq ($(DEBUG), 1)
	CFLAGS += -g -O0
else 
	CFLAGS += -O3
endif
ARCHIVE = $(PROGRAM_NAME)_$(VERSION)
LDFLAGS = -lz
OBJS = seqqs.o

.PHONY: clean all

default: all

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

seqqs.o: kseq.h kseq.h

clean: 
	rm -f $(OBJS)
	rm -f $(PROGRAM_NAME)

all: $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $? -o $(PROGRAM_NAME)
