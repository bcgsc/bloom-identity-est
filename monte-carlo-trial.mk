#!/usr/bin/make -Rrf

SHELL=/bin/bash -o pipefail 

#------------------------------------------------------------
# constants
#------------------------------------------------------------

# number of hash functions
h:=1

#------------------------------------------------------------
# filenames
#------------------------------------------------------------

genome1?=genome1.fa
genome2?=genome2.fa

#------------------------------------------------------------
# params
#------------------------------------------------------------

# percent sequence identity
pid?=0.95
# bloom filter false positive rate
fpr?=0.05
# genome sequence length
l?=100000
# kmer length
# (make k large enough that the two genomes don't share
# kmers by chance)
k:=$(shell calc-k.r --length $l)
# bloom filter size, in bytes
# (determined by false positive rate)
b:=$(shell calc-bloom-size.r --fpr $(fpr) --num_hash $h --size $l)

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: clean

default: intersection.bloom

clean:
	rm -f $(genome1) $(genome2) *.bloom.gz *.partial

vars:
	@echo 'variable settings:'
	@echo -e '\tpid=$(pid)'
	@echo -e '\tfpr=$(fpr)'
	@echo -e '\tl=$l'
	@echo -e '\tk=$k'
	@echo -e '\tb=$b'

#------------------------------------------------------------
# data generation
#------------------------------------------------------------

$(genome1):
	random-dna $l > $@

$(genome2): $(genome1)
	mutate-bases $(pid) $^ > $@

#------------------------------------------------------------
# Bloom filter rules
#------------------------------------------------------------

# build a bloom filter for a FASTA file
%.bloom: %
	abyss-bloom build -v -k$k -b$b $@ $^

# compute bloom filter intersection 
intersection.bloom: $(genome1).bloom $(genome2).bloom
	abyss-bloom intersect -v -k$k $@ $^
