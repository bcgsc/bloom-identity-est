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
bytes:=$(shell calc-bloom-size.r --fpr $(fpr) --num_hash $h --size $l)
bits:=$(shell echo '8*$(bytes)' | bc)

#------------------------------------------------------------
# helper functions
#------------------------------------------------------------

popcount=$(shell abyss-bloom info -v -k$k $(1) 2>&1 | awk -F: '/popcount/ {print $$2}')

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: clean vars print_results

default: print_results

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
	abyss-bloom build -v -k$k -b$(bytes) $@ $^

# compute bloom filter intersection 
intersection.bloom: $(genome1).bloom $(genome2).bloom
	abyss-bloom intersect -v -k$k $@ $^

#------------------------------------------------------------
# Print results to STDOUT
#------------------------------------------------------------

print_results: intersection.bloom
	@# estimated number of kmers in genome 1
	@bloom-cardinality.r \
		--size $(bits) \
		--num_hash $h \
		--popcount $(call popcount,$(genome1).bloom)
	@echo -ne '\t'
	@# estimated number of kmers in genome 2
	@bloom-cardinality.r \
		--size $(bits) \
		--num_hash $h \
		--popcount $(call popcount,$(genome2).bloom)
	@echo -ne '\t'
	@# number of bits bloom intersection
	@# estimated number of kmers in both genomes 1 and 2
	@echo -n $(call popcount,intersection.bloom)
	@echo -ne '\t'
	@intersect-cardinality.r \
		--size $(bits) \
		--num_hash $h \
		--popcount1 $(call popcount,$(genome1).bloom) \
		--popcount2 $(call popcount,$(genome2).bloom) \
		--overlap_count $(call popcount,intersection.bloom)
