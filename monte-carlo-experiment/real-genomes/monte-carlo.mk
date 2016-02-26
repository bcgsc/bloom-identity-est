#!/usr/bin/make -Rrf

SHELL=/bin/bash -o pipefail

#------------------------------------------------------------
# helper functions
#------------------------------------------------------------

popcount=$(shell zcat $(1) | abyss-bloom info -v -k$k - 2>&1 | awk -F: '/popcount/ {print $$2}')
eq=$(and $(findstring $(1),$(2)),$(findstring $(2),$(1)))
which=which $1 >/dev/null || (echo "required binary '$1' not found on PATH" && false)
rand=$(shell od -An -N8 -tu8 /dev/random)

#------------------------------------------------------------
# constants
#------------------------------------------------------------

# number of hash functions
h:=1

#------------------------------------------------------------
# params
#------------------------------------------------------------

# check for req'd params
ifndef genome1
$(error missing required arg 'genome1': FASTA file)
endif
ifndef genome1_name
$(error missing required arg 'genome1_name': output file prefix)
endif
ifndef genome2
$(error missing required arg 'genome2': FASTA file)
endif
ifndef genome2_name
$(error missing required arg 'genome2_name': output file prefix)
endif
ifndef n
$(error missing required arg 'n': number of trials)
endif
ifndef out
$(error missing required arg 'out': output text file for results)
endif

# numbero of threads per trial
j?=1
# minimum length for input fasta sequences
s?=500
# minimum seq length threshold for genome 1
s1?=$s
# minimum seq length threshold for genome 2
s2?=$s
# hash seeds for trials
hash_seeds:=$(foreach i,$(shell seq 1 $n),$(call rand))
# print headers for result row?
header?=1
# bloom filter size
b?=10M
# bloom filter size in bits
bits:=$(shell echo '$b*8' | \
	sed -e 's/k/*1024/;s/M/*1024^2/;s/G/*1024^3/' | \
	bc)
# kmer length
# (make k large enough that the two genomes don't share
# kmers by chance)
k?=$(shell calc-k.r --length $(shell fastx-length $(genome1)))
# options to 'abyss-bloom build' commands
bloom_build_opts:=-v -k$k -b$b -j$j

#------------------------------------------------------------
# special rules
#------------------------------------------------------------

.PHONY: vars print_headers run_trials

default: run_trials

clean:
	rm -f *.bloom *.partial

# check for req'd binaries on PATH
check_binaries:
	$(call which, abyss-bloom)
	$(call which, fastx-length)
	$(call which, fasta-minlen)
	$(call which, calc-k.r)
	$(call which, percent-identity-est.r)

vars:
	@echo 'variable settings:'
	@echo -e '\tk=$k'
	@echo -e '\tb=$b'
	@echo -e '\tbits=$(bits)'

#------------------------------------------------------------
# Bloom filter rules
#------------------------------------------------------------

# build bloom filter from FASTA file
$(genome1_name).seed%.bloom.gz: $(genome1)
	smartcat $(genome1) | \
		fasta-minlen $(s1) | \
		abyss-bloom build $(bloom_build_opts) -h$* $(GENOME1_BLOOM_OPT) - - | \
		gzip > $@.partial
	mv $@.partial $@

# build bloom filter from FASTA file
$(genome2_name).seed%.bloom.gz: $(genome2)
	smartcat $(genome2) | \
		fasta-minlen $(s2) | \
		abyss-bloom build $(bloom_build_opts) -h$* $(GENOME2_BLOOM_OPT) - - | \
		gzip > $@.partial
	mv $@.partial $@

# compute bitwise AND of Bloom filters
intersection.seed%.bloom.gz: \
		$(genome1_name).seed%.bloom.gz \
		$(genome2_name).seed%.bloom.gz
	zcat $^ | \
		abyss-bloom intersect -v -k$k - - - | \
		gzip > $@.partial
	mv $@.partial $@

#------------------------------------------------------------
# control rules
#------------------------------------------------------------

trial.seed%: intersection.seed%.bloom.gz
	percent-identity-est.r \
		--no_header \
		--bloom_size $(bits) \
		--kmer_size $k \
		--num_hash $h \
		--popcount1 $(call popcount,$(genome1_name).seed$*.bloom.gz) \
		--popcount2 $(call popcount,$(genome2_name).seed$*.bloom.gz) \
		--popcount3 $(call popcount,intersection.seed$*.bloom.gz) \
		>> $(out)

print_headers:
	percent-identity-est.r --header_only > $(out)

# run a set of Monte Carlo trials with different hash seeds
run_trials: check_binaries print_headers \
	$(foreach seed,$(hash_seeds),trial.seed$(seed))
