#!/usr/bin/env r

PROG_NAME='percent-identity-est.r'
DESCRIPTION = '
	Description: Estimate percent identity between
	two genomes using Bloom filter intersection popcount.
'

headers = c(
	'genome1_popcount',
	'genome2_popcount',
	'genome1_kmers_est',
	'genome2_kmers_est',
	'intersect_popcount',
	'intersect_kmers_est',
	'percent_id_est'
)

#------------------------------------------------------------
# parse command line opts
#------------------------------------------------------------

library('getopt')
spec = matrix(c(
	'bloom_size', 'b', 1, 'double', 'Bloom filter size in bits [required]',
	'header_only', 'H', 0, 'logical', 'print only column headers [false]',
	'help', 'h', 0, 'logical', 'print this help message',
	'kmer_size', 'k', 1, 'double', 'kmer size [required]',
	'no_header', 'N', 0, 'logical', 'omit column headers [false]',
	'num_hash', 'n', 0, 'double', 'number of hash functions [1]',
	'popcount1', '1', 1, 'double', 'number of true bits in Bloom filter 1 [required]',
	'popcount2', '2', 1, 'double', 'number of true bits in Bloom filter 2 [required]',
	'popcount3', '3', 1, 'double',
		 'number of true bits in bitwise AND of Bloom filters 1 and 2 [required]'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# print column headers
if (is.null(opt$no_header)) {
	cat(headers, sep='\t')
	cat('\n');
	if (!is.null(opt$header_only)) {
		q();
	}
}

if (!is.null(opt$help)
	|| is.null(opt$bloom_size)
	|| is.null(opt$kmer_size)
	|| is.null(opt$popcount1)
	|| is.null(opt$popcount2)
	|| is.null(opt$popcount3)) {

	cat(strwrap(DESCRIPTION, simplify=TRUE))
	cat('\n')
	cat(getopt(spec, command=PROG_NAME, usage=TRUE))
	q(status=1)
}

# default value
if (is.null(opt$num_hash)) { opt$num_hash = 1 }

#------------------------------------------------------------
# helper functions
#------------------------------------------------------------

#------------------------------------------------------------
# cardinality_mle
#
# Return the maximum likelihood estimate for the number
# of elements in a bloom filter.
#
# Parameters:
#
#   m: the size of the bloom filters in bits
#   h: the number of hash functions
#   t: the number of true bits in the bloom filter
#
# This is Equation (3) from:
#
# Papapetrou, Odysseas, Wolf Siberski, and Wolfgang Nejdl.
# "Cardinality estimation and dynamic length adaptation for
# bloom filters." Distributed and Parallel Databases 28.2-3
# (2010): 119-156.
#------------------------------------------------------------
cardinality_mle = function(m,h,t) {
	return(log(1-t/m)/(h*log(1-1/m)))
}

#------------------------------------------------------------
# intersect_mle
#
# Return the maximum likelihood estimate for the number
# of elements in the intersection (bitwise AND) of
# two bloom filters.
#
# Parameters:
#
#   m:   the size of the bloom filters in bits (all 3 bloom
#		 filters must be the same size)
#   h:   the number of hash functions (all 3 bloom filters
#		 must use exactly the same set of hash functions)
#   t1:  the number of true bits in bloom filter 1
#   t2:  the number of true bits in bloom filter 2
#   t3:  the nmuber of true in the bitwise AND of bloom
#        filters 1 and 2
#
# This is Equation (6) from:
#
# Papapetrou, Odysseas, Wolf Siberski, and Wolfgang Nejdl.
# "Cardinality estimation and dynamic length adaptation for
# bloom filters." Distributed and Parallel Databases 28.2-3
# (2010): 119-156.
#------------------------------------------------------------
intersect_mle = function(m,h,t1,t2,t3) {
	return((log(m-(t3*m-t1*t2)/(m-t1-t2+t3))-log(m))/(h*log(1-1/m)))
}

#------------------------------------------------------------
# percent_identity_est
#
# Return an estimate for the percent identity between
# two genomes.
#
# Parameters:
#
#     k: kmer size
#     g1: number of kmers in genome 1
#     g2: number of kmers in genome 2
#     o: number of kmers shared by both genomes
#
# Background:
#
# Assuming that the single-base differences between two genomes
# are randomly distributed, then
#
#       O = I^k
#
# where O is percent overlapping kmers, I is percent sequence
# identity, and k is kmer size.
#
# O can be calculated by:
#
#       O = o / min(g1/g2)
#
# Substituting this into the equation above and solving for
# I gives:
#
#       I = (o / min(g1/g2))^(1/k)
#------------------------------------------------------------
percent_identity_est = function(k,g1,g2,o) {
	return((o / min(g1,g2))^(1/k))
}

#------------------------------------------------------------
# do the calculations
#------------------------------------------------------------

b = opt$bloom_size
k = opt$kmer_size
h = opt$num_hash
popcount1 = opt$popcount1
popcount2 = opt$popcount2
popcount3 = opt$popcount3

vals = list()
vals[['genome1_popcount']] = opt$popcount1
vals[['genome2_popcount']] = opt$popcount2
vals[['intersect_popcount']] = opt$popcount3
vals[['genome1_kmers_est']] = cardinality_mle(b, h, popcount1);
vals[['genome2_kmers_est']] = cardinality_mle(b, h, popcount2);
vals[['intersect_kmers_est']] = intersect_mle(b, h, popcount1, popcount2, popcount3);
vals[['percent_id_est']] = percent_identity_est(k,
	vals[['genome1_kmers_est']],
	vals[['genome2_kmers_est']],
	vals[['intersect_kmers_est']]);

# round values for display

vals[['genome1_kmers_est']] = round(vals[['genome1_kmers_est']], 0);
vals[['genome2_kmers_est']] = round(vals[['genome2_kmers_est']], 0);
vals[['intersect_kmers_est']] = round(vals[['intersect_kmers_est']], 0);
vals[['percent_id_est']] = round(vals[['percent_id_est']]*100, 6);

# print column values
first_col = 1
for (h in headers) {
	if (!first_col) { cat('\t') }
	cat(toString(vals[h]))
	first_col = 0
}
cat('\n')
