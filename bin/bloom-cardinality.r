#!/usr/bin/env r

PROG_NAME='bloom-cardinality.r'
DESCRIPTION = '
	Description: Print the maximum likelihood estimate
	for the number of elements in the bloom filter.
'

#------------------------------------------------------------
# parse command line opts
#------------------------------------------------------------

library('getopt')
spec = matrix(c(
	'help', 'h', 0, 'logical', 'print this help message',
	'size', 'b', 1, 'integer', 'bloom filter size in bits [required]',
	'popcount', 'p', 1, 'double', 'number of true bits set to true [required]',
	'num_hash', 'n', 0, 'integer', 'number of hash functions [1]'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (!is.null(opt$help) || is.null(opt$size) || is.null(opt$popcount)) {
	cat(strwrap(DESCRIPTION, simplify=TRUE))
	cat('\n')
	cat(getopt(spec, command=PROG_NAME, usage=TRUE))
	q(status=1)
}

# default value
if (is.null(opt$num_hash)) { opt$num_hash = 1 }

#------------------------------------------------------------
# do the calculation
#------------------------------------------------------------

t = opt$popcount
m = opt$size
h = opt$num_hash

# NOTE: This is Equation (3) from:
# 
# Papapetrou, Odysseas, Wolf Siberski, and Wolfgang Nejdl.
# "Cardinality estimation and dynamic length adaptation for
# bloom filters." Distributed and Parallel Databases 28.2-3
# (2010): 119-156.

cardinality_mle = log(1-t/m)/(h*log(1-1/m))
cat(ceiling(cardinality_mle))
