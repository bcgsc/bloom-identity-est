#!/usr/bin/env r

PROG_NAME='intersect-cardinality.r'
DESCRIPTION = '
	Description: Print the maximum likelihood estimate
	for the number of elements represented by overlapping
	bits of two Bloom filters.
'

#------------------------------------------------------------
# parse command line opts
#------------------------------------------------------------

library('getopt')
spec = matrix(c(
	'help', 'h', 0, 'logical', 'print this help message',
	'size', 'b', 1, 'double', 'Bloom filter size in bits [required]',
	'popcount1', '1', 1, 'double', 'number true bits in Bloom filter 1 [required]',
	'popcount2', '2', 1, 'double', 'number of true bits in Bloom filter 2 [required]',
	'overlap_count', 'o', 1, 'double', 'number of true bits in the bitwise AND of Bloom filters 1 and 2 [required]',
	'num_hash', 'n', 0, 'integer', 'number of hash functions [1]'
), byrow=TRUE, ncol=5)

opt = getopt(spec)

if (!is.null(opt$help) || is.null(opt$size)
	|| is.null(opt$popcount1)
	|| is.null(opt$popcount2)
	|| is.null(opt$overlap_count)) {

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

m = opt$size
t1 = opt$popcount1
t2 = opt$popcount2
t12 = opt$overlap_count
h = opt$num_hash

# NOTE: This is Equation (6) from:
# 
# Papapetrou, Odysseas, Wolf Siberski, and Wolfgang Nejdl.
# "Cardinality estimation and dynamic length adaptation for
# bloom filters." Distributed and Parallel Databases 28.2-3
# (2010): 119-156.

cardinality_mle = (log(m-(t12*m-t1*t2)/(m-t1-t2+t12))-log(m))/(h*log(1-1/m))
cat(ceiling(cardinality_mle))
