clc;
clear all;
close all;

# ******************************************
# Davidson Test (davidsonTest.m)
#  Last Edited: Jan 26, 2017 (Shannon Houck)
#
# Tests the Davidson algorithm (davidson.m)
#  and compares results to expected values.
#
# ******************************************

# define dimensions of test matrix
dim = 6;

# define number of eigenvalues to find
eigVals = 3;

# generate random test matrix
ATest = rand(dim,dim);
ATest = ATest'+ATest; # symmetrize
ATest = ATest + diag(eye(dim,1));

# form initial search subspace
vTest = orth(rand(dim,eigVals));

# define convergence criteria
cutTest = 1e-4; # cutoff
maxTest = 30; # max iteration

printf("Mine:\n");
davidson(ATest, vTest, cutTest, maxTest);

printf("Expected:\n"),;
eigs(ATest, dim),;

