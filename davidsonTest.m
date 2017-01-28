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
dim = 5000;

# define number of eigenvalues to find
eigVals = 2;

# generate random test matrix
ATest = rand(dim,dim);
ATest = ATest'+ATest; # symmetrize
ATest = ATest + diag(eye(dim,1)) - 5000*eye(dim);

# form initial search subspace
vTest = orth(rand(dim,eigVals));

# form subspace based on A
#vTest2 = zeros(dim,eigVals);
#[aD, aIndex] = sort(diag(ATest))
#for ( vIndex = 1:columns(vTest) )
#  vTest2(aIndex(vIndex)) = 1;
#endfor

#ATest
#vTest2

# define convergence criteria
cutTest = 1e-4; # cutoff
maxTest = 300; # max iteration

printf("Mine:\n");
davidson(ATest, vTest, cutTest, maxTest);

expected = sort(eigs(ATest));
printf("Expected:   ..."),;
for( i = 1:eigVals )
  printf("%16.8f", expected(i));
endfor
printf("\n");

