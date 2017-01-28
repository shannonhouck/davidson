# ******************************************
# Davidson Solver (davidson.m)
#  Last Edited: Jan 23, 2017 (Shannon Houck)
#
# Defines a function that solves for the k
#  lowest eigenvalues and eigenvectors of a 
#  given matrix A, via the Davidson algorithm.
#
# Params:
#   A - matrix to solve for
#   rInit - initial set of guess vectors
#   threshold - cutoff for convergence
#   maxIter - maximum number of iterations
# ******************************************

function davidson( A, vInit, cutoff, maxIter )

# generate vSpace (search subspace)
# based on vInit (input vectors)
vSpace = vInit;

# size at which to collapse search subspace
collapseSize = 100;

# number of eigenvalues to solve for
k = columns(vInit);

# iterations completed
j = 0;

sig = [];

# index of last sigma added
lastSig = 1;

while ( j<maxIter )

  # form sigma vectors
  for ( i = vSpace(:, lastSig:end ) )
    sig = [sig, A*i];
  endfor
  
  # form subspace matrix
  Av = vSpace'*sig;
  # solve for k lowest eigenvalues
  [eVects, eVals] = eigs(Av, columns(Av));
  # eIndex tracks where the lowest values are
  # so we can extract the appropriate eigenvectors
  [eVals, eIndex] = sort(diag(eVals));
  # discard high values
  eVals = eVals(1:k);
  eIndex = eIndex(1:k);

  r = [];
  # compute residuals
  for ( i = 1:k )
    rNew = sig*eVects(:,eIndex(i)) - eVals(i)*vSpace*eVects(:,eIndex(i));
    r = [r, rNew];
  endfor

  printf("Iteration: %4i", j)
  for ( i = eVals )
    printf("%16.8f", i)
  endfor
  for ( i = 1:k )
    printf("%16.8f", norm(r(:,i)))
  endfor
  printf("\n")

  # check residuals for convergence (break)
  converged = true;
  for ( i = 1:k )
    norm(r(:,i));
    if( norm(r(:,i))>cutoff )
      converged = false;
    endif
  endfor

  if ( converged == true )
    printf("Done! \n");
    return;
  endif

 #collapse subspace if necessary (if)
  if( columns(vSpace) > collapseSize )
    lastSig = 1;
    vSpaceNew = [];
    for ( l = 1:k )
      newVect = zeros(rows(A),1);
      for ( i = 1:columns(vSpace) )
        newVect = newVect + eVects(i,eIndex(l))*vSpace(:,i);
      endfor
      # orthonormalize
      newVect = newVect/norm(newVect);
      vSpaceNew = [vSpaceNew, newVect];
    endfor
    vSpace = vSpaceNew;
    sig = [];

  # else, apply preconditioner to residuals (elif)
  else
    # apply preconditioner
    D = diag(diag(A));
    for ( i = 1:k )
      sNew = inv(D-eVals(i)*eye(length(D))) \ r(:,i);
      if ( norm(sNew) > 0 ) 
        # orthogonalize
        h = vSpace'*sNew;
        sNew = sNew - vSpace*h;
        # normalize
        sNew = sNew/norm(sNew);
        vSpace = [vSpace, sNew];
        lastSig++;
      endif
      # add to subspace
    endfor
    if ( columns(vSpace) > columns(A) )
      printf("...\nError: Make sure your inputs are reasonable!\n\n");
      return;
    endif
  endif

  # increase j and loop again
  j++;

endwhile

endfunction
