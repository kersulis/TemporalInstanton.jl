"""
    [L,U,Q] = luq(A,tol)
Calculate the following decomposition
A = L |Ubar  0 | Q
      |0     0 |
where Ubar is a square invertible matrix
and matrices L, Q are invertible.
Based on lu decomposition.

Inputs:
* `A`: a sparse matrix
* `tol`: used to separate zero and nonzero values

Note: Julia's lufact always does column pivoting
for sparse matrices.

Ported from Pawel Kowal:
Copyright  (c) Pawel Kowal (2006)
All rights reserved
LREM_SOLVE toolbox is available free for noncommercial academic use only.
pkowal3@sgh.waw.pl
"""
function luq(A,tol)
    A = sparse(A)
    n,m = size(A)

    ###################
    ## Special cases ##
    ###################
    if n == 0 || m == 0
        L = speye(n)
        U = A
        Q = speye(m)
        return L,U,Q
    end

    ######################
    ## LU decomposition ##
    ######################
    F = lufact(A)
    L,U,p,q = F[:L],F[:U],F[:p],F[:q]
    # Q = eye(size(A,2)) # Q'
    Q = diagm(q)'

    p = n - size(L,2)
    LL = [spzeros(n-p,p);speye(p)]
    L = [P'*L P[n-p+1:n,:]']
    U = [U;spzeros(p,m)]

    ##########################
    ## Find zero pivot rows ##
    ##########################
    if size(U,1) == 1 || size(U,2) == 1
        S = U[1,1]
    else
        S = diag(U)
    end
    I = find(abs(S).>tol)
    Jl = setdiff(1:n,I)
    Jq = setdiff(1:m,I)

    Ubar1 = U[I,I]
    Ubar2 = U[Jl,Jq]
    Qbar1 = Q[I,:]
    Lbar1 = L[:,I]

    ##########################################
    ## Eliminate nonzero elements below and ##
    ## on right of invertible block of U,   ##
    ## and update L and Q                   ##
    ##########################################
    if ~isempty(I)
        Utmp = U[I,Jq]
        X = Ubar1'\U[Jl,I]'
        Ubar2 = Ubar2 - X'*Utmp
        Lbar1 = Lbar1 + L[:,Jl]*X'

        X = Ubar1\Utmp
        Qbar1 = Qbar1 + X*Q[Jq,:]
        Utmp = []
        X = []
    end

    ################################################
    ## Find rows and cols with only zero elements ##
    ################################################
    I2 = find(maxabs(Ubar2,2).>tol);
    I5 = find(maxabs(Ubar2,1).>tol);

    I3 = Jl[I2]
    I4 = Jq[I5]
    Jq[I5] = []
    Jl[I2] = []
    U = []

    ########################################
    ## Find part of U which is not in the ##
    ## required form                      ##
    ########################################

    A = Ubar2[I2,I5]

    ####################################
    ## Perform luq decomposition on A ##
    ####################################
    L1,U1,Q1 = luq(A,tol)

    ##################
    ## Update L,U,Q ##
    ##################
    Lbar2 = L[:,I3]*L1
    Qbar2 = Q1*Q[I4,:]
    L = [Lbar1 Lbar2 L[:,Jl]]
    Q = [Qbar1; Qbar2; Q[Jq,:]]

    n1,n2,m2 = length(I), length(I3), length(I4)
    U = [Ubar1 spzeros(n1,m-n1);spzeros(n2,n1) U1 spzeros(n2,m-n1-m2);spzeros(n-n1-n2,m)]

    return L,U,Q
end

"""
    [SpLeft, SpRight] = spspaces(A,opt,tol)
 PURPOSE: finds left and right null and range space of a sparse matrix A
 INPUT:
* A                           a sparse matrix
* opt                         spaces to calculate
                               = 1: left null and range space
                               = 2: right null and range space
                               = 3: both left and right spaces
* tol                         uses the tolerance tol when calculating
                               null subspaces (optional)
  OUTPUT:
* SpLeft                      1x4 cell. SpLeft = {} if opt =2.
       SpLeft{1}               an invertible matrix Q
       SpLeft{2}               indices, I, of rows of the matrix Q that
                               span the left range of the matrix A
       SpLeft{3}               indices, J, of rows of the matrix Q that
                               span the left null space of the matrix A
                               Q(J,:)A = 0
       SpLeft{4}               inverse of the matrix Q
* SpRight                     1x4 cell. SpRight = {} if opt =1.
       SpLeft{1}               an invertible matrix Q
       SpLeft{2}               indices, I, of rows of the matrix Q that
                               span the right range of the matrix A
       SpLeft{3}               indices, J, of rows of the matrix Q that
                               span the right null space of the matrix A
                               AQ(:,J) = 0
     SpLeft{4}               inverse of the matrix Q
  COMMENTS:
   uses luq routine, that finds matrices L, U, Q such that
       A = L | U 0 | Q
             | 0 0 |
   where L, Q, U are invertible matrices, U is upper triangular. This
   decomposition is calculated using lu decomposition.
   This routine is fast, but can deliver inaccurate null and range
   spaces if zero and nonzero singular values of the matrix A are not
   well separated.
  WARNING:
   right null and range space may be very inaccurate
Copyright  (c) Pawel Kowal (2006)
All rights reserved
LREM_SOLVE toolbox is available free for noncommercial academic use only.
pkowal3@sgh.waw.pl
"""
function spspaces(A,tol=max(max(size(A))*norm(A,1)*eps,100*eps))
    L,U,Q = luq(A,tol)

    if !isempty(Q)
        QQ = inv(Q)
    else
        QQ = Q
    end
    S = maxabs(U,1)
    I = find(S.>tol)
    if !isempty(S)
        J = find(S.<=tol)
    else
        J = (1:size(S,2))'
    end
    return QQ[J,:]
end
