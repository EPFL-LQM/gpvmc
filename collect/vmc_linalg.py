from scipy.linalg import eigh
from numpy.linalg import matrix_rank
from numpy import dot,conj,amax,zeros,eye,append,shape,diag,ones
from matplotlib.mlab import find
import warnings
import code

def geneigh(A,B,tol=1e-12):
    """
    Solves the generalized eigenvalue problem also in the case where A and B share a common
    null-space. The eigenvalues corresponding to the null-space are given a Nan value.
    The null-space is defined with the tolereance tol.
    """
    # first check if there is a null-space issue
    if matrix_rank(B)==shape(B)[0]:
        return eigh(A,B)
    # first diagonalize the overlap matrix B
    Be,Bv=eigh(B)
    # rewrite the A matrix in the B-matrix eigenspace
    At=dot(conj(Bv.T),dot(A,Bv))
    Bt=diag(Be)
    # detect shared null-space. that is given by the first n null eigenvalues of B
    idx=find(Be>tol)
    idx=idx[0]
    # check that the B matrix null-space is shared by A.
    m=amax(abs(At[0:idx,:].flatten()))
    if m>tol:
        warnings.warn('Maximum non-diagonal element in A written in B null-space is bigger than the tolerance \''+str(tol)+'\'.',UserWarning)
    # diagonalize the non-null-space part of the problem
    Et,Vt=eigh(At[idx:,idx:],Bt[idx:,idx:])
    # define Ut, the change of basis in the non-truncated space
    Ut=zeros(shape(A),A.dtype)
    Ut[0:idx,0:idx]=eye(idx)
    Ut[idx:,idx:]=Vt
    U=dot(Bv,Ut)
    E=append(float('NaN')*ones(idx),Et)
    return E,U
