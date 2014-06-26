from scipy.linalg import eigh
from numpy.linalg import matrix_rank
from numpy import dot,conj,amax,argmin,zeros,eye,append,shape,diag,ones
from matplotlib.mlab import find
import copy
import warnings
import code

def argsort(seq):
    return sorted(range(len(seq)),key=seq.__getitem__)

def bunch(instat,Nsamp,indices=False):
    """
    Bunches the integers from instat into Nsamp bunches
    while trying to even out the sum of the bunches integers.
    This is not garanteed to be the best solution. It is just
    a first guess, may probably be optimized further...
    Returns the bunched integers.
    if indices==True, also returns the bunched indices corresponding
    to the original instat list.
    """
    stats=list(instat)
    args=argsort(stats)
    stats.sort()
    bunched=[[] for n in range(Nsamp)]
    argbunched=copy.deepcopy(bunched)
    n=0
    while len(stats)>2*Nsamp-1:
        bunched[n]+=[stats[0],stats[-1]]
        argbunched[n]+=[args[0],args[-1]]
        stats.pop()
        stats.pop(0)
        args.pop()
        args.pop(0)
        n+=1
        if n==Nsamp:
            n=0
    bunchedstat=[sum(b) for b in bunched]
    while len(stats)>0:
        bunched[argmin(bunchedstat)]+=[stats[-1]]
        argbunched[argmin(bunchedstat)]+=[args[-1]]
        bunchedstat[argmin(bunchedstat)]+=stats[-1]
        stats.pop()
        args.pop()
    if not indices:
        return bunched
    else:
        return bunched,argbunched

def geneigh(A,B,tol=1e-12):
    """
    Solves the generalized eigenvalue problem also in the case where A and B share a common
    null-space. The eigenvalues corresponding to the null-space are given a Nan value.
    The null-space is defined with the tolereance tol.
    """
    # first check if there is a null-space issue
    if matrix_rank(B,tol)==shape(B)[0]:
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
