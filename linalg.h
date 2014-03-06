#ifndef _LINALG_H
#define _LINALG_H
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <complex>
#include <string>
#include <iomanip>
#include <sstream>
#include "BigComplex.h"


/*! \brief various linear algebra routines
 *
 * in there:
 * int clapack_zgetrf(const enum CBLAS_ORDER Order, const int M, const int N, void *A, const int lda, int *ipiv); LU decomposition
 * int clapack_zgetri(const enum CBLAS_ORDER Order, const int N, void *A,const int lda, const int *ipiv); inverse from factorization
 */

namespace linalg
{
    /*! Compute the determinant and the inverse of a matrix A. To prevent under- or over-flow, the determinant is returned along with exp, the order of magnitude.
     * returns true if everything went fine or false if A is singular (inverse not returned).
     */
    bool DetInv(const std::complex<double> *A, std::complex<double>* I,
                size_t M, BigComplex& d);
    bool Det(const std::complex<double> *A,size_t M, BigComplex& d);
    /*! \brief Calculate the update of the determinant for several row and columns changes.
     * @param[in] A A square matrix A before rows/columns changes.
     * @param[in] Ai The inverse matrix of A.
     * @param[in] N The size of A.
     * @param[in] C The set(s) of columns to replace in A, size(C)=(N,rc,nc).
     * @param[in] ci The indices of the columns to replace in A, size(ci)=(nc,rc).
     * @param[in] rc The number of columns to replace at ones.
     * @param[in] Nc The number of sets of columns to successively replace in A (see below for in-depth explanation).
     * @param[in] R The set(s) of rows to replace in A, size(R)=(nr,rr,N).
     * @param[in] ri The indices of the rows to replace in A, size(ri)=(nr,rr).
     * @param[in] rr The number of rows to replace at ones.
     * @param[in] Nr The number of sets of rows to successively replace in A.
     * @param[out] dets the updated determinants for each sets of columns/rows, size(dets)=max(nr,nc).
     *
     * An determinant update is calculated for a replacement of rr rows and rc columns at once.
     * As it is computationally cheap to do so, one set of rows (columns) changes can be
     * combined with nc (nr) set of columns (rows). It is required that min(nc,nr)=1.
     * If the values at the crossing of rows and columns do not coincide in R and C,
     * the order in which one inserts the rows and columns matters. If nc>nr, the rows are inserted
     * first and then the columns, so the value in C prevail. If nr>nc, the columns are inserted
     * before the rows and the values in R prevail. If nr=nc=1, the columns are inserted first and the
     * rows after, so the values in R prevail.
     * Complexity: For n=max(nc,nr) and rr>0 and rc>0, the complexity is \f$O(N^2)+n*(O(N)+O((rr+rc)^3))\f$. If either rr=0 or rc=0
     * then the complexity is n*O(N). It is favourable to use this determinant update technique if rr<<N and rc<<N.
     *
     * The implementation uses the rank-(rr+rc) matrix determinant lemma:
     *
     * For a change:
     * \f[A'=A+UV\f]
     * the updated determinant i given by
     * \f[\rm{Det}(A')=\rm{Det}(1+VA^{-1}U))\rm{Det}(A)\f]
     * For rows and columns substitutions, one can define the U and V matrices as:
     * \f[U=\left[\begin{array}{cc}C'&E_R\end{array}\right]\f]
     * \f[V=\left[\begin{array}{c}E_C\\R'\end{array}\right]\f]
     * where \f$E_R\f$ and \f$E_C\f$ are sets of base vectors specifying the rows and columns to be replaced:
     * \f[E_R=\left[\begin{array}{cccc}\hat{\mathbf e}_{r_1}&\hat{\mathbf e}_{r_2}&\dots&\hat{\mathbf e}_{r_{rr}}\end{array}\right]\f]
     * \f[E_C=\left[\begin{array}{c}\hat{\mathbf e}_{c_1}^T\\\hat{\mathbf e}_{c_2}^T\\\dots\\\hat{\mathbf e}_{c_{rc}}^T\end{array}\right]\f]
     * The definition of C' and R' depends whether one inserts rows first and columns after
     * or the opposite.
     *
     * In the nc>nr case:
     * \f[C'=C-AE_C^T\f]
     * and
     * \f[R'=(1-E_C^TE_C)(R-E_R^TA)\f]
     * One needs to calculate the determinant of \f$1+VA^{-1}U\f$:
     * \f[1+VA^{-1}U=\left[\begin{array}{c|c}
     * E_CA^{-1}C & E_CA^{-1}E_R\\\hline
     * RA^{-1}C -RE_C^TE_CA^{-1}C & RA^{-1}E_R-RE_C^TE_CA^{-1}E_R\\
     * -E_R^TC+E_R^TAE_C^TE_CA^{-1}C & +E_R^TAE_C^TE_C^TA^{-1}E_R\end{array}\right]\f]
     * 
     * In the nr>nc case:
     * \f[C'=(1-E_RE_R^T)(C-AE_C^T)\f]
     * \f[R'=R-E_R^TA\f]
     * One needs to calculate the determinant of \f$1+VA^{-1}U\f$:
     * \f[1+VA^{-1}U=\left[\begin{array}{c|c}
     * E_CA^{-1}C-E_CA^{-1}E_RE_R^TC & E_CA^{-1}E_R\\
     * +E_CA^{-1}E_RE_R^TAE_C^T &\\\hline
     * RA^{-1}C-RA^{-1}E_RE_R^TC & RA^{-1}E_R\\
     * -RE_C^T+RA^{-1}E_RE_R^TAE_C^T&\end{array}\right]\f]
     */
    void DetUpdate(const std::complex<double> *A,
                   const std::complex<double> *Ai,
                   size_t N,
                   const std::complex<double> *C,
                   size_t ldc,
                   const size_t * ci,
                   const size_t * rc,
                   const size_t Nc,
                   const std::complex<double> *R,
                   size_t ldr,
                   const size_t * ri,
                   const size_t * rr,
                   const size_t Nr,
                   BigComplex *dets);
     /*! \brief Calculate the update of the inverse matrix for several row and columns changes.
     * @param[in] A A square matrix A before rows/columns changes. A is not updated.
     * @param[out] Ai In: the inverse matrix of A. Out: the updated inverse matrix
     * @param[in] N The size of A.
     * @param[in] C The set(s) of columns to replace in A, size(C)=(N,rc,nc).
     * @param[in] ci The indices of the columns to replace in A, size(ci)=(nc,rc).
     * @param[in] rc The number of columns to replace at ones.
     * @param[in] R The set(s) of rows to replace in A, size(R)=(nr,rr,N).
     * @param[in] ri The indices of the rows to replace in A, size(ri)=(nr,rr).
     * @param[in] rr The number of rows to replace at ones.
     *
     * An determinant update is calculated for a replacement of rr rows and rc columns at once.
     * The values at the crossings between rows and columns must coincide. No attempt is made to
     * prioretize rows or columns.
     * Complexity: For rr>0 and rc>0, the complexity is \f$O(N^2)+O((rr+rc)^3))\f$.
     * It is favourable to use this inverse matrix update technique if rr<<N and rc<<N.
     *
     * The implementation uses the so-called Woodbery formula:
     *
     * For a change:
     * \f[A'=A+UV\f]
     * the updated inverse matrix is
     * \f[A'^{-1}=A^{-1}-A^{-1}U(1+VA^{-1}U)^{-1}VA^{-1}\f]
     * For rows and columns substitutions, one can define the U and V matrices as:
     * \f[U=\left[\begin{array}{cc}C'&E_R\end{array}\right]\f]
     * \f[V=\left[\begin{array}{c}E_C\\R'\end{array}\right]\f]
     * where \f$E_R\f$ and \f$E_C\f$ are sets of base vectors specifying the rows and columns to be replaced:
     * \f[E_R=\left[\begin{array}{cccc}\hat{\mathbf e}_{r_1}&\hat{\mathbf e}_{r_2}&\dots&\hat{\mathbf e}_{r_{rr}}\end{array}\right]\f]
     * \f[E_C=\left[\begin{array}{c}\hat{\mathbf e}_{c_1}^T\\\hat{\mathbf e}_{c_2}^T\\\dots\\\hat{\mathbf e}_{c_{rc}}^T\end{array}\right]\f]
     * The definition of C' and R' depends whether one inserts rows first and columns after
     * or the opposite.
     *
     * Here I use:
     * \f[C'=(1-E_RE_R^T)(C-AE_C^T)\f]
     * \f[R'=R-E_R^TA\f]
     * One needs to calculate the inverse matrix of the capacitance matrix \f$\mathcal{C}=1+VA^{-1}U\f$:
     * \f[\mathcal{C}^{-1}=(1+VA^{-1}U)^{-1}=\left[\begin{array}{c|c}
     * E_CA^{-1}C-E_CA^{-1}E_RE_R^TC & E_CA^{-1}E_R\\
     * +E_CA^{-1}E_RE_R^TAE_C^T &\\\hline
     * RA^{-1}C-RA^{-1}E_RE_R^TC & RA^{-1}E_R\\
     * -RE_C^T+RA^{-1}E_RE_R^TAE_C^T&\end{array}\right]^{-1}\f]
     * And then use the formula taking advantage of the block form:
     * \f[(A')^{-1}=A^{-1}-\left[\begin{array}{cc}A^{-1}C' & A^{-1}E_R\end{array}\right]
     * \left[\begin{array}{c|c}
     * \mathcal{C}^{-1}_{\rm{rc}\times\rm{rc}} & \mathcal{C}^{-1}_{\rm{rc}\times\rm{rr}}\\\hline
     * \mathcal{C}^{-1}_{\rm{rr}\times\rm{rc}} & \mathcal{C}^{-1}_{\rm{rr}\times\rm{rr}}\end{array}\right]
     * \left[\begin{array}{c}E_CA^{-1}\\R'A^{-1}\end{array}\right]\f]
     */
     void InvUpdate(const std::complex<double> *A,
                   std::complex<double> *Ai,
                   size_t N,
                   const std::complex<double> *C,
                   const size_t *ci,
                   const size_t rc,
                   const std::complex<double> *R,
                   const size_t *ri,
                   const size_t rr, BigComplex& det);
    std::string PrintMat(const std::complex<double>* A,
                         int N, int M, int lda=0, int precision=2, bool colwise=true);
    size_t mod(int x, size_t N);//!< Safe modulo for positif AND negative integers. result in [0,N-1].
    int cmod(int x, size_t N);//!< Centered modulo. Result in [-N/2,N/2-1] for even N and in [-N/2,N/2] for odd N.
}

inline size_t linalg::mod(int x, size_t N)
{
    return x>=0 ? size_t(x)%N : size_t(int(N)-abs(x)%N)%N;
}

inline int linalg::cmod(int x, size_t N)
{
    static size_t rx;
    rx=mod(x,N);
    return rx<N/2 ? rx : -int(N)/2 + int(mod(rx,N/2));
}

#endif
