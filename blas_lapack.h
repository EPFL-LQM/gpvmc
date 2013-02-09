#ifndef _BLAS_LAPACK_H
#define _BLAS_LAPACK_H

#ifdef PROFILE
#include "Timer.h"
#endif

extern "C"{
    void zgetrf_(const int* M, const int* N, double* A, const int* LDA,
                 int* IPIV, int* info);
    void zgetri_(const int* M, double* A, const int* LDA, const int* IPIV,
                 double* work, const int* lwork, int* info);
    void zgemv_(const char* TRANS,
                const int* M,
                const int* N,
                const double* alpha,
                const double* A,
                const int* lda,
                const double* x,
                const int* incx,
                const double* beta,
                double* y,
                const int* incy);
    void zgemm_(const char* TransA,
               const char* TransB,
               const int* M,
               const int* N,
               const int* K,
               const double* alpha,
               const double* A,
               const int* lda,
               const double* B,
               const int* ldb,
               const double* beta,
               double* C,
               const int* ldc);
    void zdotu_sub_(const int* N,
                    const double* x,
                    const int* incx,
                    const double* y,
                    const int* incy,
                    double* dotu);
    void zgeru_(const int* M,
                const int* N,
                const double* alpha,
                const double* x,
                const int* incx,
                const double* y,
                const int* incy,
                double* A,
                const int* lda);
}

namespace linalg {
    //enum Transpose{ NoTrans, Trans, ConjTrans};
    void zgetrf(const int M, const int N, std::complex<double>* A,
                const int LDA, int* IPIV, int* info);
    void zgetri(const int M, std::complex<double>* A, const int LDA,
                int* IPIV, std::complex<double>* work,
                const int lwork, int* info);
    void zgemv(const char Trans,
               const int M, const int N,
               const std::complex<double>* alpha,
               const std::complex<double>* A,
               const int lda,
               const std::complex<double>* x,
               const int incx,
               const std::complex<double>* beta,
               std::complex<double>* y,
               const int incy);

    void zgemm(const char TransA,
               const char TransB,
               const int M,
               const int N,
               const int K,
               const std::complex<double>* alpha,
               const std::complex<double>* A,
               const int lda,
               const std::complex<double>* B,
               const int ldb,
               const std::complex<double>* beta,
               std::complex<double>* C,
               const int ldc);

    void zdotu_sub(const int N,
                   const std::complex<double>* x,
                   const int incx,
                   const std::complex<double>* y,
                   const int incy,
                   std::complex<double>* dotu);

    void zgeru(const int M,
               const int N,
               const std::complex<double>* alpha,
               const std::complex<double>* x,
               const int incx,
               const std::complex<double>* y,
               const int incy,
               std::complex<double>* A,
               const int lda);
}

inline void linalg::zgetrf(const int M, const int N, std::complex<double>* A,
                           const int LDA, int* IPIV, int* info)
{
    zgetrf_(&M,&N,(double*)A,&LDA,IPIV,info);
}
inline void linalg::zgetri(const int M, std::complex<double>* A,
                           const int LDA, int* IPIV,
                           std::complex<double>* work, const int lwork,
                           int* info)
{
    zgetri_(&M,(double*)A,&LDA,IPIV,(double*)work,&lwork,info);
}

inline void linalg::zgemv(const char Trans,
               const int M, const int N,
               const std::complex<double>* alpha,
               const std::complex<double>* A,
               const int lda,
               const std::complex<double>* x,
               const int incx,
               const std::complex<double>* beta,
               std::complex<double>* y,
               const int incy)
{
    zgemv_(&Trans,&M,&N,(const double*)alpha,
           (const double*)A,&lda,
           (const double*)x,&incx,
           (const double*)beta,
           (double*)y,&incy);
}

inline void linalg::zgemm(const char TransA,
               const char TransB,
               const int M,
               const int N,
               const int K,
               const std::complex<double>* alpha,
               const std::complex<double>* A,
               const int lda,
               const std::complex<double>* B,
               const int ldb,
               const std::complex<double>* beta,
               std::complex<double>* C,
               const int ldc)
{
    zgemm_(&TransA,&TransB,&M,&N,&K,
           (const double*)alpha,
           (const double*)A,
           &lda,
           (const double*)B,
           &ldb,
           (const double*)beta,
           (double*)C,
           &ldc);
}

inline void linalg::zdotu_sub(const int N,
                   const std::complex<double>* x,
                   const int incx,
                   const std::complex<double>* y,
                   const int incy,
                   std::complex<double>* dotu)
{
    zdotu_sub_(&N,(const double*)x,&incx,(const double*)y,&incy,(double*)dotu);
}

inline void linalg::zgeru(const int M,
               const int N,
               const std::complex<double>* alpha,
               const std::complex<double>* x,
               const int incx,
               const std::complex<double>* y,
               const int incy,
               std::complex<double>* A,
               const int lda)
{
    zgeru_(&M,&N,
           (const double*)alpha,
           (const double*)x,
           &incx,
           (const double*)y,
           &incy,
           (double*)A,
           &lda);
}

#endif//_BLAS_LAPACK_H
