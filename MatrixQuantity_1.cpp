#include "MatrixQuantity_1.h"
#include "FileManager.h"
#include "linalg.h"

MatrixQuantity_1::MatrixQuantity_1(const Stepper_1* stepper,
                               FileManager* fm,
                               std::string basename,
                               int M, int N)
    :VectorQuantity_1(stepper,fm,basename,M*N), m_M(M), m_N(N)
{}

void MatrixQuantity_1::save() const
{
    std::vector<std::complex<double> > val(m_vals);
    for(size_t i=0;i<val.size();++i) val[i]/=m_meas;
    m_fm->FileStream(m_basename).Write(m_M,2*m_N,(double*)&val[0]);
}

std::string MatrixQuantity_1::str() const
{
    return linalg::PrintMat(&m_vals[0],m_M,m_N,m_N,4,false);
}

