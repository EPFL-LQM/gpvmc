#include "MatrixQuantity.h"
#include "FileManager.h"
#include "linalg.h"

MatrixQuantity::MatrixQuantity(const Stepper* stepper,
                               FileManager* fm,
                               std::string basename,
                               int M, int N)
    :VectorQuantity(stepper,fm,basename,M*N), m_M(M), m_N(N)
{}

void MatrixQuantity::save() const
{
    std::vector<std::complex<double> > val(m_vals);
    for(size_t i=0;i<val.size();++i) val[i]/=m_meas;
    m_fm->FileStream(m_basename).Write(m_M,2*m_N,(double*)&val[0]);
}

std::string MatrixQuantity::str() const
{
    return linalg::PrintMat(&m_vals[0],m_M,m_N,m_N,4,false);
}

