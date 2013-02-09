#include "VectorQuantity.h"
#include "FileManager.h"

VectorQuantity::VectorQuantity(const Stepper* stepper, FileManager* fm, std::string basename, int size)
    :Quantity(stepper,fm,basename), m_vals(size,std::complex < double >(0,0))
{}

void VectorQuantity::init()
{
    Quantity::init();
    for(size_t i=0;i<m_vals.size();++i)
        m_vals[i]=0;
}

void VectorQuantity::save() const
{
    std::vector<std::complex<double> > val(m_vals);
    for(size_t i=0;i<val.size();++i) val[i]/=m_meas;
    m_fm->FileStream(m_basename).Write(val.size(),2,(double*)&val[0]);
}

std::string VectorQuantity::str() const
{
    std::ostringstream ostr;
    for(size_t i=0;i<m_vals.size();++i) ostr<<m_vals[i]<<std::endl;
    return ostr.str();
}

