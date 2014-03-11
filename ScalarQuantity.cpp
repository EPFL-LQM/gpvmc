#include "ScalarQuantity.h"
#include "FileManager.h"

void ScalarQuantity::init()
{
    Quantity::init();
    m_vals=0;
}

void ScalarQuantity::save() const
{
    std::complex<double> val(m_vals/double(m_meas));
    m_fm->FileStream(m_basename).Write(1,2,(double*)&val);
}

std::string ScalarQuantity::str() const
{
    std::ostringstream ostr;
    ostr<<m_vals/double(m_meas);
    return ostr.str();
}

