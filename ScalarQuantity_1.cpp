#include "ScalarQuantity_1.h"
#include "FileManager.h"

void ScalarQuantity_1::init()
{
    Quantity_1::init();
    m_vals=0;
}

void ScalarQuantity_1::save() const
{
    std::complex<double> val(m_vals/double(m_meas));
    m_fm->FileStream(m_basename).Write(1,2,(double*)&val);
}

std::string ScalarQuantity_1::str() const
{
    std::ostringstream ostr;
    ostr<<m_vals/double(m_meas);
    return ostr.str();
}

