#ifndef _SPINDENSITY_H
#define _SPINDENSITY_H

#include <iomanip>
#include "MatrixQuantity.h"

/*! Spin density wave single mode approximation quantity.
 *
 * The formula is:
 * \f[E_q=\frac{\langle\Psi_{GS}|[S_q^-,[\mathcal{H},S_q^+]]|\Psi_{GS}\rangle}{\langle\Psi_{GS}|S_q^-S_q^+|\Psi_{GS}\rangle}\f]
 */

class SpinDensity: public MatrixQuantity
{
    public:
        SpinDensity(const Stepper* stepper,
                    FileManager* fm,
                    std::vector<std::vector<double> > qs,
                    size_t L);
        virtual ~SpinDensity() {}
        virtual void measure();
    private:
        std::vector<std::vector<double> > m_qs;
        std::vector<std::vector<std::complex<double> > > m_ph; 
        //std::vector<std::vector<std::complex<double> > > m_xphases;
        //std::vector<std::vector<std::complex<double> > > m_yphases;
};

#endif//_SPINDENSITY_H
