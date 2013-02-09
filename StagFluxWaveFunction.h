#ifndef _STAGFLUXWAVEFUNCTION_H
#define _STAGFLUXWAVEFUNCTION_H

#include "WaveFunction.h"

/*! \brief Staggered flux wave function implementation
 * @param Lx Size of the box along x
 * @param Ly Size of the box along y
 * @param Nfsup number of spin up particles
 * @param Nfsdo number of spin down particles
 * @param phi Staggered flux along a bound
 * @param neel Neel field
 * @param bc_phase Phase shift at the x boundary condition
 * in units of pi. Boundary condition along y is taken periodic
 * (phase shift is 0). A phase shift of pi corresponds to anti-
 * periodic boundary condition.
 */
class StagFluxWaveFunction : public WaveFunction {
    public:
        StagFluxWaveFunction(size_t Lx, size_t Ly,
                             size_t Nfsup, size_t Nfsdo,
                             double phi, double neel,
                             double *bc_phase);
        virtual ~StagFluxWaveFunction();
        virtual std::complex<double>
            matrix_element(size_t f, size_t r, bool up);

    protected:
        size_t *m_qn2fock;
        size_t *m_fock2qn;
        double m_phi;
        double m_neel;
        double m_bc_phase[2];

        /*! utility function to tell whether k is inside
         * or outside of the Magmetic Zone Boundary.*/
        bool inmbz(size_t kx, size_t ky) const;
        /*! utility function to fold back k inside the
         * magnetic Brillouin zone.*/
        void mbzmod(size_t* k) const;
        std::complex<double> delta(double* k) const;
        double omega(double* k) const;
        std::complex<double> Uk(double* k,
                                bool up,
                                size_t band) const;
        std::complex<double> Vk(double* k,
                                bool up,
                                size_t band) const;
        friend std::ostream &
            operator<<(std::ostream& left,
                       const StagFluxWaveFunction & right);
};

#endif//_STAGFLUXWAVEFUNCTION_H
