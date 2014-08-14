#ifndef _SFPNPHXWAVEFUNCTION_H
#define _SFPNPHXWAVEFUNCTION_H

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
class SFpNpHxWaveFunction : public WaveFunction {
    public:
        SFpNpHxWaveFunction(FileManager* fm,
                            size_t Lx, size_t Ly,
                            size_t Nby,
                            double phi, double neel, double hx,
                            const std::vector<double>& bc_phase);
        virtual ~SFpNpHxWaveFunction();

    protected:
        size_t m_Lx;
        size_t m_Ly;
        size_t *m_qn2fock;
        size_t *m_fock2qn;
        double m_phi;
        double m_neel;
        double m_hx;
        std::vector<double> m_bc_phase;

        virtual std::complex<double>
            matrix_element(size_t f, size_t r, size_t up);
        virtual void quantum_numbers(const size_t& f, const size_t& fl, std::map<std::string,size_t>& qn);
        /*! utility function to tell whether k is inside
         * or outside of the Magmetic Zone Boundary.*/
        bool inmbz(size_t kx, size_t ky) const;
        /*! utility function to fold back k inside the
         * magnetic Brillouin zone.*/
        void mbzmod(size_t* k) const;
        int sign( double x) const;
        std::complex<double> delta(double* k) const;
        double omega(double* k, size_t band) const;
        std::complex<double> Uk(double* k,
                                bool up,
                                size_t band) const;
        std::complex<double> Vk(double* k,
                                bool up,
                                size_t band) const;
        friend std::ostream &
            operator<<(std::ostream& left,
                       const SFpNpHxWaveFunction & right);
};

#endif//_SFPNPHXWAVEFUNCTION_H
