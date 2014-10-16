#ifndef _SFPNXPHZWAVEFUNCTION_H
#define _SFPNXPHZWAVEFUNCTION_H

#include "WaveFunction.h"

/*! \brief Staggered flux plus x-Neel plus z-field wave function implementation
 */
class SFpNxpHzWaveFunction : public WaveFunction {
    public:
        SFpNxpHzWaveFunction(FileManager* fm,
                             size_t Lx, size_t Ly,
                             size_t Nby,
                             double phi, double neel,
                             double neel_exp, double hz,
                             const std::vector<double>& bc_phase);
        virtual ~SFpNxpHzWaveFunction();

    protected:
        size_t m_Lx;
        size_t m_Ly;
        size_t *m_qn2fock;
        size_t *m_fock2qn;
        double m_phi;
        double m_nx;
        double m_ne;
        double m_hz;
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
        double neelk(double* k) const;
        double omega(double* k, size_t band) const;
        std::complex<double> Uk(double* k,
                                size_t band) const;
        std::complex<double> Vk(double* k,
                                size_t band) const;
        std::complex<double> Xk(double* k,
                                size_t band) const;
        std::complex<double> Yk(double* k,
                                size_t band) const;
        friend std::ostream &
            operator<<(std::ostream& left,
                       const SFpNxpHzWaveFunction & right);
};

#endif//_SFPNPHXWAVEFUNCTION_H
