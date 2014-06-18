#ifndef _SLATERDETERMINANT_H
#define _SLATERDETERMINANT_H
#include "BigDouble.h"
#include "BigComplex.h"
#include "defs.h"
#include "Amplitude.h"

class WaveFunction;

/*! \brief scalar product \f$\langle \{R_i,\sigma_i\}|\Psi\rangle\f$
 *
 * Class to calculate the scalar product between the real space spin state
 * \f$|{R_i,\sigma_i}\rangle\f$ (represented by the class SpinState) and the
 * trial wave function \f$|\Psi\rangle\f$ (represented by classes inheriting
 * from the class WaveFunction). The class uses the determinant and inverse
 * matrix update technique.
 * For a column vector \f$u\f$ update (hop in the trial wavefunction space):
 * \f[\frac{det(A')}{det(A)}=\frac{det(A+(u-Ae_k)e_k^T)}{det(A)}
 * =e_k^TA^{-1}u\f]
 * and
 * \f[(A')^{-1}=(A+(u-Ae_k)e_k^T)^{-1}=A^{-1}-\frac{(A^{-1}u-e_k)e_k^TA^{-1}}
 * {e_k^TA^{-1}u}\f]
 * For a row vector \f$v^T\f$ update (hop in the real space):
 * \f[\frac{det(A')}{det(A)}=\frac{det(A+e_r(v^T-e_r^TA))}{det(A)}
 * =v^TA^{-1}e_r\f]
 * and
 * \f[(A')^{-1}=(A+e_r(v^T-e_r^TA))^{-1}=
 * A^{-1}-\frac{A^{-1}e_r(v^TA^{-1}-e_r^T)}
 * {v^TA^{-1}e_r}.\f]
 * For a successive \f$r\f$ row and \f$k\f$ column update by repsectively
 * \f$v^T\f$ and \f$u\f$ vectors:
 * \f[\frac{det(A)}{det(A')}=\left|\begin{array}{cc}
 * e_k^TA^{-1}u+A^{-1}_{kr}A_{rk} & A^{-1}_{kr}\\
 * v^TA^{-1}u+A_{rk}v^TA^{-1}e_r-u_r+v_k(e_k^TA^{-1}u+A_{rk}A^{-1}_{kr}) &
 * v^TA^{-1}e_r-v_kA^{-1}_{kr}\end{array}\right|,\f]
 * and for a successive \f$k\f$ column and \f$r\f$ row update by repsectively
 * \f$u\f$ and \f$v^T\f$ vectors:
 * \f[\frac{det(A)}{det(A')}=\left|\begin{array}{cc}
 * e_k^TA^{-1}u-u_rA^{-1}_{kr} & A^{-1}_{kr}\\
 * v^TA^{-1}u-v_k+A_{rk}e_k^TA^{-1}u-u_r(v^TA^{-1}e_r+A_{rk}A^{-1}_{kr})
 * & v^TA^{-1}e_r+A_{rk}A^{-1}_{kr}\end{array}\right|.\f]
 */

class SlaterDeterminant: public Amplitude {
    public:
        SlaterDeterminant(const LatticeState* sp, const WaveFunction* wav);//!< Constructor
        ~SlaterDeterminant();//!< Destructor
        /*! Initialize or reinitialize the matrices
         * Calculate determinant and inverse matrices.*/
        virtual void Init();
        /*! gives the new determinants after the row update
         * specified by r(...) for all subsequent column updates
         * specified by k(...)*/
        void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                        const std::vector<std::vector<hop_path_t> >& khop,
                        vector<BigComplex>& qs) const;

        void Update(const std::vector<hop_path_t>& rhop,
                    const std::vector<hop_path_t>& khop);

        //! returns \f$\langle \{R_i,\sigma\}|\Psi\rangle\f$.
        BigComplex Amp() const;
    private:
        const WaveFunction* m_wav;
        std::vector<std::vector<std::complex<double> > > m_mat;
        std::vector<std::vector<std::complex<double> > > m_mati;
        BigComplex m_amp;
        bool m_amp_ok;
        uint_vec_t m_N;
        size_t m_Nfl;
};

#endif//_SLATERDETERMINANT_H
