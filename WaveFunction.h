#ifndef _WAVEFUNCTION_H
#define _WAVEFUNCTION_H

#include <limits>
#include <complex>
#include <omp.h>
#include <iomanip>
#include <vector>
#include <utility>

#ifndef _HOP_PATH_T
typedef std::vector<std::pair<size_t,size_t> > hop_path_t;
#endif

/*! \brief Base class for trial wave functions.
 *
 */

class WaveFunction {

    /**************
     * protected: *
     **************/
    protected:
        /*! \brief spin up Fermi sea.
         * Fock space state index of the up spins particles.
         * m_fsup[p] is the Fock space state index of the
         * particle indexed by p.*/
        size_t* m_fsup;
        /*! \brief spin down Fermi sea.
         * See m_fsup for detailed description.*/
        size_t* m_fsdo;
        /*! \brief Many particle state written in Fock space
         * with indiced particles, up spins part. m_fockup[st]
         * is the the index of the particle in Fock state index
         * st or m_Nfsup if no particle occupies this state.*/
        size_t *m_fockup;
        /*! \brief Many particle state written in Fock space
         * with indiced particles, down spins part*/
        size_t *m_fockdo;
        /*! \brief cache of all spin up matrix elements.
         *
         * Single particle matrix elements are given by
         * \f$\langle k,\uparrow | R,\uparrow\rangle\f$
         * for k a set of quantum numbers completely describing
         * a single particle state, spin excepted. The matrix
         * must be contiguous in R, thus is indiced as
         * m_cacheup[k*Lx+Ly+Rx*Ly+Ry]*/
        std::complex<double>* m_cacheup;
        /*! \brief See description of m_cacheup.*/
        std::complex<double>* m_cachedo;
        /*! \brief Excitation subspace definition.
         * For each excited state, a hopping sequence
         * is defined in order to reach all the other
         * states. m_excup(do)[i][j] indicates the hopping
         * sequence of spins up (down) to go from state
         * i to state j.
         * This matrix of hopping sequence is build when the
         * wavefunction is created from the definition
         * of the excitation subspace (see constructor)
         */
        std::vector<std::vector<hop_path_t> > m_excup;
        /*! \brief See documentation of m_excup.
         */
        std::vector<std::vector<hop_path_t> > m_excdo;
        /*! \brief Current excitation index.*/
        size_t m_c_exc;
        size_t m_Nfsup;//!< size of up Fermi sea
        size_t m_Nfsdo;//!< size of down Fermi sea
        size_t m_Lx; //!< x size of real space lattice
        size_t m_Ly; //!< y size of real space lattice
        int m_sign; //!< fermion sign of the wavefunction

        /*! Initialize m_cacheup and m_cachedo with the derived
         * class matrix_element method. This function must be
         * explicitely called at the end of the derived class
         * constructor.*/
        void init_matrices();
        /*! Derived classes must define what are the matrix
         * elements
         * \f$\langle k,\downarrow | R,\downarrow\rangle\f$.
         * where k is the fock index associated with a
         * complete set of quantum numbers specifying
         * a single particle state.*/
        virtual std::complex<double>
            matrix_element(size_t f, size_t r,bool up)=0;

        hop_path_t get_hop_path(
                const size_t* src,
                const size_t* dest,
                size_t Nsf) const;

        void get_hop_state(hop_path_t hop,
                           const size_t* src,
                           size_t Nfs,
                           std::vector<size_t>& dest) const;

    /**************
     * public:    *
     **************/
    public:
        /*! Base constructor*/
        WaveFunction(size_t Lx, size_t Ly,
                     size_t Nfsup, size_t Nfsdo);

        virtual ~WaveFunction(); //!< Base destructor

        std::complex<double> MatEl(const size_t& f,
                                   const size_t& r,
                                   bool up) const
        {
            if(up) return m_cacheup[f*m_Lx*m_Ly+r];
            else   return m_cachedo[f*m_Lx*m_Ly+r];
        }
        size_t GetNup() const {return m_Nfsup;}
        size_t GetNdo() const {return m_Nfsdo;}
        const size_t* GetUp() const {return m_fsup;}
        const size_t* GetDo() const {return m_fsdo;}
        int GetSign() const {return m_sign;}
        size_t GetNExc() const {return m_excup.size();}

        /*! Fill up states after initialization.
         * **st are arrays
         * of Nidx*Nfs** elements specifying all the single
         * particle states.*/
        void fill(size_t *upst, size_t *dost);

        /*! Add a state to the wavefunction, defining
         * the hop path from and to the other already defined
         * states. fup and fdo are Fock states where m_Nfsup(do)
         * stands for empty in fup(do) state and any other
         * integer for occupied. The size of fup and fdo must
         * be m_Lx*m_Ly.
         */
        void add_state(const size_t* fup,
                       const size_t* fdo);

        void GetHops(std::vector<hop_path_t>& hopup,
                     std::vector<hop_path_t>& hopdo) const;

        void GetHop(size_t k,
                    hop_path_t& hopup,
                    hop_path_t& hopdo) const;

        /*! do \f$|\Psi\rangle\rightarrow\
         * c^\dagger_{k_\uparrow'}\
         * c_{k_\uparrow}c^\dagger_{k_\downarrow'}\
         * c_{k_\downarrow}\
         * |\Psi\rangle\f$,
         * update matrix elements and determinant.*/
        void hop(size_t khop);

        /* \brief Calculate the hopping fermionic sign */
        virtual int hop_sign(const hop_path_t& hopup,
                             const hop_path_t& hopdo) const;

        friend std::ostream &
            operator<<(std::ostream& left,
                       const WaveFunction & right);
};

#endif//_WAVEFUNCTION_H
