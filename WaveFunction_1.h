#ifndef _WAVEFUNCTION_1_H
#define _WAVEFUNCTION_1_H

#include <limits>
#include <complex>
#include <omp.h>
#include <iomanip>
#include <vector>
#include <utility>

#ifndef _HOP_PATH_T
typedef std::vector<std::pair<size_t,size_t> > hop_path_t;
typedef std::pair<size_t,size_t> hop_t;
#endif

class FileManager;

/*! \brief Base class for trial wave functions.
 *
 */

class WaveFunction_1 {

    /**************
     * protected: *
     **************/
    protected:
        /*! \brief Enumerated particles Fock states.
         * Fock state index in which the enumerated particles
         * are in. Organized in a block form if the spin (flavour)
         * is a good quantum number. m_fs[f,p] is the Fock state index
         * of the p'th particle of conserved spin (flavour) f.*/
        std::vector<size_t*> m_fs;
        /*! \brief Many particle state written in Fock space.
         * m_fock[f,st] is the index of the enumerated particle
         * with conserved spin (flavour) f in state st or m_Nfs[f]
         * if the state is unoccupied.*/
        std::vector<size_t*> m_fock;
        /*! \brief cache of all matrix elements.*/
        std::vector<std::complex<double>* > m_cache;
        std::vector<std::vector<std::vector<hop_path_t> > > m_exc;
        std::vector<std::vector<std::vector<int> > > m_fock_states;
        /*! \brief Current excitation index.*/
        size_t m_c_exc;
        size_t m_Nflav;//!< Number of spin flavours
        std::vector<size_t> m_Nfs;
        std::vector<size_t> m_Nby;
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
            matrix_element(size_t f, size_t r, size_t s)=0;

        hop_path_t get_hop_path(
                const size_t* src,
                const size_t* dest,
                size_t Nby,
                size_t Nfs) const;

        void get_hop_state(hop_path_t hop,
                             const size_t* src,
                             size_t Nby,
                             size_t Nfs,
                             std::vector<size_t>& dest) const;

    /**************
     * public:    *
     **************/
    public:
        /*! Base constructor*/
        WaveFunction_1(size_t Lx, size_t Ly,
                       std::vector<size_t> Nby,
                       std::vector<size_t> Nfs);

        virtual ~WaveFunction_1();

        std::complex<double> MatEl(const size_t& f,
                                   const size_t& r,
                                   const size_t& s) const
        {
            return m_cache[s][f*m_Lx*m_Ly+r];
        }
        size_t GetN(const size_t& s) const {return m_Nfs[s];}
        const size_t* GetFs(const size_t& s) const {return m_fs[s];}
        int GetSign() const {return m_sign;}
        size_t GetNExc() const {return m_exc[0].size();}

        // Writes down the states and relative fermi signs to storage
        void Save(FileManager* fm);

        /*! Fill up states after initialization.
         * **st are arrays
         * of Nidx*Nfs** elements specifying all the single
         * particle states.*/
        void fill(std::vector<size_t*> st);

        /*! Add a state to the wavefunction, defining
         * the hop path from and to the other already defined
         * states. fup and fdo are Fock states where m_Nfsup(do)
         * stands for empty in fup(do) state and any other
         * integer for occupied. The size of fup and fdo must
         * be m_Lx*m_Ly.
         */
        void AddState(const std::vector<const size_t*> f);

        void GetHops(std::vector<std::vector<hop_path_t> >& hop) const;

        void GetHop(size_t k, std::vector<hop_path_t>& hop) const;

        /*! do \f$|\Psi\rangle\rightarrow\
         * c^\dagger_{k_\uparrow'}\
         * c_{k_\uparrow}c^\dagger_{k_\downarrow'}\
         * c_{k_\downarrow}\
         * |\Psi\rangle\f$,
         * update matrix elements and determinant.*/
        void Hop(size_t khop);

        /* \brief Calculate the hopping fermionic sign */
        virtual int HopSign(const std::vector<hop_path_t>& hop) const;

        std::string Fock() const;
        friend std::ostream &
            operator<<(std::ostream& left,
                       const WaveFunction_1 & right);
};

#endif//_WAVEFUNCTION_H
