#ifndef _WAVEFUNCTION_H
#define _WAVEFUNCTION_H

#include "State.h"
#include <complex>

class FileManager;

/*! \brief Base class for trial wave functions.
 *
 */

class WaveFunction : public State{

    /**************
     * protected: *
     **************/
    protected:
        /*! \brief cache of all matrix elements.*/
        std::vector<std::vector<std::complex<double> > > m_cache;
        /*! \brief matrix of paths linking all the added states together.
         * m_exc[i][f][fl] is the hoppings particles of flavour fl must do
         * such that one goes from state i to state f.
         */
        std::vector<std::vector<std::vector<hop_path_t> > > m_exc;
        /*! \brief Current excitation index.*/
        size_t m_c_exc;

        void build_base(const std::vector<size_t>& Nby, const std::vector<size_t>& Nfs);

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

    /**************
     * public:    *
     **************/
    public:
        /*! Base constructor*/
        WaveFunction();

        virtual ~WaveFunction();

        std::complex<double> MatEl(const size_t& f,
                                   const size_t& r,
                                   const size_t& s) const
        {
            return m_cache[s][f*m_Nfs[s]+r];
        }
        size_t GetNExc() const {return m_exc.size();}

        // Writes down the states and relative fermi signs to storage
        void Save(FileManager* fm);

        /*! Add a state to the wavefunction, defining
         * the hop path from and to the other already defined
         * states. fup and fdo are Fock states where m_Nfsup(do)
         * stands for empty in fup(do) state and any other
         * integer for occupied. The size of fup and fdo must
         * be m_Lx*m_Ly.
         */
        void AddState(const std::vector<uint_vec_t>& f);

        /*! \brief hop to khop'th state
         */
        void Hop(size_t khop);

        const std::vector<std::vector<hop_path_t> >& GetHops() const;

        const std::vector<hop_path_t>& GetHop(size_t k) const;

};

#endif//_WAVEFUNCTION_H
