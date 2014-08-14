#ifndef _LATTICESTATE_H
#define _LATTICESTATE_H

#include "defs.h"
#include "State.h"
#include <iostream>

class Lattice;
class FileManager;

class LatticeState: public State{
    private:
        const Lattice* m_lattice;
        uint_vec_t m_Nifs; //!< Site internal degrees of freedom per flavour
    public:
        LatticeState(FileManager* fm,
                     const Lattice* lattice,
                     const uint_vec_t& Npt,
                     const uint_vec_t& Nifs);
        size_t GetNsites() const;
        const uint_vec_t& GetNifs() const;
        const Lattice* GetLattice() const;
        /*! \brief Initialize at random the state
         */
        void RanInit(const std::vector<std::vector<size_t> >& pop);
        /*! \brief Initialize from a fock state (0=empty, 1 occupied)
         */
        void FockInit(const std::vector<std::vector<int> >& fock);
        /*! \brief st[f][:] lists the occupied internal states at vertex v.
         */
        void GetLatOc(size_t v,
                      std::vector<uint_vec_t>& st) const;
        /*! \brief Translate the fock index fidx into the corresponding quantum numbers.
         */
        void Fock2QN(size_t fidx,size_t flav, uint_vec_t& st) const;

        friend std::ostream& operator<<(std::ostream& out,const LatticeState& lst);
};

#endif//_LATTICESTATE_H
