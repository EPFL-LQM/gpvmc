#ifndef _LATTICESTATE_H
#define _LATTICESTATE_H

#include "defs.h"
#include "State.h"
#include <iostream>

class Lattice;

class LatticeState: public State{
    private:
        const Lattice* m_lattice;
        uint_vec_t m_Nifs; //!< Site internal degrees of freedom per flavour
    public:
        LatticeState(const Lattice* lattice,
                       const uint_vec_t& Npt,
                       const uint_vec_t& Nifs);
        size_t GetNsites() const;
        const uint_vec_t& GetNifs() const;
        const Lattice* GetLattice() const;
        /*! \brief Initialize at random the state
         */
        void RanInit();
        /*! \brief st[f][:] lists the occupied internal states at vertex v.
         */
        void GetLatOc(size_t v,
                      std::vector<uint_vec_t>& st) const;

        friend std::ostream& operator<<(std::ostream& out,const LatticeState& lst);
};

#endif//_LATTICESTATE_H
