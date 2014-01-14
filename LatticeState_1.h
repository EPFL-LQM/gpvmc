#ifndef _LATTICESTATE_1_H
#define _LATTICESTATE_1_H

#include "defs.h"
#include "Lattice.h"
#include "State_1.h"
#include <iostream>

class LatticeState_1: public State_1{
    private:
        const Lattice* m_lattice;
        uint_vec_t m_Nifs; //!< Site internal degrees of freedom per flavour
    public:
        LatticeState_1(const Lattice* lattice,
                       const uint_vec_t& Npt,
                       const uint_vec_t& Nifs);
        size_t GetNsites() const;
        const Lattice* GetLattice() const;
        /*! \brief st[f][:] lists the occupied internal states at vertex v.
         */
        void GetLatOc(size_t v,
                      std::vector<uint_vec_t>& st);

        friend std::ostream& operator<<(std::ostream& out,const LatticeState_1& lst);
};

#endif//_LATTICESTATE_1_H
