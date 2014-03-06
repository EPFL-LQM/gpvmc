#ifndef _SQUARELATIICE_H
#define _SQUARELATTICE_H

#include <string>
#include "Lattice.h"

class SquareLattice: public Lattice {
    public:
        SquareLattice(size_t Lx, size_t Ly);
        virtual std::string str(std::vector<std::string> st=std::vector<std::string>()) const;
};

#endif//_SQUARELATTICE_H
