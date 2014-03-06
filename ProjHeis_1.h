#ifndef _PROJHEIS_H
#define _PROJHEIS_H

#include "MatrixQuantity_1.h"
#include "Lattice.h"

class ProjHeis_1: public MatrixQuantity_1
{
    public:
        ProjHeis_1(const Stepper_1* stepper,
                 FileManager* fm, const Lattice* lat,
                 double Bx=0.0);
        virtual ~ProjHeis_1() {}
        virtual void measure();
        std::string str() const;
    private:
        const Lattice* m_lat;
        double m_Bx;
};

#endif//_PROJHEIS_H
