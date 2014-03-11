#ifndef _PROJHEIS_H
#define _PROJHEIS_H

#include "MatrixQuantity.h"

class Lattice;

class ProjHeis: public MatrixQuantity
{
    public:
        ProjHeis(const Stepper* stepper,
                 FileManager* fm, const Lattice* lat,
                 double Bx=0.0);
        virtual ~ProjHeis() {}
        virtual void measure();
        std::string str() const;
    private:
        const Lattice* m_lat;
        double m_Bx;
};

#endif//_PROJHEIS_H
