#ifndef _PROJHEIS_H
#define _PROJHEIS_H

#include "MatrixQuantity.h"

class ProjHeis: public MatrixQuantity
{
    public:
        ProjHeis(const Stepper* stepper,
                 FileManager* fm,
                 double jr=0);
        virtual ~ProjHeis() {}
        virtual void measure();
        std::string str() const;
    private:
        double m_jr;
};

#endif//_PROJHEIS_H
