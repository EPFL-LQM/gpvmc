#ifndef _JASTROW_H
#define _JASTROW_H

#include "defs.h"
#include "Amplitude.h"

class JastrowPotential;

class Jastrow: public Amplitude {
    public:
        Jastrow(const LatticeState* st,const JastrowPotential* pot);
        virtual void Init();
        virtual double Jas() const;

        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<double>& js) const;
        virtual void Update(const std::vector<hop_path_t>& rhop);
    protected:
        double m_jassum;
        const JastrowPotential* m_jaspot;
};

#endif//_JASTROW_H
