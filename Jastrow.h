#ifndef _JASTROW_H
#define _JASTROW_H

#include "defs.h"
#include "Amplitude.h"
#include "BigDouble.h"

class JastrowPotential;

class Jastrow: public Amplitude {
    public:
        Jastrow(const LatticeState* st,const JastrowPotential* pot);
        virtual void Init();
        virtual BigDouble Jas() const;

        virtual void VirtUpdate(const std::vector<std::vector<hop_path_t> >& rhop,
                                std::vector<BigDouble>& js) const;
        virtual void Update(const std::vector<hop_path_t>& rhop);
    protected:
        double m_jassum;
        const JastrowPotential* m_jaspot;
};

#endif//_JASTROW_H
