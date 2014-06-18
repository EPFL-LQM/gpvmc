#ifndef _AMPLITUDE_H
#define _AMPLITUDE_H

class LatticeState;

class Amplitude {
    public:
        Amplitude(const LatticeState* latstate)
            :m_latstate(latstate) {}
        virtual void Init()=0;
    protected:
        const LatticeState* m_latstate;
};

#endif//_AMPLITUDE_H
