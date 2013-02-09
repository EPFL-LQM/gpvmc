#ifndef _RANGEN_H
#define _RANGEN_H

#include "Timer.h"

/*! \brief Random number generator. Uses gsl, mkl, acml or std depending of compile-time flags
 *
 */


struct state_t;

class RanGen
{
    private:
        static RanGen m_ran;
        state_t* m_rng;
        RanGen& operator=(RanGen&);
        RanGen(RanGen&);
        RanGen();
    public:
        ~RanGen();
        static void srand(size_t seed);
        static double uniform();
};

#endif//_RANGEN_H
