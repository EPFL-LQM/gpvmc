#include "RanGen.h"

#ifdef USE_RNG_MKL
    #ifdef __cplusplus
        extern "C" {
    #endif//__cplusplus
    #include <mkl_vsl.h>
    #ifdef __cplusplus
        }
    #endif//__cplusplus
#elif defined USE_RNG_GSL
    #include <gsl/gsl_rng.h>
#elif defined USE_RNG_ACML
        namespace acml {
            #include <acml.h>
    }
#elif defined USE_RNG_STD
    #include <cstdlib>
#endif//USE_RNG_MKL

struct state_t{
#ifdef USE_RNG_MKL
        VSLStreamStatePtr state;
#elif defined USE_RNG_ACML
        int *state;
#elif defined USE_RNG_GSL
        gsl_rng* state;
#endif
};


RanGen::RanGen()
{
    m_rng=new state_t;
#ifdef USE_RNG_MKL
    vslNewStream(&m_rng->state,VSL_BRNG_MT19937,1);
#elif defined USE_RNG_ACML
    m_rng->state=new int[650];
#elif defined USE_RNG_GSL
    m_rng->state=gsl_rng_alloc(gsl_rng_ranlux);
    gsl_rng_set(m_rng->state,1);
#endif
}

RanGen::~RanGen()
{
#ifdef USE_RNG_MKL
    vslDeleteStream(&m_rng->state);
#elif defined USE_RNG_ACML
    delete [] m_rng->state;
#elif defined USE_RNG_MKL
    gsl_rng_free(m_rng->state);
#endif
    delete m_rng;
}

void RanGen::srand(size_t seed)
{
#ifdef USE_RNG_MKL
    vslDeleteStream(&m_ran.m_rng->state);
    vslNewStream(&m_ran.m_rng->state,VSL_BRNG_MT19937,seed);
#elif defined USE_RNG_ACML
    int sseed(seed),info(0),lseed(1), lstate(650);
    acml::drandinitialize(3,0,&sseed,&lseed,m_ran.m_rng->state,&lstate,&info);
#elif defined USE_RNG_GSL
    gsl_rng_set(m_ran.m_rng->state,seed);
#elif defined USE_RNG_STD
    std::srand(seed);
#endif
}

double RanGen::uniform()
{
    double out(0);
#ifdef USE_RNG_MKL
    vdRngUniform(VSL_METHOD_DUNIFORM_STD_ACCURATE,m_ran.m_rng->state,1,&out,0,1);
#elif defined USE_RNG_ACML
    int info;
    acml::dranduniform(1,0.0,1.0,m_ran.m_rng->state,&out,&info);
#elif defined USE_RNG_GSL
    out=gsl_rng_uniform(m_ran.m_rng->state);
#elif defined USE_RNG_STD
    out=((double)std::rand())/((double)RAND_MAX);
#endif
    return out;
}

RanGen RanGen::m_ran;
