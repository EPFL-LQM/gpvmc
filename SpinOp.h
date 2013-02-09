#ifndef _SPINOP_H
#define _SPINOP_H

#include "BigComplex.h"
#include <vector>
#include <complex>

class Amplitude;

//! total Sz conserving operators
/*! return \f$\langle\alpha|\hat{O}|\alpha'\rangle\langle \psi|\alpha'\rangle\f$
 * with \f$|\alpha'\rangle\f$ specified by swapping coordinates for instance.
 */
namespace SpinOp{
    BigComplex SipSjm(const Amplitude* amp,
                      const size_t& ix, const size_t& iy,
                      const size_t& jx, const size_t& jy);
    BigComplex SimSjp(const Amplitude* amp,
                      const size_t& ix, const size_t& iy,
                      const size_t& jx, const size_t& jy);
    BigComplex SizSjz(const Amplitude* amp,
                      const size_t& ix, const size_t& iy,
                      const size_t& jx, const size_t& jy);
    BigComplex Siz(const Amplitude* amp,
                   const size_t& ix, const size_t& i);
}


#endif//_SPINOP_H
