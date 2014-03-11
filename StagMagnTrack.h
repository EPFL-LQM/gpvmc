#ifndef _STAGMAGNTRACK_H
#define _STAGMAGNTRACK_H
#include "Quantity.h"
#include <vector>
#include <complex>

/*!\brief class to track down the evolution of the staggered magnetization
 *
 */
class StagMagnTrack: public Quantity {
    public:
        StagMagnTrack(const Stepper* stepper, FileManager* fm);
        virtual ~StagMagnTrack();
        virtual void init();
        virtual void save() const;
        virtual void measure();
    protected:
        std::vector<std::complex<double> > m_vals;
};

#endif
