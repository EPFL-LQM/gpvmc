#ifndef _OVERLPATRACK_H
#define _OVERLAPTRACK_H
#include "Quantity.h"
#include "BigDouble.h"
#include <vector>

/*!\brief class to track down the evolution of the MC weight
 *
 */
class OverlapTrack: public Quantity {
    public:
        OverlapTrack(const Stepper* stepper, FileManager* fm);
        virtual ~OverlapTrack();
        virtual void init();
        virtual void save() const;
        virtual void measure();
    protected:
        std::vector<BigDouble> m_vals;
};

#endif
