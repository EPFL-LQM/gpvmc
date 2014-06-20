#ifndef _JASTROWTRACK_H
#define _JASTROWTRACK_H
#include "Quantity.h"
#include "BigDouble.h"
#include <vector>

/*!\brief class to track down the evolution of the MC weight
 *
 */
class JastrowTrack: public Quantity {
    public:
        JastrowTrack(const Stepper* stepper, FileManager* fm);
        virtual ~JastrowTrack();
        virtual void init();
        virtual void save() const;
        virtual void measure();
    protected:
        std::vector<BigDouble> m_vals;
};

#endif
