#ifndef _QUANTITY_H
#define _QUANTITY_H
#include <string>

class Stepper;
class FileManager;

/*! \brief Base class for quantities to be measured.
 *
 */

class Quantity {
    //for operators that are diagonal in real spin space.
    public:
        Quantity(const Stepper* stepper, FileManager* fm, std::string basename);
        virtual ~Quantity();
        virtual void init();
        virtual void measure();
        virtual void save() const=0;
    protected:
        int m_comm_size;
        int m_comm_rank;
        int m_meas;
        const Stepper* m_stepper;
        FileManager* m_fm;
        std::string m_basename;
};

#endif//_QUANTITY_H
