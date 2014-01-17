#ifndef _QUANTITY_1_H
#define _QUANTITY_1_H
#include <string>

class Stepper_1;
class FileManager;

/*! \brief Base class for quantities to be measured.
 *
 */

class Quantity_1 {
    //for operators that are diagonal in real spin space.
    public:
        Quantity_1(const Stepper_1* stepper, FileManager* fm, std::string basename);
        virtual ~Quantity_1();
        virtual void init();
        virtual void measure();
        virtual void save() const=0;
    protected:
        int m_comm_size;
        int m_comm_rank;
        int m_meas;
        const Stepper_1* m_stepper;
        FileManager* m_fm;
        std::string m_basename;
};

#endif//_QUANTITY_H
