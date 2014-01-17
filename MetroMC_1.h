#ifndef _METROMC_1_H
#define _METROMC_1_H
#include <iostream>
#include <vector>
#include "BigDouble.h"

class Quantity_1;
class Stepper_1;
class FileManager;

class MetroMC_1 {
    protected:
        int m_comm_size;
        int m_comm_rank;
        BigDouble m_state_weight;
        double m_rwtimer;
        double m_gtimer;
        size_t m_steps;
        size_t m_rejection;
        Stepper_1* m_stepper;
        FileManager* m_fm;
        std::vector<Quantity_1*> m_quantity;
    public:
        //! \brief The constructor takes the stepper kernel. 
        MetroMC_1(Stepper_1* step, FileManager* fm);
        ~MetroMC_1() {}
        void AddQuantity(Quantity_1* quantity) {m_quantity.push_back(quantity);}
        std::vector<Quantity_1*>& GetQuantities() {return m_quantity;}
        const std::vector<Quantity_1*>& GetQuantities() const {return m_quantity;}
        double GetTimer() {return m_rwtimer;}
        /*! \brief Perform a random.
         * Perform a random walk of len steps and measure 
         * the quantity each meas steps.
         * If meas==0, no measurement is done.
         */
        void Walk(const size_t& len, size_t meas);
        /*! \bief Perform a step of the random walk.
         */
        void Step(bool meas);
        /*! \brief Return the proportion of rejected states.
         */
        double Rejection() const;
        /*! \brief YesNo random variable. p(Yes) = in
         */
        bool YesNo(const BigDouble& in);
};

#endif//_METROMC_H
