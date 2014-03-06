#ifndef _METROMC_H
#define _METROMC_H
#include <iostream>
#include <vector>
#include "BigDouble.h"

class Quantity;
class Stepper;
class FileManager;

class MetroMC {
    protected:
        int m_comm_size;
        int m_comm_rank;
        BigDouble m_state_weight;
        double m_rwtimer;
        double m_gtimer;
        size_t m_steps;
        size_t m_rejection;
        Stepper* m_stepper;
        FileManager* m_fm;
        std::vector<Quantity*> m_quantity;
    public:
        //! \brief The constructor takes the stepper kernel. 
        MetroMC(Stepper* step, FileManager* fm);
        ~MetroMC() {}
        void AddQuantity(Quantity* quantity) {m_quantity.push_back(quantity);}
        std::vector<Quantity*> GetQuantities() {return m_quantity;}
        const std::vector<Quantity*> GetQuantities() const {return m_quantity;}
        double GetTimer() {return m_rwtimer;}
        /*! \brief Perform a random.
         * Perform a random walk of len steps and measure 
         * the quantity each meas steps.
         * If meas==0, no measurement is done.
         */
        void Walk(const size_t& len, size_t meas);
        /*! \brief Perform a step of the random walk.
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
