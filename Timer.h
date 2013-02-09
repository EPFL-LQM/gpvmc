#ifndef _TIMER_H
#define _TIMER_H

#include <string>
#include <map>
#include <iostream>
#ifdef USEPARA
#include <omp.h>
#endif

/*! \brief A singleton class to add and manage timers
 *
 */

class Timer
{
    private:
        class tim{
            public:
                double m_t0;
                double m_t;
                std::map<std::string,tim*> m_timtree;
                tim(): m_t0(0), m_t(0){}
                ~tim();
                tim* get_timer(std::string timer);
                std::string report();
                double get_self_time();
        };

        static Timer s_ti;
        std::map<std::string,tim> m_timers;
        Timer& operator=(Timer&);
        Timer(Timer&);
        Timer() {}
    public:
        static void tic(std::string timer);
        static double toc(std::string timer);
        static const Timer& timer();
        static double timer(std::string ti);
        static void reset(std::string ti);
        /*! \brief report timings gathered in standard output.
         * In MPI programs, this must be collectively called.
         */
        static std::string report();
};

#endif//_TIMER_H
