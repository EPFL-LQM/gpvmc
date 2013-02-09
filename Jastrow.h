#ifndef _JASTROW_H
#define _JASTROW_H

#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

#ifndef _HOP_PATH_T
typedef std::vector<std::pair<size_t,size_t> > hop_path_t;
#endif

class SpinState;

class Jastrow {
    public:
        Jastrow(SpinState* sp);
        virtual ~Jastrow();
        double Jas() const {return m_jas;}
        virtual void hop(const hop_path_t& hopup,
                         const hop_path_t& hopdo);
        virtual double virtualhop(const hop_path_t& hopup,
                                  const hop_path_t& hopdo) const;
        virtual void init();
        friend std::ostream& operator<<(std::ostream& os,const Jastrow& j);
    protected:
        virtual double jpot(const size_t* r)=0;
        const SpinState* m_sp;
        double* m_pot;
        double m_jas;
};

#endif//_JASTROW_H
