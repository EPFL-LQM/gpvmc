#ifndef _SPINSTATE_H
#define _SPINSTATE_H
#include <cmath>
#include <cstring>
#include <complex>
#include <vector>
#include <utility>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>
#include "RanGen.h"
#include "BigComplex.h"

#ifndef _HOP_PATH_T
typedef std::vector<std::pair<size_t,size_t> > hop_path_t;
#endif

enum occu_t{
    UP=1,
    DOWN=0,
    EMPTY=2,
    DO=3
};

class Jastrow;

/*! \brief real space spin state class
 *
 */

class SpinState {
    private:
        size_t* m_up;
        size_t* m_do;
        size_t* m_latupid;
        size_t* m_latdoid;
        occu_t* m_latoc;
        size_t m_L;
        size_t m_Nup;
        size_t m_Ndo;
        Jastrow* m_jas;
        bool m_own_amp;
        int m_sign;
        int cdagsign(size_t rf, bool up);
        size_t count_nn() const;
        SpinState(const SpinState& sp);

    public:
        SpinState(size_t L, size_t Nup, size_t Ndown,
                  bool neel=false, bool doccu=false);
        ~SpinState();
        void Init(bool neel=false, bool doccu=false);
        void Init(size_t * state);

        size_t GetNup() const {return m_Nup;}
        size_t GetNdo() const {return m_Ndo;}
        size_t GetL() const {return m_L;}

        void SetJas(Jastrow* jas) {m_jas=jas;}
        const Jastrow* GetJas() const {return m_jas;}
        double Jas() const;

        const size_t* GetUp() const {return m_up;}
        const size_t* GetDo() const {return m_do;}

        //const size_t* GetLatUpId() const {return m_latupid;}
        //const size_t* GetLatDoId() const {return m_latdoid;}
        //const occu_t* GetLatOc() const {return m_latoc;}
        occu_t GetLatOc(size_t x, size_t y) const
        {
            return m_latoc[y*m_L+x];
        }
        size_t GetLatUpId(size_t x, size_t y) const
        {
            return m_latupid[y*m_L+x];
        }
        size_t GetLatDoId(size_t x, size_t y) const
        {
            return m_latdoid[y*m_L+x];
        }
        int GetSign() const {return m_sign;}

        void hop(const hop_path_t& hopup,
                 const hop_path_t& hopdo);
        int hop_sign(const hop_path_t& hopup,
                     const hop_path_t& hopdo) const;
        friend std::ostream& operator<<(std::ostream& out,const SpinState& sp);
};

#endif//_SPINSTATE_H
