#include "Amplitude.h"
#include <limits>
#include <complex>
#include <string>
#include <sstream>
#include <iomanip>
#include "linalg.h"
#include "Timer.h"
#include "blas_lapack.h"
#ifdef USEPARA
#include <omp.h>
#endif

using namespace std;

#define EPS -8
#define DEPS 1e-8

Amplitude::Amplitude(SpinState* sp, WaveFunction* wav)
    :m_sp(sp), m_wav(wav),
    m_matup(0),m_matiup(0),
    m_matdo(0),m_matido(0),
    m_amp(0), m_amp_ok(false), m_Nup(0), m_Ndo(0)
{
    if(sp->GetNup()!=wav->GetNup() || sp->GetNdo()!=wav->GetNdo())
#ifdef EXCEPT
        throw(std::logic_error("Amplitude::"
                    "Amplitude(const SpinState*, const WaveFunction*):"
                    "cannot define amplitude with different number of "
                    "spins up or down in wave-function and spin state."));
#else
        abort();
#endif
    m_Nup=m_sp->GetNup();
    m_Ndo=m_sp->GetNdo();
    m_matup=new std::complex<double>[m_Nup*m_Nup];
    m_matiup=new std::complex<double>[m_Nup*m_Nup];
    m_matdo=new std::complex<double>[m_Ndo*m_Ndo];
    m_matido=new std::complex<double>[m_Ndo*m_Ndo];
    Init();
}

Amplitude::~Amplitude()
{
    if(m_matup) delete [] m_matup;
    if(m_matiup) delete [] m_matiup;
    if(m_matdo) delete [] m_matdo;
    if(m_matido) delete [] m_matido;
}

void Amplitude::Init()
{
#ifdef PROFILE
    Timer::tic("Amplitude::Init");
#endif
#ifdef DEBUG
    std::cout<<"Amplitude::Init: has been called"<<std::endl;
#endif
    for(size_t rup=0;rup<m_Nup;++rup)
        for(size_t fup=0;fup<m_Nup;++fup)
            m_matup[rup*m_Nup+fup]=m_wav->MatEl(m_wav->GetUp()[fup],m_sp->GetUp()[rup],true);
    for(size_t rdo=0;rdo<m_Ndo;++rdo)
        for(size_t fdo=0;fdo<m_Ndo;++fdo)
            m_matdo[rdo*m_Ndo+fdo]=m_wav->MatEl(m_wav->GetDo()[fdo],m_sp->GetDo()[rdo],false);
    BigComplex dup(0.0), ddo(0.0);
    linalg::DetInv(m_matup, m_matiup,m_Nup, dup);
    linalg::DetInv(m_matdo, m_matido,m_Ndo, ddo);
    if(dup.exp()<EPS || ddo.exp()<EPS || dup==0.0 || ddo==0.0){
        m_amp_ok=false;
        m_amp=0.0;
    } else {
        // real space state sign omitted since it is canceled
        // in this calculation with the sign of the measured
        // quantity matrix elements
        m_amp=dup*ddo;
        m_amp_ok=true;
    }
#ifdef PROFILE
    Timer::toc("Amplitude::Init");
#endif
}

BigComplex Amplitude::Amp() const {return m_amp*m_wav->GetSign();}


BigComplex Amplitude::VirtColUpdate(const size_t& idup,
                                    const size_t& fup,
                                    const size_t& iddo,
                                    const size_t& fdo) const
{
    if(!m_amp_ok){
#ifdef EXCEPT
        throw(std::logic_error("Amplitude::VirtColUpdate:"
                          " Should not be called from"
                          " a state without overlap."));
#else
        abort();
#endif
    }
#ifdef PROFILE
    Timer::tic("Amplitude::VirtColUpdate");
#endif
    const std::complex<double> *akup=&m_matiup[idup*m_Nup];
    const std::complex<double> *akdo=&m_matido[iddo*m_Ndo];
    std::complex<double>* uup=new std::complex<double>[m_Nup];
    std::complex<double>* udo=new std::complex<double>[m_Ndo];
    std::complex<double> qup,qdo;
    for(size_t r=0;r<m_Nup;++r)
        uup[r]=m_wav->MatEl(fup,m_sp->GetUp()[r],true);
    for(size_t r=0;r<m_Ndo;++r)
        udo[r]=m_wav->MatEl(fdo,m_sp->GetDo()[r],false);
    linalg::zdotu_sub(m_Nup,akup,1,uup,1,&qup);
    linalg::zdotu_sub(m_Ndo,akdo,1,udo,1,&qdo);
    delete [] uup;
    delete [] udo;
#ifdef PROFILE
    Timer::toc("Amplitude::VirtColUpdate");
#endif
    hop_path_t hopup(1,pair<size_t,size_t>(idup,fup));
    hop_path_t hopdo(1,pair<size_t,size_t>(iddo,fdo));
    return m_amp*qup*qdo*m_wav->GetSign()
        *m_wav->hop_sign(hopup,hopdo);
}

void Amplitude::ColUpdate(const size_t& idup,
                          const size_t& iddo)
{
    if(!m_amp_ok){
        Init();
        return;
    }
#ifdef PROFILE
    Timer::tic("Amplitude::ColUpdate");
#endif
    size_t fup=m_wav->GetUp()[idup];
    size_t fdo=m_wav->GetDo()[iddo];
    std::complex<double> *akup=new std::complex<double>[2*m_Nup];
    memcpy(akup, &m_matiup[idup*m_Nup],m_Nup*sizeof(std::complex<double>));;
    std::complex<double> *akdo=new std::complex<double>[2*m_Ndo];
    memcpy(akdo, &m_matido[iddo*m_Ndo],m_Ndo*sizeof(std::complex<double>));
    std::complex<double> qup,qdo, one(1,0), zer(0,0);
    std::complex<double> *bup=&akup[m_Nup];
    std::complex<double> *bdo=&akdo[m_Ndo];
    for(size_t r=0;r<m_Nup;++r)
        m_matup[r*m_Nup+idup]=m_wav->MatEl(fup,m_sp->GetUp()[r],true);
    for(size_t r=0;r<m_Ndo;++r)
        m_matdo[r*m_Ndo+iddo]=m_wav->MatEl(fdo,m_sp->GetDo()[r],false);
    linalg::zdotu_sub(m_Nup,akup,1,&m_matup[idup],m_Nup,&qup);
    linalg::zdotu_sub(m_Ndo,akdo,1,&m_matdo[iddo],m_Ndo,&qdo);
    if(norm(qup*qdo)<DEPS){
        m_amp_ok=false;
        m_amp=0;
        delete [] akup;
        delete [] akdo;
#ifdef PROFILE
        Timer::toc("Amplitude::ColUpdate");
#endif
        return;
    }
    m_amp*=qup*qdo;
    qup=-1.0/qup;
    qdo=-1.0/qdo;
    linalg::zgemv('T',m_Nup,m_Nup,
                &one,m_matiup,m_Nup,&m_matup[idup],m_Nup,&zer,bup,1);
    linalg::zgemv('T',m_Ndo,m_Ndo,
                &one,m_matido,m_Ndo,&m_matdo[iddo],m_Ndo,&zer,bdo,1);
    bup[idup]-=1;
    bdo[iddo]-=1;
    linalg::zgeru(m_Nup,m_Nup,
                &qup,akup,1,bup,1,m_matiup,m_Nup);
    linalg::zgeru(m_Ndo,m_Ndo,
                &qdo,akdo,1,bdo,1,m_matido,m_Ndo);
    delete [] akup;
    delete [] akdo;
#ifdef PROFILE
    Timer::toc("Amplitude::ColUpdate");
#endif
}

BigComplex Amplitude::VirtRowUpdate(const size_t& idup,
                                    const size_t& rup,
                                    const size_t& iddo,
                                    const size_t& rdo) const
{
    if(!m_amp_ok){
#ifdef EXCEPT
        throw(std::logic_error("Amplitude::VirtRowUpdate:"
                          " Should not be called from"
                          " a state without overlap."));
#else
        abort();
#endif
    }
#ifdef PROFILE
    Timer::tic("Amplitude::VirtRowUpdate");
#endif
    const std::complex<double>* arup=&m_matiup[idup];
    const std::complex<double>* ardo=&m_matido[iddo];
    std::complex<double>* vup=new std::complex<double>[m_Nup];
    std::complex<double>* vdo=new std::complex<double>[m_Ndo];
    std::complex<double> qup,qdo;
    for(size_t f=0;f<m_Nup;++f)
        vup[f]=m_wav->MatEl(m_wav->GetUp()[f],rup,true);
    for(size_t f=0;f<m_Ndo;++f)
        vdo[f]=m_wav->MatEl(m_wav->GetDo()[f],rdo,false);
    linalg::zdotu_sub(m_Nup,vup,1,arup,m_Nup,&qup);
    linalg::zdotu_sub(m_Ndo,vdo,1,ardo,m_Ndo,&qdo);
    delete [] vup;
    delete [] vdo;
#ifdef PROFILE
    Timer::toc("Amplitude::VirtRowUpdate");
#endif
    return m_amp*qup*qdo*m_wav->GetSign();
}

void Amplitude::RowUpdate(const size_t& idup,
                          const size_t& iddo)
{
    if(!m_amp_ok){
        Init();
        return;
    }
#ifdef PROFILE
    Timer::tic("Amplitude::RowUpdate");
#endif
    std::complex<double>* arup=new std::complex<double>[2*m_Nup];
    for(size_t r=0;r<m_Nup;++r) arup[r]=m_matiup[r*m_Nup+idup];
    std::complex<double>* ardo=new std::complex<double>[2*m_Ndo];
    for(size_t r=0;r<m_Ndo;++r) ardo[r]=m_matido[r*m_Ndo+iddo];
    std::complex<double>* bup=&arup[m_Nup];
    std::complex<double>* bdo=&ardo[m_Ndo];
    std::complex<double> qup,qdo,one(1,0),zer(0,0);
    for(size_t f=0;f<m_Nup;++f)
        m_matup[idup*m_Nup+f]=m_wav->MatEl(m_wav->GetUp()[f],m_sp->GetUp()[idup],true);
    for(size_t f=0;f<m_Ndo;++f)
        m_matdo[iddo*m_Ndo+f]=m_wav->MatEl(m_wav->GetDo()[f],m_sp->GetDo()[iddo],false);
    linalg::zdotu_sub(m_Nup,&m_matup[idup*m_Nup],1,arup,1,&qup);
    linalg::zdotu_sub(m_Ndo,&m_matdo[iddo*m_Ndo],1,ardo,1,&qdo);
    if(norm(qup*qdo)<DEPS){
        m_amp_ok=false;
        m_amp=0;
        delete [] arup;
        delete [] ardo;
#ifdef PROFILE
        Timer::toc("Amplitude::RowUpdate");
#endif
        return;
    }
    m_amp*=qup*qdo;
    qup=-1.0/qup;
    qdo=-1.0/qdo;
    linalg::zgemv('N',m_Nup,m_Nup,
                &one,m_matiup,m_Nup,&m_matup[idup*m_Nup],1,&zer,bup,1);
    linalg::zgemv('N',m_Ndo,m_Ndo,
                &one,m_matido,m_Ndo,&m_matdo[iddo*m_Ndo],1,&zer,bdo,1);
    bup[idup]-=1;
    bdo[iddo]-=1;
    linalg::zgeru(m_Nup,m_Nup,
                &qup,bup,1,arup,1,m_matiup,m_Nup);
    linalg::zgeru(m_Ndo,m_Ndo,
                &qdo,bdo,1,ardo,1,m_matido,m_Ndo);
    delete [] arup;
    delete [] ardo;
#ifdef PROFILE
    Timer::toc("Amplitude::RowUpdate");
#endif
}

void Amplitude::VirtUpdate(
        const std::vector<hop_path_t> rhopup,
        const std::vector<hop_path_t> rhopdo,
        const std::vector<hop_path_t> khopup,
        const std::vector<hop_path_t> khopdo,
        BigComplex* qs) const
{
    if(!m_amp_ok)
#ifdef EXCEPT
        throw(std::logic_error("Amplitude::VirtUpdate:"
                               " Should not be called from"
                               " a state without overlap."));
#else
        abort();
#endif
    size_t Nr=rhopup.size();
    size_t Nk=khopup.size();
    if(Nr*Nk==0)
#ifdef EXCEPT
        throw(std::logic_error("Amplitude::VirtUpdate:"
                               " the condition min(Nr)=1 "
                               " and min(Nk)=1"
                               " must be fullfilled."));
#else
        abort();
#endif
#ifdef PROFILE
    Timer::tic("Amplitude::VirtUpdate");
#endif
    size_t NN=max(Nr,Nk);
    // Create vectors of block coordinates
    vector<size_t> ridxup(Nr,0);
    vector<size_t> ridxdo(Nr,0);
    vector<size_t> kidxup(Nk,0);
    vector<size_t> kidxdo(Nk,0);
    for(size_t n=1;n<Nr;++n){
        ridxup[n]=ridxup[n-1]+rhopup[n-1].size();
        ridxdo[n]=ridxdo[n-1]+rhopdo[n-1].size();
    }
    for(size_t n=1;n<Nk;++n){
        kidxup[n]=kidxup[n-1]+khopup[n-1].size();
        kidxdo[n]=kidxdo[n-1]+khopdo[n-1].size();
    }
    size_t Nrup=ridxup.back()+rhopup.back().size();
    size_t Nrdo=ridxdo.back()+rhopdo.back().size();
    size_t Nkup=kidxup.back()+khopup.back().size();
    size_t Nkdo=kidxdo.back()+khopdo.back().size();
    // allocate block matrices of change.
    complex<double> *vup=new complex<double>[m_Nup*Nrup+
                                             m_Ndo*Nrdo+
                                             Nkup*m_Nup+
                                             Nkdo*m_Ndo];
    std::complex<double> *vdo=&vup[m_Nup*Nrup];
    std::complex<double> *uup=&vup[m_Nup*Nrup+
                                   m_Ndo*Nrdo];
    std::complex<double> *udo=&vup[m_Nup*Nrup+
                                   m_Ndo*Nrdo+
                                   Nkup*m_Nup];
    /* copy column and row changes as they would
     * happen separatly (crossings between row
     * and column changes are treated later
     * below).
     */
    for(size_t n=0;n<Nr;++n){
        for(size_t r=0;r<rhopup[n].size();++r){
            for(size_t f=0;f<m_Nup;++f)
                vup[(ridxup[n]+r)*m_Nup+f]=
                    m_wav->MatEl(m_wav->GetUp()[f],
                                 rhopup[n][r].second,true);
        }
        for(size_t r=0;r<rhopdo[n].size();++r)
            for(size_t f=0;f<m_Ndo;++f)
                vdo[(ridxdo[n]+r)*m_Ndo+f]=
                    m_wav->MatEl(m_wav->GetDo()[f],
                                 rhopdo[n][r].second,false);
    }
    for(size_t n=0;n<Nk;++n){
        for(size_t r=0;r<khopup[n].size();++r){
            for(size_t f=0;f<m_Nup;++f){
                uup[f*Nkup+kidxup[n]+r]=
                    m_wav->MatEl(khopup[n][r].second,
                                 m_sp->GetUp()[f],true);
            }
        }
        for(size_t r=0;r<khopdo[n].size();++r)
            for(size_t f=0;f<m_Ndo;++f)
                udo[f*Nkdo+kidxdo[n]+r]=
                    m_wav->MatEl(khopdo[n][r].second,
                                 m_sp->GetDo()[f],false);
    }
    BigComplex* qup=new BigComplex[2*NN];
    BigComplex* qdo=&qup[NN];
    /*
     * ugly shit, I need to implement a contiguous container
     * for the hop_path_t.
     */
    vector<size_t> rupid(Nrup);
    vector<size_t> rdoid(Nrdo);
    vector<size_t> kupid(Nkup);
    vector<size_t> kdoid(Nkdo);
    vector<size_t> kup(Nkup);
    vector<size_t> kdo(Nkdo);
    vector<size_t> rupr(Nr);
    vector<size_t> rdor(Nr);
    vector<size_t> kupr(Nk);
    vector<size_t> kdor(Nk);
    size_t kupmax(0),kdomax(0);
    for(size_t nr=0;nr<Nr;++nr){
        rupr[nr]=rhopup[nr].size();
        for(size_t rr=0;rr<rhopup[nr].size();++rr)
            rupid[ridxup[nr]+rr]=rhopup[nr][rr].first;
        rdor[nr]=rhopdo[nr].size();
        for(size_t rr=0;rr<rhopdo[nr].size();++rr)
            rdoid[ridxdo[nr]+rr]=rhopdo[nr][rr].first;
    }
    for(size_t nk=0;nk<Nk;++nk){
        kupr[nk]=khopup[nk].size();
        if(kupr[nk]>kupmax) kupmax=kupr[nk];
        for(size_t rk=0;rk<khopup[nk].size();++rk){
            kupid[kidxup[nk]+rk]=khopup[nk][rk].first;
            kup[kidxup[nk]+rk]=khopup[nk][rk].second;
        }
        kdor[nk]=khopdo[nk].size();
        if(kdor[nk]>kdomax) kdomax=kdor[nk];
        for(size_t rk=0;rk<khopdo[nk].size();++rk){
            kdoid[kidxdo[nk]+rk]=khopdo[nk][rk].first;
            kdo[kidxdo[nk]+rk]=khopdo[nk][rk].second;
        }
    }
    /* treat row and columns crossings and then update
     * the determinants. As linalg::DetUpdate is efficient
     * for large multiple updates, I use NN=max(Nk,Nr) for
     * simultaneous update.
     */
    if(Nk==NN){
        for(size_t nr=0;nr<Nr;++nr){
            for(size_t rr=0;rr<rhopup[nr].size();++rr)
                for(size_t nk=0;nk<Nk;++nk)
                    for(size_t rk=0;rk<khopup[nk].size();++rk)
                        uup[rupid[ridxup[nr]+rr]*Nkup+
                            kidxup[nk]+rk]=
                            m_wav->MatEl(khopup[nk][rk].second,
                                         rhopup[nr][rr].second,
                                         true);
            for(size_t rr=0;rr<rhopdo[nr].size();++rr)
                for(size_t nk=0;nk<Nk;++nk)
                    for(size_t rk=0;rk<khopdo[nk].size();++rk)
                        udo[rdoid[ridxdo[nr]+rr]*Nkdo+
                            kidxdo[nk]+rk]=
                            m_wav->MatEl(khopdo[nk][rk].second,
                                         rhopdo[nr][rr].second,
                                         false);
#ifdef PROFILE
            Timer::tic("Amplitude::VirtUpdate/DetUpdate");
#endif
            /*
             * ugly again, should modify for better
             * interface.
             */
            linalg::DetUpdate(m_matup,m_matiup,m_Nup,
                              &vup[ridxup[nr]*m_Nup],m_Nup,
                              &rupid[ridxup[nr]],&rupr[nr],1,
                              uup,Nkup,&kupid[0],&kupr[0],Nk,
                              qup);
            linalg::DetUpdate(m_matdo,m_matido,m_Ndo,
                              &vdo[ridxdo[nr]*m_Ndo],m_Ndo,
                              &rdoid[ridxdo[nr]],&rdor[nr],1,
                              udo,Nkdo,&kdoid[0],&kdor[0],Nk,
                              qdo);
#ifdef PROFILE
            Timer::toc("Amplitude::VirtUpdate/DetUpdate");
#endif
            for(size_t nk=0;nk<Nk;++nk){
                qs[nk*Nr+nr]=qup[nk]*qdo[nk]*m_amp*
                    m_wav->GetSign()*
                    m_wav->hop_sign(khopup[nk],khopdo[nk]);
            }
            for(size_t rr=0;rr<rhopup[nr].size();++rr)
                for(size_t nk=0;nk<Nk;++nk)
                    for(size_t rk=0;rk<khopup[nk].size();++rk)
                        uup[rupid[ridxup[nr]+rr]*Nkup+
                            kidxup[nk]+rk]=
                            m_wav->MatEl(khopup[nk][rk].second,
                                     m_sp->GetUp()[
                                         rupid[ridxup[nr]+rr]],
                                     true);
            for(size_t rr=0;rr<rhopdo[nr].size();++rr)
                for(size_t nk=0;nk<Nk;++nk)
                    for(size_t rk=0;rk<khopdo[nk].size();++rk)
                        udo[rdoid[ridxdo[nr]+rr]*Nkdo+
                            kidxdo[nk]+rk]=
                            m_wav->MatEl(khopdo[nk][rk].second,
                                     m_sp->GetDo()[
                                         rdoid[ridxdo[nr]+rr]],
                                     false);
        }
    } else {
        for(size_t nk=0;nk<Nk;++nk){
            for(size_t rk=0;rk<khopup[nk].size();++rk)
                for(size_t nr=0;nr<Nr;++nr)
                    for(size_t rr=0;rr<rhopup[nr].size();++rr)
                        vup[(ridxup[nr]+rr)*m_Nup+
                            kupid[kidxup[nk]+rk]]=
                            m_wav->MatEl(khopup[nk][rk].second,
                                         rhopup[nr][rr].second,
                                         true);
            for(size_t rk=0;rk<khopdo[nk].size();++rk)
                for(size_t nr=0;nr<Nr;++nr)
                    for(size_t rr=0;rr<rhopdo[nr].size();++rr)
                        vdo[(ridxdo[nr]+rr)*m_Ndo+
                            kdoid[kidxdo[nk]+rk]]=
                            m_wav->MatEl(khopdo[nk][rk].second,
                                         rhopdo[nr][rr].second,
                                         false);
#ifdef PROFILE
            Timer::tic("Amplitude::VirtUpdate/DetUpdate");
#endif
            linalg::DetUpdate(m_matup,m_matiup,m_Nup,
                              vup,m_Nup,&rupid[0],&rupr[0],Nr,
                              &uup[kidxup[nk]],Nkup,&kupid[kidxup[nk]],
                              &kupr[nk],1,qup);
            linalg::DetUpdate(m_matdo,m_matido,m_Ndo,
                              vdo,m_Ndo,&rdoid[0],&rdor[0],Nr,
                              &udo[kidxdo[nk]],Nkdo,&kdoid[kidxdo[nk]],
                              &kdor[nk],1,qdo);
#ifdef PROFILE
            Timer::toc("Amplitude::VirtUpdate/DetUpdate");
#endif
            for(size_t nr=0;nr<Nr;++nr){
                qs[nk*Nr+nr]=qup[nr]*qdo[nr]*m_amp*
                    m_wav->GetSign()*
                    m_wav->hop_sign(khopup[nk],khopdo[nk]);
            }
            for(size_t rk=0;rk<khopup[nk].size();++rk)
                for(size_t nr=0;nr<Nr;++nr)
                    for(size_t rr=0;rr<rhopup[nr].size();++rr)
                        vup[(ridxup[nr]+rr)*m_Nup+
                            kupid[kidxup[nk]+rk]]=
                            m_wav->MatEl(m_wav->GetUp()
                                    [kupid[kidxup[nk]+rk]],
                                         rhopup[nr][rr].second,
                                         true);
            for(size_t rk=0;rk<khopdo[nk].size();++rk)
                for(size_t nr=0;nr<Nr;++nr)
                    for(size_t rr=0;rr<rhopdo[nr].size();++rr)
                        vdo[(ridxdo[nr]+rr)*m_Ndo+
                            kdoid[kidxdo[nk]+rk]]=
                            m_wav->MatEl(m_wav->GetDo()
                                    [kdoid[kidxdo[nk]+rk]],
                                         rhopdo[nr][rr].second,
                                         false);
        }
    }
    delete [] vup;
    delete [] qup;
#ifdef PROFILE
    Timer::toc("Amplitude::VirtUpdate");
#endif
}

void Amplitude::Update(const hop_path_t rhopup,
                       const hop_path_t rhopdo,
                       const hop_path_t khopup,
                       const hop_path_t khopdo)
{
    if(!m_amp_ok){
#ifdef EXCEPT
        throw(std::logic_error("Amplitude::Update:"
                               " Should not be called from"
                               " a state without overlap."));
#else
        cerr<<"Amplitude::Update:"
              " Should not be called from"
              " a state without overlap."<<endl;
        abort();
#endif
    }
#ifdef PROFILE
    Timer::tic("Amplitude::Update");
#endif
    std::complex<double> *vup=new std::complex<double>[
                                        m_Nup*rhopup.size()+
                                        m_Ndo*rhopdo.size()+
                                        m_Nup*khopup.size()+
                                        m_Ndo*khopdo.size()];
    std::complex<double> *vdo=&vup[m_Nup*rhopup.size()];
    for(size_t r=0;r<rhopup.size();++r)
        for(size_t f=0;f<m_Nup;++f)
            vup[r*m_Nup+f]=m_wav->MatEl(m_wav->GetUp()[f],
                    m_sp->GetUp()[rhopup[r].first],true);
    for(size_t r=0;r<rhopdo.size();++r)
        for(size_t f=0;f<m_Ndo;++f)
            vdo[r*m_Ndo+f]=m_wav->MatEl(m_wav->GetDo()[f],
                    m_sp->GetDo()[rhopdo[r].first],false);
    std::complex<double> *uup=&vup[m_Nup*rhopup.size()+
                                   m_Ndo*rhopdo.size()];
    std::complex<double> *udo=&vup[m_Nup*rhopup.size()+
                                   m_Ndo*rhopdo.size()+
                                   m_Nup*khopup.size()];
    for(size_t r=0;r<khopup.size();++r)
        for(size_t f=0;f<m_Nup;++f)
            uup[f*khopup.size()+r]=m_wav->MatEl(
                    m_wav->GetUp()[khopup[r].first],
                    m_sp->GetUp()[f],true);
    for(size_t r=0;r<khopdo.size();++r)
        for(size_t f=0;f<m_Ndo;++f)
            udo[f*khopdo.size()+r]=m_wav->MatEl(
                    m_wav->GetDo()[khopdo[r].first],
                    m_sp->GetDo()[f],false);
    BigComplex qup,qdo;
    /* again bad interface, need contiguous hop_path
     */
    vector<size_t> kupid(khopup.size());
    vector<size_t> kdoid(khopdo.size());
    vector<size_t> rupid(rhopup.size());
    vector<size_t> rdoid(rhopdo.size());
    for(size_t r=0;r<khopup.size();++r)
        kupid[r]=khopup[r].first;
    for(size_t r=0;r<khopdo.size();++r)
        kdoid[r]=khopdo[r].first;
    for(size_t r=0;r<rhopup.size();++r)
        rupid[r]=rhopup[r].first;
    for(size_t r=0;r<rhopdo.size();++r)
        rdoid[r]=rhopdo[r].first;
#ifdef PROFILE
    Timer::tic("Amplitude::Update/InvUpdate");
#endif
    linalg::InvUpdate(m_matup,m_matiup,m_Nup,
                      vup,&rupid[0],rhopup.size(),
                      uup,&kupid[0],khopup.size(),qup);
    linalg::InvUpdate(m_matdo,m_matido,m_Ndo,
                      vdo,&rdoid[0],rhopdo.size(),
                      udo,&kdoid[0],khopdo.size(),qdo);
#ifdef PROFILE
    Timer::toc("Amplitude::Update/InvUpdate");
#endif
    m_amp*=qup*qdo;
    delete [] vup;
#ifdef PROFILE
    Timer::toc("Amplitude::Update");
#endif
}

