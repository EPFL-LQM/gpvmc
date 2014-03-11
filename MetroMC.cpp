#ifdef USEMPI
#include <mpi.h>
#endif
#ifdef USEPARA
#include <omp.h>
#endif
#include <vector>
#include "MetroMC.h"
#include "FileManager.h"
#include "RanGen.h"
#include "Stepper.h"
#include "Quantity.h"
#include <unistd.h>

using namespace std;

MetroMC::MetroMC(Stepper* step, FileManager* fm)
    :m_comm_size(0), m_comm_rank(0), m_state_weight(-1), m_rwtimer(0), m_gtimer(0), m_steps(0), m_rejection(0), m_stepper(step), m_fm(fm), m_quantity(0)
{
#ifdef USEMPI
    MPI_Comm_size(MPI_COMM_WORLD,&m_comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&m_comm_rank);
#endif
    m_state_weight=m_stepper->weight();
}

bool MetroMC::YesNo(const BigDouble& in)
{
    double din(in);
#ifdef DEBUG
    if(din>1.0){
        return RanGen::uniform()<1;
    } else {
        double r=RanGen::uniform();
        if(r < din && din < 1e-6){
            cout<<"MetroMC::TesNo: Warning, highly improbable event: rand="<<r<<", rho="<<din<<endl;
        }
        return r < din;
    }
#else
    return din>1.0 ? RanGen::uniform()<1 : RanGen::uniform() < din;
#endif
}

void MetroMC::Walk(const size_t& len, size_t meas)
{
    double rwtimei, gtimei(0);
    m_gtimer=0;
    if(meas!=0)
        for(size_t q=0;q<m_quantity.size();++q)
            m_quantity[q]->init();
    for(size_t s=0;s<len;++s){
#ifdef USEPARA
        rwtimei=omp_get_wtime();
#else
        rwtimei=clock();
#endif
        if(meas) gtimei=rwtimei;
        Step(meas);
#ifdef USEPARA
        m_rwtimer+=omp_get_wtime()-rwtimei;
#else
        m_rwtimer+=(clock()-rwtimei)/double(CLOCKS_PER_SEC);
#endif
        if(meas && (meas==1 || s%meas==0)){
            for(size_t q=0; q<m_quantity.size();++q)
               m_quantity[q]->measure(); 
        }
#ifdef USEPARA
        if(meas) m_gtimer+=omp_get_wtime()-gtimei;
#else
        if(meas) m_gtimer+=(clock()-gtimei)/double(CLOCKS_PER_SEC);
#endif
        if(meas && len>=100 && s%(len/100)==0){
            double done=double(s)/len;
            double finishes=m_gtimer*len/s;
            if(s==0) finishes=-1;
            int comm_rank=0;
#ifdef USEMPI
            int mess=m_fm->message_monitor;
            MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
            //cout<<"rank "<<comm_rank<<": sends message_monitor"<<endl;
            MPI_Send(&mess,1,MPI_INT,0,m_fm->message_comm,MPI_COMM_WORLD);
#endif
            m_fm->Monitor(comm_rank,done,finishes);
        }
    }
}

void MetroMC::Step(bool meas)
{
    BigDouble new_weight(0);
    new_weight=m_stepper->trystep();
    //if(meas) std::cout<<sqrt(new_weight/m_state_weight*pow(10,new_exp-m_state_weight_exp))<<std::endl;
    if(YesNo(new_weight/m_state_weight)){
        m_stepper->step();
        m_state_weight=new_weight;
    } else {
        //std::cout<<"refused: "<<new_weight/m_state_weight<<std::endl;
        if(meas) m_rejection++;
    }
    if(meas) m_steps++;
}

double MetroMC::Rejection() const
{
    if(m_steps)
        return double(m_rejection)/m_steps;
    else
        return 0;
}

