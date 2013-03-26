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
#include "Amplitude.h"// remove when cleaning up is complete

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
    return din>1.0 ? RanGen::uniform()<1 : RanGen::uniform() < din;
}

void MetroMC::Walk(const size_t& len, size_t meas, bool silent, int num_rep)
{
    double rwtimei, gtimei(0);
    m_gtimer=0;
    if(meas!=0) for(size_t q=0;q<m_quantity.size();++q)
        m_quantity[q]->init();
    for(size_t s=0;s<len;++s){
#ifdef USEPARA
        rwtimei=omp_get_wtime();
#else
        rwtimei=clock();
#endif
        if(meas && !silent) gtimei=rwtimei;
        Step(meas);
#ifdef USEPARA
        m_rwtimer+=omp_get_wtime()-rwtimei;
#else
        m_rwtimer+=(clock()-rwtimei)/double(CLOCKS_PER_SEC);
#endif
        if(meas && s%meas==0){
            for(size_t q=0; q<m_quantity.size();++q)
               m_quantity[q]->measure(); 
        }
#ifdef USEPARA
        if(!silent && meas) m_gtimer+=omp_get_wtime()-gtimei;
#else
        if(!silent && meas) m_gtimer+=(clock()-gtimei)/double(CLOCKS_PER_SEC);
#endif
        if(!silent && meas && len>=100 && s%(len/100)==0){
#ifdef USEMPI
            int mess=m_fm->message_monitor;
            double done=double(s)/len;
            double finishes=m_gtimer*len/s;
            MPI_Send(&mess,1,MPI_INT,0,0,MPI_COMM_WORLD);
            MPI_Send(&done,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            MPI_Send(&finishes,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
            MPI_Send(&num_rep,1,MPI_INT,0,0,MPI_COMM_WORLD);
#else
            vector<int> r(1,0);
            vector<double> d(1,double(s)/len);
            vector<double> f(1,m_gtimer*len/s);
            vector<int> n(1,num_rep);
            m_fm->Monitor(r,d,f,n);
#endif
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
    return double(m_rejection)/m_steps;
}

