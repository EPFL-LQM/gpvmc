#ifdef USEMPI
#include <mpi.h>
#endif
#include <iostream>
#include "SpinState.h"
#include "Amplitude.h"
#include "StagFluxGroundState.h"
#include "StagFluxTransExciton.h"
#include "StagFluxLongExciton.h"
#include "FullSpaceStepper.h"
#include "ProjHeis.h"
#include "StaggMagnJastrow.h"
#include "SpinDensity.h"
#include "FileManager.h"
#include "MetroMC.h"

using namespace std;

int main(int argc, char* argv[])
{
#ifdef USEMPI
    MPI_Init(&argc,&argv);
#endif
    size_t N=256;
    //int seed=1351781469;
    int seed=time(NULL);
    size_t L=4;
    RanGen::srand(seed);
    size_t q[2]={L/2,0};
    /*vector<vector<double> > qs(L/4*5+1,vector<double>(2,0));
    for(size_t qi=0;qi<L/4;++qi){
        qs[qi][0]=(L/4.0+qi)/L;
        qs[qi][1]=(L/4.0-qi)/L;
    }
    for(size_t qi=L/4;qi<L/4*3;++qi){
        qs[qi][0]=(L/2.0-(qi-L/4.0))/L;
    }
    for(size_t qi=L/4*3;qi<L/4*5+1;++qi){
        qs[qi][0]=(qi-L/4*3)/L;
        qs[qi][1]=(qi-L/4*3)/L;
    }*/
    //SpinState sp(L,L*L/2+1,L*L/2-1,false);
    SpinState sp(L,L*L/2,L*L/2,false);
    double bc[2]={1.0,1.0};
    FileManager fm;
    //StagFluxTransExciton ex(L,L,0.1*M_PI,0.0,bc,q);
    //StagFluxGroundState ex(L,L,0.1*M_PI,0.0,bc);
    StagFluxLongExciton ex(L,L,0.1*M_PI,0.0,bc,q);
    StaggMagnJastrow jas(&sp,0.08);
    Amplitude amp(&sp,&ex);
    for(size_t k=0;k<ex.GetNExc();++k){
        hop_path_t rhopup(1),rhopdo(1),khopup,khopdo;
        ex.GetHop(k,khopup,khopdo);
        ex.hop(k);
        amp.Update(rhopup,rhopdo,khopup,khopdo);
        cout<<"sign="<<ex.GetSign()<<endl;
        cout<<ex.Fock()<<endl;
        cout<<amp.Amp()<<endl;
    }
    /*FullSpaceStepper step(&amp);
    ProjHeis ph(&step,&fm);
    //SpinDensity sd(&step,&fm,qs,L);
    MetroMC mc(&step,&fm);
    mc.AddQuantity(&ph);
    //mc.AddQuantity(&sd);
    cout<<sp<<endl;
    Timer::tic("random walk");
    mc.Walk(N*L*L,0);
    mc.Walk(N*L*L,L*L,true);
    Timer::toc("random walk");
    //cout<<ph.Val(0,0)/N<<endl;
    //cout<<sd.str()<<endl;
    cout<<ph.str()<<endl;
    cout<<sp<<endl;
    //cout<<Timer::report()<<endl;*/
#ifdef USEMPI
    MPI_Finalize();
#endif
}
