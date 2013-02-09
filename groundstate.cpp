#ifdef USEMPI
#include <mpi.h>
#endif
#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <vector>
#include "MetroMC.h"
#include "SpinState.h"
#include "StaggMagnJastrow.h"
#include "StagJastrow.h"
#include "StagFluxGroundState.h"
#include "FullSpaceStepper.h"
#include "Amplitude.h"
#include "RanGen.h"
#include "ProjHeis.h"
#include "SpinDensity.h"
#include "FileManager.h"
#include "ArgParse.h"
#include "SpinSpinCorr.h"
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

using namespace std;

int main(int argc, char* argv[])
{
    // Initialize
#ifdef USEMPI
    MPI_Init(&argc,&argv);
#endif
    int comm_size(1);
    int comm_rank(0);
#ifdef USEMPI
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
#endif
    signal(SIGTERM,FileManager::EmergencyClose);
    RanGen::srand(time(NULL)+100*comm_rank);
    Timer::tic("main");

    // Parse input arguments
    ArgParse arg(argc,argv);
    size_t L=arg.i("L");
    double phi=arg.d("phi");
    double neel=arg.d("neel");
    double jastrow=arg.d("jastrow");
    double N=arg.i("N");
    size_t s=arg.i("s");
    size_t saves=arg.i("saves");
    int prefix=arg.i("prefix");
    string dir=arg.s("dir");
    string qfilestr=arg.s("qfile");
    double therm=arg.i("therm");
    double phase_shift[2]={arg.d("phase_shift_x"),arg.d("phase_shift_y")};
    bool jas_stagmagn=arg.b("jas_onebodystag");
    bool jas_stag=arg.b("jas_twobodystag");
    double jr=arg.d("Jr");

    // Setup calculation parameters
    FileManager fm(dir,prefix);
    phi*=M_PI;
    double rej(0);
    std::vector<std::vector<double> > qs;
    std::ifstream qfile(qfilestr.c_str());
    if(!qfilestr.empty() && qfile.is_open()){
        while(qfile.good()){
            double qx,qy;
            qfile>>qx>>qy;
            if(qfile.good()){
                qs.push_back(std::vector<double>(2,0));
                qs.back()[0]=qx*2*M_PI;
                qs.back()[1]=qy*2*M_PI;
            }
        }
    }
    size_t* Rs=new size_t[2*(L/2+1)*(L/2+1)];
    for(size_t x=0;x<=L/2;++x){
        for(size_t y=0;y<=L/2;++y){
            Rs[(x*(L/2+1)+y)*2]=x;
            Rs[(x*(L/2+1)+y)*2+1]=y;
        }
    }
    fm.FileAttribute("L",L);
    fm.FileAttribute("N",N);
    fm.FileAttribute("neel",neel);
    fm.FileAttribute("phi",phi/M_PI);
    fm.FileAttribute("jastrow",jastrow);
    fm.FileAttribute("s",s);
    fm.FileAttribute("saves",saves);
    fm.FileAttribute("phasex",phase_shift[0]);
    fm.FileAttribute("phasey",phase_shift[1]);
    fm.FileAttribute("Jr",jr);
    if(jas_stag) fm.FileAttribute("twobodystagjastrow",1);
    if(jas_stagmagn) fm.FileAttribute("onebodystagjastrow",1);

    fm.StatPerSample()=saves/s;
#ifdef USEMPI
    if(comm_rank){
#endif
        // Setup calculation
        SpinState sp(L,L*L/2,L*L/2);
        Jastrow* jas=0;
        if(jastrow!=0){
            if(jas_stagmagn)
                jas=new StaggMagnJastrow(&sp,jastrow);
            else if(jas_stag)
                jas=new StagJastrow(&sp,jastrow);
        }
        StagFluxGroundState wav(L,L,phi,neel,phase_shift);
        Amplitude amp(&sp,&wav);
        FullSpaceStepper step(&amp);
        ProjHeis heisen(&step,&fm,jr);
        SpinSpinCorr sscorr(&step,&fm,(L/2+1)*(L/2+1),Rs);
        SpinDensity spindense(&step,&fm,qs,L);
        MetroMC varmc(&step,&fm);
        varmc.AddQuantity(&heisen);
        varmc.AddQuantity(&sscorr);
        if(qs.size()) varmc.AddQuantity(&spindense);
        while(amp.Amp()==0) sp.Init();

        // Start calculation: thermalize
        varmc.Walk(int(therm*L*L),0);

        int stat_count=0;

        // Calculation
        for(size_t m=0;m<saves;++m){
            Timer::tic("main/ranwalk");
            varmc.Walk(L*L*N/saves,L*L);
            Timer::toc("main/ranwalk");
            if(qs.size()) spindense.save();
            heisen.save();
            sscorr.save();
            rej=varmc.Rejection();
            std::stringstream ahead;
            if(stat_count>=fm.StatPerSample()) stat_count=0;
            fm.DataAttribute("rej",rej);
            fm.DataAttribute("time",Timer::timer("main/ranwalk"));
            fm.DataAttribute("statistics",(stat_count+1)*N/saves);
#ifdef USEMPI
            int stop(0), mess(fm.message_save);
            if(m!=saves-1) stop=comm_rank;
            MPI_Send(&mess,1, MPI_INT,0,0,MPI_COMM_WORLD);
            MPI_Send(&stop,1,MPI_INT,0,0,MPI_COMM_WORLD);
#endif
            fm.Write();
            ++stat_count;
        }

        if(jas) delete jas;
#ifdef USEMPI
    } else {
        fm.MainLoop();
    }
#endif
    // Output
    Timer::toc("main");
    if(comm_rank==0){
        for(int r=0;r<comm_size;++r){
            std::cout<<"############################################"<<std::endl;
            std::cout<<"rank="<<r<<endl;
            if(r==0)
                std::cout<<Timer::report();
#ifdef USEMPI
            else {
                int len(0);
                MPI_Recv(&len,1,MPI_INT,r,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                char *mes=new char[len+1];
                MPI_Recv(mes,len+1,MPI_CHAR,r,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                cout<<mes<<endl;
                delete [] mes;
            }
#endif
            std::cout<<"############################################"<<std::endl;
        }
        std::cout<<"calc id: "<<fm.Prefix()<<std::endl;
    } else {
#ifdef USEMPI
        string rep=Timer::report();
        int len=rep.size();
        MPI_Send(&len,1,MPI_INT,0,1,MPI_COMM_WORLD);
        char* mes=new char[len+1];
        memcpy(mes,rep.c_str(),(len+1)*sizeof(char));
        MPI_Send(mes,len+1,MPI_CHAR,0,1,MPI_COMM_WORLD);
        delete [] mes;
#endif
    }
    delete [] Rs;
#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}

