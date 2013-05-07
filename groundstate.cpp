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
#include "StatSpinStruct.h"
#include "SpinDensity.h"
#include "FileManager.h"
#include "ArgParse.h"
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

using namespace std;

void setup_params(int argc,char* argv[],
                  size_t& L,
                  size_t& N,
                  size_t& s,
                  size_t& saves,
                  int& prefix,
                  int& therm,
                  int& verbose,
                  double& phi,
                  double& neel,
                  double& jastrow,
                  double* phase_shift,
                  double& jr,
                  double& cutoff,
                  bool& jas_stagmagn,
                  bool& jas_stag,
                  string& dir);

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

    // calculation parameters:
    size_t L,N,s,saves;
    int prefix,therm,verbose;
    double phi,neel,jastrow,phase_shift[2],jr,cutoff;
    bool jas_stagmagn,jas_stag;
    string dir;
    setup_params(argc,argv,
                 L,
                 N,
                 s,
                 saves,
                 prefix,
                 therm,
                 verbose,
                 phi,
                 neel,
                 jastrow,
                 phase_shift,
                 jr,
                 cutoff,
                 jas_stagmagn,
                 jas_stag,
                 dir);

    // Setup calculation parameters
    FileManager fm(dir,prefix);
    phi*=M_PI;
    double rej(0);
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
        wav.save(&fm);
        Amplitude amp(&sp,&wav);
        FullSpaceStepper step(&amp);
        ProjHeis heisen(&step,&fm,jr);
        StatSpinStruct stat(&step,&fm);
        MetroMC varmc(&step,&fm);
        varmc.AddQuantity(&heisen);
        varmc.AddQuantity(&stat);
        while(amp.Amp()==0){
            sp.Init();
            amp.Init();
        }

        // Start calculation: thermalize
        varmc.Walk(int(therm*L*L),0);

        int stat_count=0;

        // Calculation
        for(size_t m=0;m<saves;++m){
            Timer::tic("main/ranwalk");
            varmc.Walk(L*L*N/saves,L*L,!verbose,m);
            Timer::toc("main/ranwalk");
            heisen.save();
            stat.save();
            rej=varmc.Rejection();
            std::stringstream ahead;
            if(stat_count>=fm.StatPerSample()) stat_count=0;
            fm.DataAttribute("rej",rej);
            fm.DataAttribute("time",Timer::timer("main/ranwalk"));
            fm.DataAttribute("statistics",(stat_count+1)*N/saves);
#ifdef USEMPI
            int stop(0), mess(fm.message_save);
            if(m!=saves-1) stop=comm_rank;
            MPI_Send(&mess,1, MPI_INT,0,fm.message_comm,MPI_COMM_WORLD);
            MPI_Send(&stop,1,MPI_INT,0,fm.message_save,MPI_COMM_WORLD);
#endif
            fm.Write();
            ++stat_count;
        }

        if(jas) delete jas;
#ifdef USEMPI
    } else {
        fm.MainLoop(verbose);
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
#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}

void setup_params(int argc,char* argv[],
                  size_t& L,
                  size_t& N,
                  size_t& s,
                  size_t& saves,
                  int& prefix,
                  int& therm,
                  int& verbose,
                  double& phi,
                  double& neel,
                  double& jastrow,
                  double* phase_shift,
                  double& jr,
                  double& cutoff,
                  bool& jas_stagmagn,
                  bool& jas_stag,
                  string& dir)
{
    int comm_rank(0);
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
#endif
    // Parse input arguments
    if(comm_rank==0){
        ArgParse arg(argc,argv);
        L=arg.i("L");
        phi=arg.d("phi");
        neel=arg.d("neel");
        jastrow=arg.d("jastrow");
        N=arg.i("N");
        s=arg.i("s");
        saves=arg.i("saves");
        prefix=arg.i("prefix");
        dir=arg.s("dir");
        therm=arg.i("therm");
        phase_shift[0]=arg.d("phase_shift_x");phase_shift[1]=arg.d("phase_shift_y");
        jas_stagmagn=arg.b("jas_onebodystag");
        jas_stag=arg.b("jas_twobodystag");
        jr=arg.d("Jr");
        verbose=arg.i("verbose");
        cutoff=arg.d("cutoff");
    }
#ifdef USEMPI
    MPI_Bcast(&L,sizeof(size_t),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&phi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&neel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jastrow,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&N,sizeof(size_t),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&s,sizeof(size_t),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&saves,sizeof(saves),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&prefix,1,MPI_INT,0,MPI_COMM_WORLD);
    // send string "dir":
    int strlen;
    if(comm_rank==0)
        strlen=dir.size();
    MPI_Bcast(&strlen,1,MPI_INT,0,MPI_COMM_WORLD);
    char* dir_c_str=new char[strlen+1];
    if(comm_rank==0)
        memcpy(dir_c_str,dir.c_str(),(strlen+1));
    MPI_Bcast(dir_c_str,strlen+1,MPI_CHAR,0,MPI_COMM_WORLD);
    dir=string(dir_c_str);
    delete [] dir_c_str;
    MPI_Bcast(&therm,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(phase_shift,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jas_stagmagn,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jas_stag,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&verbose,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&cutoff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
}
