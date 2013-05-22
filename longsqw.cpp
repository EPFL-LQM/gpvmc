#ifdef USEMPI
#include <mpi.h>
#endif
#include <ctime>
#include <vector>
#include "MetroMC.h"
#include "SpinState.h"
#include "StagFluxLongExciton.h"
#include "FullSpaceStepper.h"
#include "ProjHeis.h"
#include "Amplitude.h"
#include "RanGen.h"
#include "FileManager.h"
#include "StaggMagnJastrow.h"
#include "StagJastrow.h"
#include "ArgParse.h"
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <hdf5.h>
#include <hdf5_hl.h>

using namespace std;

void setup_params(int argc,char* argv[],
                  size_t& L,
                  size_t& N,
                  size_t& s,
                  size_t& saves,
                  size_t* Q,
                  int& prefix,
                  int& therm,
                  int& verbose,
                  double& phi,
                  double& neel,
                  double& jastrow,
                  double* phase_shift,
                  double& jr,
                  bool& jas_stagmagn,
                  bool& jas_stag,
                  string& dir,
                  string& spinstste);

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
    int seed=time(NULL)+100*comm_rank;
    RanGen::srand(seed);
    Timer::tic("main");

    size_t L,N,s,saves,Q[2];
    int prefix,therm,verbose;
    double phi,neel,jastrow,phase_shift[2],jr;
    bool jas_stagmagn,jas_stag;
    string dir, spinstate;
    setup_params(argc,argv,
                 L,
                 N,
                 s,
                 saves,
                 Q,
                 prefix,
                 therm,
                 verbose,
                 phi,
                 neel,
                 jastrow,
                 phase_shift,
                 jr,
                 jas_stagmagn,
                 jas_stag,
                 dir,spinstate);

    // Setup calculation parameters
    FileManager fm(dir,prefix);
    MPI_Barrier(MPI_COMM_WORLD);
    if(comm_rank==0) std::cout<<Q[0]<<" "<<Q[1]<<std::endl;
    phi*=M_PI;
    fm.FileAttribute("L",L);
    fm.FileAttribute("N",N);
    fm.FileAttribute("neel",neel);
    fm.FileAttribute("phi",phi/M_PI);
    fm.FileAttribute("jastrow",jastrow);
    fm.FileAttribute("s",s);
    fm.FileAttribute("saves",saves);
    fm.FileAttribute("qx",Q[0]);
    fm.FileAttribute("qy",Q[1]);
    fm.FileAttribute("phasex",phase_shift[0]);
    fm.FileAttribute("phasey",phase_shift[1]);
    fm.FileAttribute("Jr",jr);
    fm.FileAttribute("channel","long");
    fm.FileAttribute("spinstate",spinstate);
    if(jas_stag) fm.FileAttribute("twobodystagjastrow",1);
    if(jas_stagmagn) fm.FileAttribute("onebodystagjastrow",1);

    fm.StatPerSample()=saves/s;

#ifdef USEMPI
    if(comm_rank){
#endif
        // Setup calculation
        double rej=0;
        SpinState sp(L,L*L/2,L*L/2);
        if(spinstate!=""){
            char* ist=new char[L*L];
#ifdef USEMPI
            if(comm_rank==1){
#endif
                hid_t ifile=H5Fopen(spinstate.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
#ifdef USEMPI
                H5LTread_dataset_char(ifile,"/rank-1",ist);
#else
                H5LTread_dataset_char(ifile,"/rank-0",ist);
#endif
                sp.Init(ist);
#ifdef USEMPI
                for(int r=2;r<comm_size;++r){
                    ostringstream rst;
                    rst<<"/rank-"<<r;
                    H5LTread_dataset_char(ifile,rst.str().c_str(),ist);
                    MPI_Send(ist,L*L,MPI_CHAR,r,0,MPI_COMM_WORLD);
                }
                H5Fclose(ifile);
            } else {
                MPI_Recv(ist,L*L,MPI_CHAR,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                sp.Init(ist);
            }
#endif
            delete [] ist;
        } else {
           sp.Init();
        } 
        StagFluxLongExciton wav(L,phi,neel,phase_shift,Q);
        wav.save(&fm);
        Jastrow* jas=0;
        if(jastrow!=0){
            if(jas_stagmagn)
                jas=new StaggMagnJastrow(&sp,jastrow);
            else if(jas_stag)
                jas=new StagJastrow(&sp,jastrow);
        }
        Amplitude amp(&sp,&wav);
        if(spinstate==""){
            amp.Init();
            while(amp.Amp()==0.0){
                sp.Init();
                amp.Init();
            }
        } else {
           amp.Init();
           if(amp.Amp()==0.0){
              cerr<<comm_rank<<"bad starting point!"<<endl;
           }
        } 
        FullSpaceStepper step(&amp);
        ProjHeis seen(&step,&fm,jr);
        MetroMC varmc(&step,&fm);
        varmc.AddQuantity(&seen);

        // Start calculation: thermalize
        if(spinstate==""){
            varmc.Walk(int(therm*L*L),0);
        }
        
        int stat_count=0;

        // Calculation
        for(size_t m=0;m<saves;++m){
            Timer::tic("main/randwalk");
            varmc.Walk(L*L*N/saves,L*L,!verbose,m);
            Timer::toc("main/randwalk");
            rej=varmc.Rejection();
            seen.save();
            if(stat_count>=fm.StatPerSample()) stat_count=0;
            fm.DataAttribute("rej",rej);
            fm.DataAttribute("time",Timer::timer("randwalk"));
            fm.DataAttribute("statistics",(stat_count+1)*N/saves);
#ifdef USEMPI
            int stop(0), mess(fm.message_save);
            if(m!=saves-1) stop=comm_rank;
            MPI_Send(&mess,1,MPI_INT,0,fm.message_comm,MPI_COMM_WORLD);
            MPI_Send(&stop,1,MPI_INT,0,fm.message_save,MPI_COMM_WORLD);
#endif
            fm.Write();
            ++stat_count;
        }
        sp.save(&fm);
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
        std::cout<<"output prefix="<<fm.Prefix()<<std::endl;
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
                  size_t* Q,
                  int& prefix,
                  int& therm,
                  int& verbose,
                  double& phi,
                  double& neel,
                  double& jastrow,
                  double* phase_shift,
                  double& jr,
                  bool& jas_stagmagn,
                  bool& jas_stag,
                  string& dir,
                  string& spinstate)
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
        Q[0]=arg.i("qx");Q[1]=arg.i("qy");
        therm=arg.i("therm");
        phase_shift[0]=arg.d("phase_shift_x");phase_shift[1]=arg.d("phase_shift_y");
        jas_stagmagn=arg.b("jas_onebodystag");
        jas_stag=arg.b("jas_twobodystag");
        jr=arg.d("Jr");
        verbose=arg.i("verbose");
        if(arg.HasArg("spinstate")){
            spinstate=arg.s("spinstate");
        } else {
            spinstate="";
        }
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
    // send string "spinstate":
    if(comm_rank==0)
        strlen=spinstate.size();
    MPI_Bcast(&strlen,1,MPI_INT,0,MPI_COMM_WORLD);
    char* spinstate_c_str=new char[strlen+1];
    if(comm_rank==0)
        memcpy(spinstate_c_str,spinstate.c_str(),(strlen+1));
    MPI_Bcast(spinstate_c_str,strlen+1,MPI_CHAR,0,MPI_COMM_WORLD);
    spinstate=string(spinstate_c_str);
    delete [] spinstate_c_str;
    MPI_Bcast(Q,2*sizeof(size_t),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&therm,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(phase_shift,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jas_stagmagn,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jas_stag,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Bcast(&jr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&verbose,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
}
