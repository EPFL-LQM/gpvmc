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
    int seed=time(NULL)+100*comm_rank;
    RanGen::srand(seed);
    Timer::tic("main");

    // Parse input arguments
    ArgParse arg(argc,argv);
    size_t L=arg.i("L");
    double phi=arg.d("phi");
    double neel=arg.d("neel");
    double jastrow=arg.d("jastrow");
    size_t N=arg.i("N");
    size_t s=arg.i("s");
    size_t saves=arg.i("saves");
    int prefix=arg.i("prefix");
    string dir=arg.s("dir");
    size_t Q[2]={arg.i("qx"),arg.i("qy")};
    int therm=arg.i("therm");
    double phase_shift[2]={arg.d("phase_shift_x"),arg.d("phase_shift_y")};
    bool jas_stagmagn=arg.b("jas_onebodystag");
    bool jas_stag=arg.b("jas_twobodystag");
    double jr=arg.d("Jr");
    bool verbose=arg.b("verbose");

    // Setup calculation parameters
    FileManager fm(dir,prefix);
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
    fm.FileAttribute("longitudinal",1);
    if(jas_stag) fm.FileAttribute("twobodystagjastrow",1);
    if(jas_stagmagn) fm.FileAttribute("onebodystagjastrow",1);

    fm.StatPerSample()=saves/s;

#ifdef USEMPI
    if(comm_rank){
#endif
        // Setup calculation
        double rej=0;
        SpinState sp(L,L*L/2,L*L/2);
        StagFluxLongExciton wav(L,L,phi,neel,phase_shift,Q);
        Jastrow* jas=0;
        if(jastrow!=0){
            if(jas_stagmagn)
                jas=new StaggMagnJastrow(&sp,jastrow);
            else if(jas_stag)
                jas=new StagJastrow(&sp,jastrow);
        }
        Amplitude amp(&sp,&wav);
        while(amp.Amp()==0.0){
            sp.Init();
            amp.Init();
        }
        FullSpaceStepper step(&amp);
        ProjHeis seen(&step,&fm,jr);
        MetroMC varmc(&step,&fm);
        varmc.AddQuantity(&seen);

        // Start calculation: thermalize
        varmc.Walk(int(therm*L*L),0);
        
        int stat_count=0;

        // Calculation
        for(size_t m=0;m<saves;++m){
            Timer::tic("main/ranwalk");
            varmc.Walk(L*L*N/saves,L*L,!verbose);
            Timer::toc("main/ranwalk");
            rej=varmc.Rejection();
            seen.save();
            if(stat_count>=fm.StatPerSample()) stat_count=0;
            fm.DataAttribute("rej",rej);
            fm.DataAttribute("time",Timer::timer("randwalk"));
            fm.DataAttribute("statistics",(stat_count+1)*N/saves);
#ifdef USEMPI
            int stop(0), mess(fm.message_save);
            if(m!=saves-1) stop=comm_rank;
            MPI_Send(&mess,1,MPI_INT,0,0,MPI_COMM_WORLD);
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

