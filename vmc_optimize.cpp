#ifdef USEMPI
#include <mpi.h>
#endif
#include <ctime>
#include <vector>
#include "MetroMC.h"
#include "LatticeState.h"
#include "SquareLattice.h"
#include "StagFluxTransExciton.h"
#include "StagFluxLongExciton.h"
#include "StagFluxGroundState.h"
#include "SFpNpHxWaveFunction.h"
#include "SFpNpHxGroundState.h"
#include "SFpNpHxExciton.h"
#include "LatticeStepper.h"
#include "ProjHeis.h"
#include "StatSpinStruct.h"
#include "StagMagn.h"
#include "Magnetization.h"
#include "OverlapTrack.h"
#include "StagMagnTrack.h"
#include "Amplitude.h"
#include "RanGen.h"
#include "FileManager.h"
#include "ArgParse.h"
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <gsl/gsl_multimin.h>

using namespace std;

struct Params_s{
    map<string,bool> bomap;
    map<string,size_t> simap;
    map<string,int> inmap;
    map<string,double> domap;
    map<string,string> stmap;
};

void wait_and_see(Params_s * params)
{
    int count(-1);
    int comm_size(1);
    int comm_rank(0);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);

    map<string,bool> bomap=((Params_s*)params)->bomap;
    map<string,size_t> simap=((Params_s*)params)->simap;
    map<string,int> inmap=((Params_s*)params)->inmap;
    map<string,double> domap=((Params_s*)params)->domap;
    map<string,string> stmap=((Params_s*)params)->stmap;
    MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
    while(count>=0){
        double neel, phi,hx;
        MPI_Bcast(&phi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&neel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&hx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        domap["phi"]=phi;
        domap["neel"]=neel;
        domap["hx"]=hx;
        // Setup calculation parameters
        ostringstream ostr;
        ostr<<inmap["prefix"]<<"-"<<count;
        FileManager fm(stmap["dir"],ostr.str());
        fm.Verbose()=inmap["verbose"];
        fm.MonitorTotal()=simap["samples"]*simap["samples_saves"];
        fm.StatPerSample()=simap["samples_saves"];
        for(map<string,bool>::iterator it=bomap.begin();it!=bomap.end();++it)
            fm.FileAttribute(it->first,it->second);
        for(map<string,size_t>::iterator it=simap.begin();it!=simap.end();++it)
            fm.FileAttribute(it->first,it->second);
        for(map<string,int>::iterator it=inmap.begin();it!=inmap.end();++it)
            fm.FileAttribute(it->first,it->second);
        for(map<string,double>::iterator it=domap.begin();it!=domap.end();++it)
            fm.FileAttribute(it->first,it->second);
        for(map<string,string>::iterator it=stmap.begin();it!=stmap.end();++it)
            fm.FileAttribute(it->first,it->second);
        fm.FileAttribute("gitversion",GIT_SHA1);
        double mean_energy(0), std_energy(0);

        // Setup calculation
        double rej=0;
        size_t L=simap["L"];
        vector<double> phase_shift(2);
        phase_shift[0]=domap["phase_shift_x"];
        phase_shift[1]=domap["phase_shift_y"];
        SquareLattice slat(L,L);
        LatticeState* sp(0);
        sp=new LatticeState(&slat,{L*L},{2});
        WaveFunction* wav(0);
        wav=new SFpNpHxGroundState(L,L,phi,neel,hx,phase_shift);
        wav->Save(&fm);
        Amplitude amp(sp,wav);
        vector<vector<size_t> > pop;
        pop.push_back({L*L/2,L*L/2});
        while(amp.Amp()==0.0){
            sp->RanInit(pop);
            amp.Init();
        }
        LatticeStepper step(&amp);
        MetroMC varmc(&step,&fm);
        ProjHeis* heisen=new ProjHeis(&step,&fm,&slat,domap["Bx"]);
        varmc.AddQuantity(heisen);
        if(bomap["meas_magnetization"] && (domap["Bx"]!=0 || !bomap["Sztot_conserved"])){
            Magnetization* magn=new Magnetization(&step,&fm);
            varmc.AddQuantity(magn);
        }
        // Start calculation: thermalize
        if(simap["therm"]){
            Timer::tic("main/thermalize");
            varmc.Walk(size_t(simap["therm"]*L*L),0);
            Timer::toc("main/thermalize");
        }
        fm.MonitorTotal()=simap["samples"]*simap["samples_saves"];         
        // Calculation
        for(size_t sample=0;sample<simap["samples"];++sample){
            for(size_t m=0;m<simap["samples_saves"];++m){
                fm.MonitorCompletion()=double(sample*simap["samples_saves"]+m)/(simap["samples"]*simap["samples_saves"]);
                Timer::tic("main/ranwalk");
                varmc.Walk(simap["meas_interv"]*simap["samples_saves_stat"],simap["meas_interv"]);
                Timer::toc("main/ranwalk");
                rej=varmc.Rejection();
                for(size_t qu=0;qu<varmc.GetQuantities().size();++qu)
                    varmc.GetQuantities()[qu]->save();
                fm.DataAttribute("rej",rej);
                fm.DataAttribute("time",Timer::timer("randwalk"));
                fm.DataAttribute("statistics",(m+1)*simap["samples_saves_stat"]);
                int mess(fm.message_save);
                //cout<<"rank "<<comm_rank<<": sends message_save"<<endl;
                MPI_Send(&mess,1,MPI_INT,0,fm.message_comm,MPI_COMM_WORLD);
                fm.Write();
            }
        }
        int mess=fm.message_loop, stop=1;
        //cout<<"rank "<<comm_rank<<": sends message_loop"<<endl;
        MPI_Send(&mess,1,MPI_INT,0,fm.message_comm,MPI_COMM_WORLD);
        MPI_Send(&stop,1,MPI_INT,0,fm.message_loop,MPI_COMM_WORLD);
        //sp->save(&fm);
        delete wav;
        delete sp;
        double en=real(heisen->Val(0,0))/heisen->GetNmeas();
        MPI_Gather(&en,1,MPI_DOUBLE,NULL,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        for(size_t qu=0;qu<varmc.GetQuantities().size();++qu)
            delete varmc.GetQuantities()[qu];
        MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
    }
}

double var_energy(const gsl_vector * x, void * params)
{
    static int count(0);
    int comm_size(1);
    int comm_rank(0);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);

    map<string,bool> bomap=((Params_s*)params)->bomap;
    map<string,size_t> simap=((Params_s*)params)->simap;
    map<string,int> inmap=((Params_s*)params)->inmap;
    map<string,double> domap=((Params_s*)params)->domap;
    map<string,string> stmap=((Params_s*)params)->stmap;

    //ask waiting processes to start calculating
    MPI_Bcast(&count,1,MPI_INT,0,MPI_COMM_WORLD);
    double phi,neel,hx;
    phi=gsl_vector_get(x,0)/4.0*M_PI;
    neel=gsl_vector_get(x,1)*2;
    hx=gsl_vector_get(x,2);
    //send variational parameters value
    MPI_Bcast(&phi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);//phi
    MPI_Bcast(&neel,1,MPI_DOUBLE,0,MPI_COMM_WORLD);//neel
    MPI_Bcast(&hx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);//hx
    // Setup calculation parameters
    ostringstream ostr;
    ostr<<inmap["prefix"]<<"-"<<count;
    FileManager fm(stmap["dir"],ostr.str());
    fm.Verbose()=inmap["verbose"];
    fm.MonitorTotal()=simap["samples"]*simap["samples_saves"];
    fm.StatPerSample()=simap["samples_saves"];
    for(map<string,bool>::iterator it=bomap.begin();it!=bomap.end();++it)
        fm.FileAttribute(it->first,it->second);
    for(map<string,size_t>::iterator it=simap.begin();it!=simap.end();++it)
        fm.FileAttribute(it->first,it->second);
    for(map<string,int>::iterator it=inmap.begin();it!=inmap.end();++it)
        fm.FileAttribute(it->first,it->second);
    for(map<string,double>::iterator it=domap.begin();it!=domap.end();++it)
        fm.FileAttribute(it->first,it->second);
    for(map<string,string>::iterator it=stmap.begin();it!=stmap.end();++it)
        fm.FileAttribute(it->first,it->second);
    fm.FileAttribute("gitversion",GIT_SHA1);
    double mean_energy(0), std_energy(0);

    fm.MainLoop();
    double dummy(0);
    vector<double> energies(comm_size,0);
    MPI_Gather(&dummy,1,MPI_DOUBLE,&energies[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(size_t n=1;n<comm_size;++n){
        mean_energy+=energies[n];
    }
    mean_energy/=(comm_size-1.0);
    for(size_t n=1;n<comm_size;++n){
        std_energy+=pow(energies[n]-mean_energy,2);
    }
    std_energy=sqrt(std_energy/(comm_size-2.0));
    ++count;
    cout<<"evaluation: "<<mean_energy<<" +/- "<<std_energy<<endl;
    return mean_energy;
}

int main(int argc, char* argv[])
{
    // Initialize
    MPI_Init(&argc,&argv);
    int comm_size,comm_rank;
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
    //int wait=0;
    //if(comm_rank!=0){
    //    cout<<getpid()<<endl;
    //    wait=1;
    //    while(wait)
    //        sleep(5);
    //}
    MPI_Barrier(MPI_COMM_WORLD);
    signal(SIGTERM,FileManager::EmergencyClose);
    Timer::tic("main");

    // calculation parameters:
    map<string,bool> bomap;
    map<string,size_t> simap;
    map<string,int> inmap;
    map<string,double> domap;
    map<string,string> stmap;
    simap["L"]=8;
    simap["samples"]=10;
    simap["samples_saves"]=1;
    simap["samples_saves_stat"]=100;
    simap["qx"]=0;
    simap["qy"]=0;
    simap["meas_interv"]=0;
    inmap["prefix"]=-1;
    inmap["therm"]=100;
    inmap["verbose"]=1;
    inmap["seed"]=time(NULL);
    domap["phi"]=0.085;
    domap["neel"]=0.0;
    domap["hx"]=0.0;
    domap["step_phi"]=0.05;
    domap["step_neel"]=0.025;
    domap["step_hx"]=0.05;
    domap["phase_shift_x"]=1.0;
    domap["phase_shift_y"]=1.0;
    domap["jr"]=0.0;
    domap["Bx"]=0.0;
    domap["tolerance"]=1e-4;
    bomap["track_stagmagn"]=false;
    bomap["track_overlap"]=false;
    bomap["meas_projheis"]=true;
    bomap["meas_stagmagn"]=true;
    bomap["meas_statspinstruct"]=true;
    bomap["meas_magnetization"]=false;
    bomap["Sztot_conserved"]=true;
    stmap["dir"]=".";
    stmap["spinstate"]="";
    stmap["channel"]="groundstate";
    ArgParse arg(argc,argv);
    arg.SetupParams(bomap,simap,inmap,domap,stmap);
    if(!simap["meas_interv"])
        simap["meas_interv"]=pow(simap["L"],2);
    if(!comm_rank) cout<<"seed="<<inmap["seed"]<<endl;
    RanGen::srand(inmap["seed"]+100*comm_rank);
    Params_s params;
    params.bomap=bomap;
    params.simap=simap;
    params.inmap=inmap;
    params.domap=domap;
    params.stmap=stmap;
    if(comm_rank){
        wait_and_see(&params);
    } else {
        const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer *s=NULL;
        gsl_vector *x,*ss;
        gsl_multimin_function minexc_func;
        size_t iter=0;
        int status;
        double size;
        x=gsl_vector_alloc(3);
        gsl_vector_set(x,0,domap["phi"]*4);
        gsl_vector_set(x,1,domap["neel"]/2);
        gsl_vector_set(x,2,domap["hx"]);// scale params such that they have approximate same scale
        ss=gsl_vector_alloc(3);
        gsl_vector_set(ss,0,domap["step_phi"]*4);
        gsl_vector_set(ss,1,domap["step_neel"]/2);
        gsl_vector_set(ss,2,domap["step_hx"]);
        minexc_func.n=3;
        minexc_func.f=var_energy;
        minexc_func.params=&params;
        s=gsl_multimin_fminimizer_alloc(T,3);
        gsl_multimin_fminimizer_set(s,&minexc_func,x,ss);
        do {
            iter++;
            status=gsl_multimin_fminimizer_iterate(s);
            if(status)
                break;
            size=gsl_multimin_fminimizer_size(s);
            status=gsl_multimin_test_size(size,domap["tolerance"]);
            if (status==GSL_SUCCESS){
                cout<<"converged to minimum at"<<endl;
            }
            cout<<iter<<" "<<gsl_vector_get(s->x,0)/4.0<<" "<<gsl_vector_get(s->x,1)*2<<" "<<gsl_vector_get(s->x,2)<<" = "<<s->fval<<" size = "<<size<<endl;
        } while(status==GSL_CONTINUE && iter<100);
        gsl_vector_free(x);
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free(s);
        int stop(-1);
        MPI_Bcast(&stop,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}

