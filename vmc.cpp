#ifdef USEMPI
#include <mpi.h>
#endif
#include <ctime>
#include <vector>
#include "MetroMC.h"
#include "SpinState.h"
#include "StagFluxTransExciton.h"
#include "StagFluxLongExciton.h"
#include "StagFluxGroundState.h"
#include "FullSpaceStepper.h"
#include "ProjHeis.h"
#include "StatSpinStruct.h"
#include "StagMagn.h"
#include "OverlapTrack.h"
#include "StagMagnTrack.h"
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
    domap["jastrow"]=0.0;
    domap["phase_shift_x"]=1.0;
    domap["phase_shift_y"]=1.0;
    domap["jr"]=0.0;
    domap["cutoff"]=0.0;
    bomap["jas_onebodystag"]=true;
    bomap["jas_twobodystag"]=false;
    bomap["track_stagmagn"]=false;
    bomap["track_overlap"]=false;
    bomap["meas_projheis"]=true;
    bomap["meas_stagmagn"]=true;
    bomap["meas_statspinstruct"]=true;
    stmap["dir"]=".";
    stmap["spinstate"]="";
    stmap["channel"]="groundstate";
    ArgParse arg(argc,argv);
    arg.SetupParams(bomap,simap,inmap,domap,stmap);
    if(!simap["meas_interv"])
        simap["meas_interv"]=pow(simap["L"],2);
    if(!comm_rank) cout<<"seed="<<inmap["seed"]<<endl;
    RanGen::srand(inmap["seed"]+100*comm_rank);
    // Setup calculation parameters
    FileManager fm(stmap["dir"],inmap["prefix"]);
    fm.Verbose()=inmap["verbose"];
    fm.MonitorTotal()=simap["samples"]*simap["samples_saves"];
    fm.StatPerSample()=simap["samples_saves"];
    domap["phi"]*=M_PI;
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

#ifdef USEMPI
    if(comm_rank){
#endif
        // Setup calculation
        double rej=0;
        size_t L=simap["L"], Q[2];
        double neel=domap["neel"], phi=domap["phi"], cutoff=domap["cutoff"];
        double phase_shift[2];
        phase_shift[0]=domap["phase_shift_x"];
        phase_shift[1]=domap["phase_shift_y"];
        Q[0]=simap["qx"];
        Q[1]=simap["qy"];
        SpinState* sp(0);
        if(stmap["channel"]=="groundstate" || stmap["channel"]=="long"){
            sp=new SpinState(L,L*L/2,L*L/2);
        } else if(stmap["channel"]=="trans"){
            sp=new SpinState(L,L*L/2+1,L*L/2-1);
        }
        if(stmap["spinstate"]!=""){
            char* ist=new char[L*L];
#ifdef USEMPI
            if(comm_rank==1){
#endif
                hid_t ifile=H5Fopen(stmap["spinstate"].c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
#ifdef USEMPI
                H5LTread_dataset_char(ifile,"/rank-1",ist);
#else
                H5LTread_dataset_char(ifile,"/rank-0",ist);
#endif
                sp->Init(ist);
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
                sp->Init(ist);
            }
#endif
            delete [] ist;
        } else {
            sp->Init();
        }
        double bdwd=2*sqrt(1+neel*neel);
        WaveFunction* wav(0);
        if(stmap["channel"]=="groundstate"){
            wav=new StagFluxGroundState(L,L,phi,neel,phase_shift);
        } else if(stmap["channel"]=="trans"){
            wav=new StagFluxTransExciton(L,L,phi,neel,phase_shift,Q,cutoff*bdwd);
        } else if(stmap["channel"]=="long"){
            wav=new StagFluxLongExciton(L,phi,neel,phase_shift,Q);
        }
        wav->save(&fm);
        Jastrow* jas=0;
        if(domap["jastrow"]!=0){
            if(bomap["jas_onebodystag"])
                jas=new StaggMagnJastrow(sp,domap["jastrow"]);
            else if(bomap["jas_twobodystag"])
                jas=new StagJastrow(sp,domap["jastrow"]);
        }
        Amplitude amp(sp,wav);
        if(stmap["spinstate"]==""){
            while(amp.Amp()==0.0){
                sp->Init();
                amp.Init();
            }
        } else {
            amp.Init();
        }
        FullSpaceStepper step(&amp);
        MetroMC varmc(&step,&fm);
        if(bomap["meas_projheis"]){
            ProjHeis* heisen=new ProjHeis(&step,&fm,domap["jr"]);
            varmc.AddQuantity(heisen);
        }
        if(stmap["channel"]=="groundstate"){
            if(bomap["meas_stagmagn"]){
                StagMagn* stagsz=new StagMagn(&step,&fm);
                varmc.AddQuantity(stagsz);
            }
            if(bomap["meas_statspinstruct"]){
                StatSpinStruct* stat=new StatSpinStruct(&step,&fm);
                varmc.AddQuantity(stat);
            }
        }
        if(bomap["track_stagmagn"]){
            StagMagnTrack* stt=new StagMagnTrack(&step,&fm);
            varmc.AddQuantity(stt);
        }
        if(bomap["track_overlap"]){
            OverlapTrack* ovt=new OverlapTrack(&step,&fm);
            varmc.AddQuantity(ovt);
        }

        // Start calculation: thermalize
        if(simap["therm"]){
            Timer::tic("main/thermalize");
            varmc.Walk(int(simap["therm"]*L*L),0);
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
#ifdef USEMPI
                int mess(fm.message_save);
                //cout<<"rank "<<comm_rank<<": sends message_save"<<endl;
                MPI_Send(&mess,1,MPI_INT,0,fm.message_comm,MPI_COMM_WORLD);
#endif
                fm.Write();
            }
        }
#ifdef USEMPI
        int mess=fm.message_loop, stop=1;
        //cout<<"rank "<<comm_rank<<": sends message_loop"<<endl;
        MPI_Send(&mess,1,MPI_INT,0,fm.message_comm,MPI_COMM_WORLD);
        MPI_Send(&stop,1,MPI_INT,0,fm.message_loop,MPI_COMM_WORLD);
#endif
        sp->save(&fm);
        if(jas) delete jas;
        delete wav;
        delete sp;
        for(size_t qu=0;qu<varmc.GetQuantities().size();++qu)
            delete varmc.GetQuantities()[qu];
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

