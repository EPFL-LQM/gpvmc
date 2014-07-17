#include "ArgParse.h"
#include <stdexcept>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif
#include <getopt.h>

using namespace std;

void print_help(const map<std::string,double>& double_opts,
                const map<std::string,int>& int_opts,
                const map<std::string,size_t>& uint_opts,
                const map<std::string,bool>& bool_opts,
                const map<std::string,std::string>& string_opts)
{
    cout<<"Available options:"<<endl;
    map<string,double>::const_iterator doiter=double_opts.begin();
    while(doiter!=double_opts.end()){
        cout<<"--"<<doiter->first<<"=<double type> (default "<<doiter->second<<")"<<endl;
        doiter++;
    }
    map<string,int>::const_iterator initer=int_opts.begin();
    while(initer!=int_opts.end()){
        cout<<"--"<<initer->first<<"=<int type> (default "<<initer->second<<")"<<endl;
        initer++;
    }
    map<string,size_t>::const_iterator siiter=uint_opts.begin();
    while(siiter!=uint_opts.end()){
        cout<<"--"<<siiter->first<<"=<size_t type> (default "<<siiter->second<<")"<<endl;
        siiter++;
    }
    map<string,bool>::const_iterator boiter=bool_opts.begin();
    while(boiter!=bool_opts.end()){
        cout<<"--"<<boiter->first<<endl;
        boiter++;
    }
    map<string,string>::const_iterator stiter=string_opts.begin();
    while(stiter!=string_opts.end()){
        cout<<"--"<<stiter->first<<"=<string type> (default \""<<stiter->second<<"\")"<<endl;
        stiter++;
    }
}

int ArgParse(int argc, char* argv[],
             map<std::string,double>& double_opts,
             map<std::string,int>& int_opts,
             map<std::string,size_t>& uint_opts,
             map<std::string,bool>& bool_opts,
             map<std::string,std::string>& string_opts)
{
    int rank=0;
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
    vector<string> keys,types;
    map<string,double>::iterator doiter=double_opts.begin();
    while(doiter!=double_opts.end()){
        keys.push_back(doiter->first);
        types.push_back("double");
        doiter++;
    }
    map<string,int>::iterator iiter=int_opts.begin();
    while(iiter!=int_opts.end()){
        keys.push_back(iiter->first);
        types.push_back("int");
        iiter++;
    }
    map<string,size_t>::iterator uiiter=uint_opts.begin();
    while(uiiter!=uint_opts.end()){
        keys.push_back(uiiter->first);
        types.push_back("uint");
        uiiter++;
    }
    map<string,bool>::iterator boiter=bool_opts.begin();
    while(boiter!=bool_opts.end()){
        keys.push_back(boiter->first);
        types.push_back("bool");
        boiter++;
    }
    map<string,string>::iterator striter=string_opts.begin();
    while(striter!=string_opts.end()){
        keys.push_back(striter->first);
        types.push_back("string");
        striter++;
    }
    bool help=false;
    if(rank==0){
        vector<struct option> options(keys.size()+1);
        for(size_t k=0;k<keys.size();++k){
            options[k].name=keys[k].c_str();
            if(types[k]=="bool")
                options[k].has_arg=no_argument;
            else
                options[k].has_arg=required_argument;
            options[k].flag=NULL;
            options[k].val=k;
        }
        options.back().name="help";
        options.back().has_arg=no_argument;
        options.back().flag=NULL;
        options.back().val=keys.size();

        struct option sigil;
        options.push_back(sigil);

        options.back().name=NULL;
        options.back().has_arg=0;
        options.back().flag=NULL;
        options.back().val=0;
        
        if(argc==1){
            print_help(double_opts,int_opts,uint_opts,bool_opts,string_opts);
            help=true;
        }
        int optidx;
        while(getopt_long(argc,argv,"",&options[0],&optidx)!=-1){
            if(optidx==keys.size()){
                print_help(double_opts,int_opts,uint_opts,bool_opts,string_opts);
                help=true;
            } else {
                if(types[optidx]=="double")
                    double_opts[keys[optidx]]=stod(optarg);
                if(types[optidx]=="int")
                    int_opts[keys[optidx]]=stoi(optarg);
                if(types[optidx]=="uint")
                    uint_opts[keys[optidx]]=stoul(optarg);
                if(types[optidx]=="string")
                    string_opts[keys[optidx]]=optarg;
                if(types[optidx]=="bool")
                    bool_opts[keys[optidx]]=true;
            }
        }
    }
#ifdef USEMPI
    MPI_Bcast(&help,sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
    for(size_t k=0;k<keys.size();++k){
        if(types[k]=="double")
            MPI_Bcast(&double_opts[keys[k]],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(types[k]=="int")
            MPI_Bcast(&int_opts[keys[k]],1,MPI_INT,0,MPI_COMM_WORLD);
        if(types[k]=="uint")
            MPI_Bcast(&uint_opts[keys[k]],sizeof(size_t),MPI_BYTE,0,MPI_COMM_WORLD);
        if(types[k]=="bool")
            MPI_Bcast(&bool_opts[keys[k]],sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
        if(types[k]=="string"){
            int strlen=string_opts[keys[k]].size()+1;
            MPI_Bcast(&strlen,1,MPI_INT,0,MPI_COMM_WORLD);
            char* c_str=new char[strlen];
            if(rank==0)
                memcpy(c_str,string_opts[keys[k]].c_str(),strlen);
            MPI_Bcast(c_str,strlen,MPI_CHAR,0,MPI_COMM_WORLD);
            string_opts[keys[k]]=string(c_str);
        }
    }
#endif//USEMPI
    return help;
}
