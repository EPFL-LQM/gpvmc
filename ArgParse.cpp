#include "ArgParse.h"
#include <stdexcept>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

ArgParse::ArgParse(int argc, char* argv[])
{
    for(int a=1;a<argc;++a){
        if(argv[a][0]=='-'){
            string in(argv[a]),trig,val;
            trig=in.substr(1,in.find('=')-1);
            if(in.find('=')==string::npos) {
                val='1';//true
            } else {
                val=in.substr(in.find('=')+1);
            }
            m_streams[trig]=val;
        } else {
            ifstream fin(argv[a]);
            string in,trig,val;
            while(getline(fin,in)) {
                if(in.size()==0 || in[0]=='#') continue;
                trig=in.substr(0,in.find('='));
                if(in.find('=')==string::npos) {
                    val='1';
                } else {
                    val=in.substr(in.find('=')+1);
                }
                m_streams[trig]=val;
            }
        }
    }
}

bool ArgParse::HasArg(const string& str) const
{
    return m_streams.find(str)!=m_streams.end();
}    

double ArgParse::d(const string& str) const
{
    double out;
    try {
        out=stod(m_streams.at(str));
        return out;
    } catch(out_of_range err) {
        cerr<<"No member \""<<str<<"\" found in arguments"<<endl;
        throw err;
    }
}

int ArgParse::i(const string& str) const
{
    int out;
    try {
        out=stoi(m_streams.at(str));
        return out;
    } catch (out_of_range err) {
        cerr<<"No member \""<<str<<"\" found in arguments"<<endl;
        throw err;
    }
}

bool ArgParse::b(const string& str) const
{
    bool out;
    try {
        out=stoi(m_streams.at(str));
        return out;
    } catch (out_of_range err) {
        cerr<<"No member \""<<str<<"\" found in arguments"<<endl;
        throw err;
    }
}
 string ArgParse::s(const string& str) const
{
    string out;
    try {
        out=m_streams.at(str);
        return out;
    } catch (out_of_range err) {
        cerr<<"No member \""<<str<<"\" found in arguments"<<endl;
        throw err;
    }
}

void ArgParse::SetupParams(std::map<std::string,bool>& map_bool,
                           std::map<std::string,size_t>& map_size_t,
                           std::map<std::string,int>& map_int,
                           std::map<std::string,double>& map_double,
                           std::map<std::string,std::string>& map_string)
{
    int comm_rank(0);
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
#endif
    if(comm_rank==0){
        map<string,bool>::iterator itb=map_bool.begin();
        while(itb!=map_bool.end()){
            if(HasArg(itb->first))
                itb->second=b(itb->first);
            itb++;
        }
        map<string,size_t>::iterator its=map_size_t.begin();
        while(its!=map_size_t.end()){
            if(HasArg(its->first))
                its->second=i(its->first);
            its++;
        }
        map<string,int>::iterator iti=map_int.begin();
        while(iti!=map_int.end()){
            if(HasArg(iti->first))
                iti->second=i(iti->first);
            iti++;
        }
        map<string,double>::iterator itd=map_double.begin();
        while(itd!=map_double.end()){
            if(HasArg(itd->first))
                itd->second=d(itd->first);
            itd++;
        }
        map<string,string>::iterator itst=map_string.begin();
        while(itst!=map_string.end()){
            if(HasArg(itst->first))
                itst->second=s(itst->first);
            itst++;
        }
    }
#ifdef USEMPI
    map<string,bool>::iterator itb=map_bool.begin();
    while(itb!=map_bool.end()){
        MPI_Bcast(&(itb->second),sizeof(bool),MPI_BYTE,0,MPI_COMM_WORLD);
        itb++;
    }
    map<string,size_t>::iterator its=map_size_t.begin();
    while(its!=map_size_t.end()){
        MPI_Bcast(&(its->second),sizeof(size_t),MPI_BYTE,0,MPI_COMM_WORLD);
        its++;
    }
    map<string,int>::iterator iti=map_int.begin();
    while(iti!=map_int.end()){
        MPI_Bcast(&(iti->second),1,MPI_INT,0,MPI_COMM_WORLD);
        iti++;
    }
    map<string,double>::iterator itd=map_double.begin();
    while(itd!=map_double.end()){
        MPI_Bcast(&(itd->second),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        itd++;
    }
    map<string,string>::iterator itst=map_string.begin();
    while(itst!=map_string.end()){
        int strlen;
        if(comm_rank==0)
            strlen=itst->second.size();
        MPI_Bcast(&strlen,1,MPI_INT,0,MPI_COMM_WORLD);
        char* c_str=new char[strlen+1];
        if(comm_rank==0)
            memcpy(c_str,itst->second.c_str(),strlen+1);
        MPI_Bcast(c_str,strlen+1,MPI_CHAR,0,MPI_COMM_WORLD);
        itst->second=string(c_str);
        delete [] c_str;
        itst++;
    }
#endif//USEMPI
}
