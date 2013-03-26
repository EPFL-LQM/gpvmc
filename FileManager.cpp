#ifdef USEMPI
#include <mpi.h>
#endif
#include <dirent.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <signal.h>
//#include <sys/sendfile.h>

#include "FileManager.h"
using namespace std;

FileManager::FileManager(const string& dir, const int& num)
    :m_dir(dir)
{
    if(num>=0){
        m_num=num;
    } else {
        int rank(0);
#ifdef USEMPI
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
        ostringstream dirstr;
        dirstr<<dir<<"/.runfile";
        if(rank==0){
            ofstream orunfile;
            try{
                ifstream irunfile;
                irunfile.open(dirstr.str().c_str());
                if(!irunfile.is_open()){
                    ifstream::failure e("Cannot open runfile");
                    throw e;
                }
                irunfile>>m_num;
                ++m_num;
                irunfile.close();
            } catch(ifstream::failure e){
                m_num=0;
            }
            try{
                ofstream orunfile;
                orunfile.open(dirstr.str().c_str());
                orunfile<<m_num<<endl;
                orunfile.close();
            } catch(ofstream::failure e){
                cerr<<"error writing runfile: "<<e.what()<<endl;
#ifdef EXCEPT
                throw e;
#else
                abort();
#endif
            }
        }
#ifdef USEMPI
        MPI_Bcast(&m_num,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
    }
}

#ifdef USEMPI
void FileManager::MainLoop()
{
    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(rank!=0) return;
    vector<MPI_Request> requests(size-1);
    vector<int> ready(size-1);
    vector<double> done(size-1);
    vector<double> tot(size-1);
    vector<int> num(size-1);
    for(int p=0;p<size-1;++p) ready[p]=p+1;
    bool loop=true;
    for(int p=1;p<size;++p)
        MPI_Irecv(&ready[p-1],1,MPI_INT,p,0,MPI_COMM_WORLD,&requests[p-1]);
    while(loop){
        MPI_Status status;
        int isready;
        MPI_Waitany(size-1,&requests[0],&isready,&status);
        int mess=ready[isready];
        isready+=1;
        if(mess==message_monitor){
            MPI_Recv(&done[isready-1],1,MPI_DOUBLE,isready,0,MPI_COMM_WORLD,&status);
            MPI_Recv(&tot[isready-1],1,MPI_DOUBLE,isready,0,MPI_COMM_WORLD,&status);
            MPI_Recv(&num[isready-1],1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
            ready[isready-1]=isready;
            Monitor(ready,done,tot,num);
        } else {
            MPI_Recv(&ready[isready-1],1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
            Write(isready);
        }
        loop=false;
        for(int p=1;p<size;++p)
            if(ready[p-1]) loop=true;
        if(loop){ 
            MPI_Irecv(&ready[isready-1],1,MPI_INT,isready,0,
                      MPI_COMM_WORLD,&requests[isready-1]);
        }
    }
}
#endif

void FileManager::Monitor(const vector<int>& ranks,
                          const vector<double>& percents,
                          const vector<double>& total_time,
                          const vector<int>& num_rep)
{
    // each ranks take 12 character long
    // get terminal width:
    size_t colwid=80;
    size_t nrperline=colwid/12;
    size_t nlines=ranks.size()/nrperline+1;
    size_t nsaves=m_fileattr["saves"];
    ostringstream ostr;
    for(size_t nl=0;nl<nlines;++nl){
        size_t r0=nl*nrperline;
        size_t re=min(ranks.size(),nrperline*(nl+1));
        if(nl==0){
            for(size_t r=r0;r<re;++r) ostr<<"____________";
        }
        ostr<<endl;
        for(size_t r=r0;r<re;++r){
            ostr.unsetf(ios_base::right);
            ostr.setf(ios_base::left);
            ostr<<"r-"<<setw(4)<<setfill(' ')<<ranks[r];
            ostr.unsetf(ios_base::left);
            ostr.setf(ios_base::right);
            ostr<<setw(3)<<setfill(' ')<<int((num_rep[r]+percents[r])/nsaves*100.0)<<"% |";
        }
        ostr<<endl;
        for(size_t r=r0;r<re;++r){
            int d=floor(total_time[r]*nsaves/3600/24);
            int h=floor(total_time[r]*nsaves/3600-d*24);
            int m=floor(total_time[r]*nsaves/60-d*24*60-h*60);
            int s=floor(total_time[r]*nsaves-d*24*3600-h*3600-m*60);
            ostr<<setw(2)<<d<<"-"<<setw(2)<<h<<":"<<setw(2)<<m<<":"<<setw(2)<<s<<"|";
        }
        ostr<<endl;
        for(size_t r=r0;r<re;++r) ostr<<"____________";
    }
    cout<<ostr.str()<<endl;
}

void FileManager::FileAttribute(string attr, double val)
{
    m_fileattr[attr]=val;
}

void FileManager::FileAttribute(string attr, string val)
{
    m_file_str_attr[attr]=val;
}

void FileManager::DataAttribute(string attr, double val)
{
    m_dataattr[attr]=val;
}

int& FileManager::StatPerSample()
{
    return m_stat_per_sample;
}

void FileManager::Write(int isready)
{
#ifdef USEMPI
    int rank(0),size(1);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Status status;
    if(rank==0){
        // Get number of quantities
        int Nq;
        MPI_Recv(&Nq,1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
        for(int nq=0;nq<Nq;++nq){
            // Get name of quantity
            char strin[128];
            int len;
            MPI_Recv(&len,1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
            MPI_Recv(strin,len,MPI_CHAR,isready,0,MPI_COMM_WORLD,&status);
            // Open corresponding file
            ostringstream fn;
            fn<<m_dir<<"/"<<m_num<<'-'<<strin<<".h5";
            struct stat sb;
            hid_t fout;
            if(stat(fn.str().c_str(),&sb)==-1){
                fout=H5Fcreate(fn.str().c_str(),
                               H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT);
                // set file attributes
                map<string,double>::iterator atit=m_fileattr.begin();
                while(atit!=m_fileattr.end()){
                    H5LTset_attribute_double(fout,"/",atit->first.c_str(),
                                             &(atit->second),1);
                    atit++;
                }
                map<string,string>::iterator stratit=m_file_str_attr.begin();
                while(stratit!=m_file_str_attr.end()){
                    H5LTset_attribute_string(fout,"/",stratit->first.c_str(),
                                             stratit->second.c_str());
                    stratit++;
                }
                H5LTset_attribute_string(fout,"/","type",strin);
            } else {
                // make backup of file in case process is killed while writing
                /*int src,dest;
                src=open(fn.str().c_str(),O_RDONLY,0);
                ostringstream fnbak;
                fnbak<<fn.str()<<".bak";
                dest=open(fnbak.str().c_str(),O_WRONLY | O_CREAT | O_TRUNC, 0644);
                struct stat sbb;
                fstat(src,&sbb);
                cout<<"copy file"<<endl;
                sendfile(dest,src,0,sbb.st_size);
                close(src);
                close(dest);*/
                // now open file for reading
                fout=H5Fopen(fn.str().c_str(),
                             H5F_ACC_RDWR,H5P_DEFAULT);
            }
            if(fout<0){
                ostringstream err;
                err<<"Could not create nor open file \""<<fn.str()<<"\".";
#ifdef EXCEPT
                throw(std::runtime_error(err.str()));
#else
                cerr<<err.str()<<endl;
                abort();
#endif
            }
            // get group for rank
            ostringstream gout;
            gout<<"/rank-"<<isready;
            // get current stat
            int statistics;
            MPI_Recv(&statistics,1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
            ostringstream dout;
            H5G_info_t info;
            if(H5Lexists(fout,gout.str().c_str(),H5P_DEFAULT)<=0){
                H5Gcreate2(fout,gout.str().c_str(),
                        H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
                dout<<gout.str()<<"/data-0";
            } else {
                H5Gget_info_by_name(fout,gout.str().c_str(),
                                    &info,H5P_DEFAULT);
                // get dataset number
                if(statistics==1)
                    dout<<gout.str()<<"/data-"<<info.nlinks;
                else
                    dout<<gout.str()<<"/data-"<<(info.nlinks-1);
            }
            // get data and writeo
            int dims[2];
            MPI_Recv(dims,2,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
            double *buff=new double[dims[0]*dims[1]];
            MPI_Recv(buff,dims[0]*dims[1],MPI_DOUBLE,isready,0,
                     MPI_COMM_WORLD,&status);
            hsize_t hdims[2]={dims[0],dims[1]};
            herr_t e;
            if(H5Lexists(fout,dout.str().c_str(),H5P_DEFAULT)<=0)
                e=H5LTmake_dataset_double(fout,dout.str().c_str(),
                                          2,hdims,buff);
            else{
                hid_t did=H5Dopen2(fout,dout.str().c_str(),H5P_DEFAULT);
                hid_t dt=H5Dget_type(did);
                e=H5Dwrite(did,dt,H5S_ALL,H5S_ALL,H5P_DEFAULT,buff);
                e=H5Tclose(dt);
                e=H5Dclose(did);
            }
            if(e<0){
                ostringstream err;
                err<<"Could not write dataset in file \""<<fn.str()<<"\".";
#ifdef EXCEPT
                throw(std::runtime_error(err.str()));
#else
                cerr<<err.str()<<endl;
                abort();
#endif
            }
            delete [] buff;
            // set data attributes
            int Na;
            double da;
            MPI_Recv(&Na,1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
            for(int a=0;a<Na;++a){
                MPI_Recv(&len,1,MPI_INT,isready,0,MPI_COMM_WORLD,&status);
                MPI_Recv(strin,len,MPI_CHAR,isready,0,MPI_COMM_WORLD,&status);
                MPI_Recv(&da,1,MPI_DOUBLE,isready,0,MPI_COMM_WORLD,&status);
                H5LTset_attribute_double(fout,dout.str().c_str(),strin,&da,1);
            }
            e=H5Fclose(fout);
        }
    } else {
        map<string,MatStream >::iterator it=m_streams.begin();
        int Nq=m_streams.size();
        MPI_Send(&Nq,1,MPI_INT,0,0,MPI_COMM_WORLD);
        while(it!=m_streams.end()){
            // send name of quantity.
            int len=it->first.size()+1;
            MPI_Send(&len,1,MPI_INT,0,0,MPI_COMM_WORLD);
            char *qname=new char[len];
            memcpy(qname,it->first.c_str(),len*sizeof(char));
            MPI_Send(qname,len,MPI_CHAR,0,0,
                     MPI_COMM_WORLD);
            delete [] qname;
            // send stat
            MPI_Send(&it->second.m_stat,1,MPI_INT,0,0,MPI_COMM_WORLD);
            if(it->second.m_stat==StatPerSample())
                it->second.m_stat=0;
            int dims[2]={it->second.m_nrow,it->second.m_ncol};
            MPI_Send(dims,2,MPI_INT,0,0,MPI_COMM_WORLD);
            MPI_Send(it->second.m_mat,dims[0]*dims[1],MPI_DOUBLE,0,0,
                     MPI_COMM_WORLD);
            // send data attributes
            map<string,double>::iterator ait=m_dataattr.begin();
            len=m_dataattr.size();
            MPI_Send(&len,1,MPI_INT,0,0,MPI_COMM_WORLD);
            while(ait!=m_dataattr.end()){
                len=ait->first.size()+1;
                MPI_Send(&len,1,MPI_INT,0,0,MPI_COMM_WORLD);
                char *aname=new char[len];
                memcpy(aname,ait->first.c_str(),len*sizeof(char));
                MPI_Send(aname,len,MPI_CHAR,0,0,MPI_COMM_WORLD);
                delete [] aname;
                MPI_Send(&(ait->second),1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
                ait++;
            }
            it++;
        }
    }
#else//!USEMPI
    map<string,MatStream >::iterator it=m_streams.begin();
    while(it!=m_streams.end()){
        ostringstream fn;
        fn<<m_dir<<"/"<<m_num<<"-"<<it->first<<".h5";
        struct stat sb;
        hid_t fout;
        if(stat(fn.str().c_str(),&sb)==-1){
            fout=H5Fcreate(fn.str().c_str(),
                           H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT);
            map<string,double>::iterator atit=m_fileattr.begin();
            while(atit!=m_fileattr.end()){
                H5LTset_attribute_double(fout,"/",atit->first.c_str(),
                                         &(atit->second),1);
                atit++;
            }
            map<string,string>::iterator stratit=m_file_str_attr.begin();
            while(stratit!=m_file_str_attr.end()){
                H5LTset_attribute_string(fout,"/",stratit->first.c_str(),
                                         stratit->second.c_str());
                stratit++;
            }
            H5LTset_attribute_string(fout,"/","type",it->first.c_str());
        } else {
            fout=H5Fopen(fn.str().c_str(),
                         H5F_ACC_RDWR,H5P_DEFAULT);
        }
        if(fout<0){
            ostringstream err;
            err<<"Could not create nor open file \""<<fn.str()<<"\".";
#ifdef EXCEPT
            throw(std::runtime_error(err.str()));
#else
            cerr<<err.str()<<endl;
            abort();
#endif
        }
        H5G_info_t info;
        ostringstream dout;
        dout<<"/rank-0";
        if(H5Lexists(fout,"/rank-0",H5P_DEFAULT)<=0){
            H5Gcreate2(fout,"/rank-0",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
            dout<<"/data-0";
        } else {
            H5Gget_info_by_name(fout,"/rank-0",&info,H5P_DEFAULT);
            if(it->second.m_stat==1)
                dout<<"/data-"<<info.nlinks;
            else
                dout<<"/data-"<<(info.nlinks-1);
        }
        if(it->second.m_stat==StatPerSample())
            it->second.m_stat=0;
        hsize_t hdims[2]={it->second.m_nrow, it->second.m_ncol};
        herr_t e;
        if(H5Lexists(fout,dout.str().c_str(),H5P_DEFAULT)<=0)
            e=H5LTmake_dataset_double(fout,dout.str().c_str(),2,hdims,
                                      it->second.m_mat);
        else {
            hid_t did=H5Dopen2(fout, dout.str().c_str(),H5P_DEFAULT);
            hid_t dt=H5Dget_type(did);
            e=H5Dwrite(did,dt,H5S_ALL,H5S_ALL,H5P_DEFAULT,it->second.m_mat);
            e=H5Tclose(dt);
            e=H5Dclose(did);
        }
        if(e<0){
            ostringstream err;
            err<<"Could not write dataset in file \""<<fn.str()<<"\".";
#ifdef EXCEPT
            throw(std::runtime_error(err.str()));
#else
            cerr<<err.str()<<endl;
            abort();
#endif
        }
        map<string,double>::iterator datit=m_dataattr.begin();
        while(datit!=m_dataattr.end()){
            H5LTset_attribute_double(fout,dout.str().c_str(),
                                     datit->first.c_str(),&datit->second,1);
            datit++;
        }
        e=H5Fclose(fout);
        it++;
    }
#endif//USEMPI
}

void FileManager::EmergencyClose(int signum)
{
    signal(signum,FileManager::EmergencyClose);
    H5close();
    exit(0);
}

MatStream& FileManager::FileStream(string basename)
{
    return m_streams[basename];
}
