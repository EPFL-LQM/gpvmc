#ifdef USEMPI
#include <mpi.h>
#endif
#include <dirent.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <limits>
#include <iomanip>
#include <cmath>
#include <unistd.h>
#include <signal.h>
//#include <sys/sendfile.h>

#include "FileManager.h"
using namespace std;

FileManager::FileManager(const string& dir, const int& num)
    :m_dir(dir), m_compl(0), m_verbose(1), m_total(1)
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
    vector<MPI_Request> requests(size);
    vector<int> finished(size,0);
    vector<int> mess(size);
    vector<double> done(size-1);
    vector<double> tot(size-1);
    vector<int> num(size-1);
    bool loop=true;
    for(int p=1;p<size;++p)
        MPI_Irecv(&mess[p],1,MPI_INT,p,message_comm,MPI_COMM_WORLD,&requests[p]);
    while(loop){
        MPI_Status status;
        int isready;
        MPI_Waitany(size-1,&requests[1],&isready,&status);
        isready+=1;
        if(mess[isready]==message_monitor){
            //cout<<"rank 0 recieved message_monitor from"<<isready<<endl;
            Monitor(isready);
        } else if(mess[isready]==message_save){
            //cout<<"rank 0 recieved message_save from"<<isready<<endl;
            Write(isready);
        } else if(mess[isready]==message_loop){
            //cout<<"rank 0 received message_loop from"<<isready<<endl;
            MPI_Recv(&finished[isready],1,MPI_INT,isready,message_loop,MPI_COMM_WORLD,&status);
            loop=false;
            for(int p=1;p<size;++p)
                if(!finished[p]) loop=true;
        }
        if(!finished[isready]){ 
            MPI_Irecv(&mess[isready],1,MPI_INT,isready,message_comm,
                      MPI_COMM_WORLD,&requests[isready]);
        }
    }
}
#endif

void FileManager::Monitor(int rank,double done, double tottime)
{
    int comm_rank=0;
    int comm_size=1;
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
#endif
    static vector<double> done_vec(comm_size,0), time_vec(comm_size,0), compl_vec(comm_size,0);
#ifdef USEMPI
    if(comm_rank!=0){
        MPI_Send(&done,1,MPI_DOUBLE,0,message_monitor,MPI_COMM_WORLD);
        MPI_Send(&tottime,1,MPI_DOUBLE,0,message_monitor,MPI_COMM_WORLD);
        MPI_Send(&m_compl,1,MPI_DOUBLE,0,message_monitor,MPI_COMM_WORLD);
    } else {
        MPI_Recv(&done_vec[rank],1,MPI_DOUBLE,rank,message_monitor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&time_vec[rank],1,MPI_DOUBLE,rank,message_monitor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(&compl_vec[rank],1,MPI_DOUBLE,rank,message_monitor,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
#else
    done_vec[0]=done;
    time_vec[0]=tottime;
    compl_vec[0]=m_compl;
#endif
    if(comm_rank==0){
        if(m_verbose>1){
            ostringstream ostr;
            // each ranks take 12 character long
            // get terminal width:
            size_t colwid=80;
            size_t nrperline=colwid/12;
            size_t nlines=(comm_size-1)/nrperline+1;
            for(size_t nl=0;nl<nlines;++nl){
                size_t r0=nl*nrperline+(comm_size>0);
                size_t re=min(size_t(comm_size),nrperline*(nl+1));
                if(nl==0){
                    for(size_t r=r0;r<re;++r) ostr<<"____________";
                }
                ostr<<endl;
                for(size_t r=r0;r<re;++r){
                    ostr.unsetf(ios_base::right);
                    ostr.setf(ios_base::left);
                    ostr<<"r-"<<setw(4)<<setfill(' ')<<r+1;
                    ostr.unsetf(ios_base::left);
                    ostr.setf(ios_base::right);
                    ostr<<setw(3)<<setfill(' ')<<100*(compl_vec[r]+done_vec[r]/m_total)<<"% |";
                }
                ostr<<endl;
                for(size_t r=r0;r<re;++r){
                    if(time_vec[r]>0){
                        int d=floor(time_vec[r]*m_total/3600/24);
                        int h=floor(time_vec[r]*m_total/3600-d*24);
                        int m=floor(time_vec[r]*m_total/60-d*24*60-h*60);
                        int s=floor(time_vec[r]*m_total-d*24*3600-h*3600-m*60);
                        ostr<<setw(2)<<d<<"-"<<setw(2)<<h<<":"<<setw(2)<<m<<":"<<setw(2)<<s<<"|";
                    } else {
                        ostr<<"    N/A    |";
                    }
                }
                ostr<<endl;
                for(size_t r=r0;r<re;++r) ostr<<"____________";
            }
            cout<<ostr.str()<<endl;
        } else if(m_verbose==1){
            static int mean_adv=0;
            double new_mav=0;
            double minav=1;
            double maxav=0;
            double tav=-1,tmax=-1,tmin=numeric_limits<double>::infinity();
            int far(0),slr(0);
            int nnna=0;
            for(size_t r=(comm_size>1);r<size_t(comm_size);++r){
                double av=compl_vec[r]+done_vec[r]/m_total;
                new_mav+=av;
                if(time_vec[r]>=0){
                    tav+=time_vec[r]*m_total;
                    if(av<minav) minav=av;
                    if(av>maxav) maxav=av;
                    if(time_vec[r]*m_total<tmin){
                        tmin=time_vec[r]*m_total;
                        far=r;
                    }
                    if(time_vec[r]*m_total>tmax){
                        tmax=time_vec[r]*m_total;
                        slr=r;
                    }
                    ++nnna;
                }
            }
            if(comm_size>1)
                new_mav/=(comm_size-1);
            tav/=nnna;
            if(floor(new_mav*100)!=mean_adv){
                mean_adv=floor(new_mav*100);
                int dav=floor(tav/3600/24);
                int hav=floor(tav/3600-24*dav);
                int mav=floor(tav/60-24*60*dav-60*hav);
                int sav=floor(tav-dav*24*3600-hav*3600-mav*60);
                int dmin=floor(tmin/3600/24);
                int hmin=floor(tmin/3600-24*dmin);
                int mmin=floor(tmin/60-24*60*dmin-60*hmin);
                int smin=floor(tmin-dmin*24*3600-hmin*3600-mmin*60);
                int dmax=floor(tmax/3600/24);
                int hmax=floor(tmax/3600-24*dmax);
                int mmax=floor(tmax/60-24*60*dmax-60*hmax);
                int smax=floor(tmax-dmax*24*3600-hmax*3600-mmax*60);
                ostringstream out;
                out<<"Adv: (mean: "<<floor(100*new_mav)
                   <<"%, max: "<<floor(maxav*100)
                   <<"%, min: "<<floor(minav*100)
                   <<"%), ";
                if(tav>=0){
                    out<<"Time: (mean: "<<dav<<"-"<<hav<<":"<<mav<<":"<<sav
                       <<", max: "<<dmax<<"-"<<hmax<<":"<<mmax<<":"<<smax
                       <<", min: "<<dmin<<"-"<<hmin<<":"<<mmin<<":"<<smin
                       <<"), ";
                } else {
                    out<<"Time: (mean: N/A, max: N/A, min: N/A), ";
                }
                out<<"Rnk: (slowest: "<<slr<<", fastest: "<<far<<")";
                cout<<out.str()<<endl;
            }
        }
    }
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
        MPI_Recv(&Nq,1,MPI_INT,isready,message_save,MPI_COMM_WORLD,&status);
        for(int nq=0;nq<Nq;++nq){
            // Get name of quantity
            char strin[128];
            int len;
            MPI_Recv(&len,1,MPI_INT,isready,message_save,MPI_COMM_WORLD,&status);
            MPI_Recv(strin,len,MPI_CHAR,isready,message_save,MPI_COMM_WORLD,&status);
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
                time_t ti;
                time(&ti);
                struct tm * timeinfo;
                timeinfo=localtime(&ti);
                H5LTset_attribute_string(fout,"/","date",asctime(timeinfo));
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
#endif//EXCEPT
            }
            // get group for rank
            ostringstream gout;
            gout<<"/rank-"<<isready;
            // get current stat
            int statistics;
            MPI_Recv(&statistics,1,MPI_INT,isready,message_save,MPI_COMM_WORLD,&status);
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
            MPI_Recv(dims,2,MPI_INT,isready,message_save,MPI_COMM_WORLD,&status);
            double *buff=new double[dims[0]*dims[1]];
            MPI_Recv(buff,dims[0]*dims[1],MPI_DOUBLE,isready,message_save,
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
#endif//EXCEPT
            }
            delete [] buff;
            // set data attributes
            int Na;
            double da;
            MPI_Recv(&Na,1,MPI_INT,isready,message_save,MPI_COMM_WORLD,&status);
            for(int a=0;a<Na;++a){
                MPI_Recv(&len,1,MPI_INT,isready,message_save,MPI_COMM_WORLD,&status);
                MPI_Recv(strin,len,MPI_CHAR,isready,message_save,MPI_COMM_WORLD,&status);
                MPI_Recv(&da,1,MPI_DOUBLE,isready,message_save,MPI_COMM_WORLD,&status);
                H5LTset_attribute_double(fout,dout.str().c_str(),strin,&da,1);
            }
            e=H5Fclose(fout);
        }
    } else {
        map<string,MatStream >::iterator it=m_streams.begin();
        int Nq=m_streams.size();
        MPI_Send(&Nq,1,MPI_INT,0,message_save,MPI_COMM_WORLD);
        while(it!=m_streams.end()){
            // send name of quantity.
            int len=it->first.size()+1;
            MPI_Send(&len,1,MPI_INT,0,message_save,MPI_COMM_WORLD);
            char *qname=new char[len];
            memcpy(qname,it->first.c_str(),len*sizeof(char));
            MPI_Send(qname,len,MPI_CHAR,0,message_save,
                     MPI_COMM_WORLD);
            delete [] qname;
            // send stat
            MPI_Send(&it->second.m_stat,1,MPI_INT,0,message_save,MPI_COMM_WORLD);
            if(it->second.m_stat==StatPerSample())
                it->second.m_stat=0;
            int dims[2]={it->second.m_nrow,it->second.m_ncol};
            MPI_Send(dims,2,MPI_INT,0,message_save,MPI_COMM_WORLD);
            MPI_Send(it->second.m_mat,dims[0]*dims[1],MPI_DOUBLE,0,message_save,
                     MPI_COMM_WORLD);
            // send data attributes
            map<string,double>::iterator ait=m_dataattr.begin();
            len=m_dataattr.size();
            MPI_Send(&len,1,MPI_INT,0,message_save,MPI_COMM_WORLD);
            while(ait!=m_dataattr.end()){
                len=ait->first.size()+1;
                MPI_Send(&len,1,MPI_INT,0,message_save,MPI_COMM_WORLD);
                char *aname=new char[len];
                memcpy(aname,ait->first.c_str(),len*sizeof(char));
                MPI_Send(aname,len,MPI_CHAR,0,message_save,MPI_COMM_WORLD);
                delete [] aname;
                MPI_Send(&(ait->second),1,MPI_DOUBLE,0,message_save,MPI_COMM_WORLD);
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
            time_t ti;
            time(&ti);
            struct tm* timeinfo;
            timeinfo=localtime(&ti);
            H5LTset_attribute_string(fout,"/","date",asctime(timeinfo));
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
#endif//EXCEPT
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
#endif//EXCEPT
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

hid_t FileManager::WriteSimple(string filename)
{
    ostringstream fn;
    fn<<m_dir<<"/"<<m_num<<"-"<<filename<<".h5";
    hid_t fout=H5Fcreate(fn.str().c_str(),
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
    H5LTset_attribute_string(fout,"/","type",filename.c_str());
    return fout;
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
