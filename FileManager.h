#ifndef _FILEMANAGER_H
#define _FILEMANAGER_H
#include <hdf5.h>
#include <iostream>
#include <fcntl.h>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <vector>

class MatStream;

/*! \brief class to keep track and perform output.
 *
 */

class FileManager
{
    private:
        int m_num;
        std::string m_suffix;
        std::string m_dir;
        double m_compl;
        int m_verbose;
        int m_total;
        std::map<std::string,MatStream> m_streams;
        std::map<std::string,double> m_fileattr;
        std::map<std::string,std::string> m_file_str_attr;
        std::map<std::string,double> m_dataattr;
        int m_stat_per_sample;
    public:
        FileManager(const std::string& dir="", const int& num=-1);
        enum {message_comm=0,message_monitor=1,message_save=2,message_loop=3};
        MatStream& FileStream(std::string basename);
        void FileAttribute(std::string attr, double val);
        void FileAttribute(std::string attr, std::string val);
        void DataAttribute(std::string attr, double val);
        int& StatPerSample();
        int& Verbose() {return m_verbose;}
        int Prefix() {return m_num;}
#ifdef USEMPI
        void MainLoop();
#endif
        void Monitor(int rank,double percents=0, double total_time=0);
        double& MonitorCompletion() {return m_compl;}
        int& MonitorTotal() {return m_total;}
        void Write(int isready=0);
        static void EmergencyClose(int signum);
        hid_t WriteSimple(std::string filename);
};

class MatStream
{
    public:
        int m_nrow;
        int m_ncol;
        double *m_mat;
        int m_stat;
        MatStream()
            :m_nrow(0), m_ncol(0), m_mat(0),m_stat(0) {}
        ~MatStream() {delete [] m_mat;}
        void Write(const int& m, const int& n, const double* mat)
        {
            if(m_mat) {
                if(!m_stat){
                    delete [] m_mat;
                    m_mat=new double[m*n];
                    m_nrow=m;
                    m_ncol=n;
                    std::memcpy(m_mat,mat,m*n*sizeof(double));
                    m_stat=1;
                } else {
                    if(m!=m_nrow || n!=m_ncol){
#ifdef EXCEPT
                        throw std::logic_error("MatStream::Write: sizes must match.");
#else
                        std::cerr<<"MatStream::Write: sizes must match."<<std::endl;
                        abort();
#endif
                    }
                    for(int a=0;a<m_nrow*m_ncol;++a)
                        m_mat[a]=(m_stat*m_mat[a]+mat[a])/(m_stat+1);
                    ++m_stat;
                }
            } else {
                m_mat=new double[m*n];
                m_nrow=m;
                m_ncol=n;
                std::memcpy(m_mat,mat,m*n*sizeof(double));
                m_stat=1;
            }
        }
        int& Stat()
        {
            return m_stat;
        }
};

#endif//_FILEMANAGER_H
