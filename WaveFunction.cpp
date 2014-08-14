#include "WaveFunction.h"
#include "FileManager.h"
#include <hdf5_hl.h>
#include <hdf5.h>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

WaveFunction::WaveFunction(FileManager* fm)
    :State(fm)
{}

void WaveFunction::build_base(const uint_vec_t& Npt, const uint_vec_t& Nfs)
{
    State::build_base(Npt,Nfs);
    m_cache=vector<vector<complex<double> > >(m_Nfl);
    for(size_t fl=0;fl<m_Nfl;++fl)
        m_cache[fl]=vector<complex<double> >(size_t(pow(m_Nfs[fl],2)));
}

WaveFunction::~WaveFunction()
{}

void WaveFunction::init_matrices()
{
    for(size_t fl=0;fl<m_Nfl;++fl){
        for(size_t f=0;f<m_Nfs[fl];++f){
            for(size_t R=0;R<m_Nfs[fl];++R){
                m_cache[fl][f*m_Nfs[fl]+R]=
                    this->matrix_element(f,R,fl);
            }
        }
    }
}

void WaveFunction::AddState(const vector<uint_vec_t>& f)
{
    if(!m_exc.size()){
        for(size_t fl=0;fl<m_Nfl;++fl){
            size_t sti=0;
            for(size_t s=0;s<m_Nfs[fl];++s){
                if(f[fl][s]!=m_Npt[fl]){
                    m_fock[fl][s]=sti;
                    m_part[fl][sti]=s;
                    ++sti;
                }
            }
        }
        m_c_exc=0;
    }
    hop_path_t em;// empty hop
    // span hop path matrix with empty hops
    m_exc.push_back(vector<vector<hop_path_t> >(m_exc.size()+1,vector<hop_path_t>(m_Nfl,em)));
    for(size_t i=0;i<m_exc.size()-1;++i)
        m_exc[i].push_back(vector<hop_path_t>(m_Nfl,em));
    for(size_t flav=0;flav<m_Nfl;++flav){
        uint_vec_t state(m_Nfs[flav],0);
        for(size_t i=0;i<m_exc.size()-1;++i){
            get_hop_state(m_exc[0][i][flav],m_fock[flav],m_Npt[flav],state);
            get_hop_path(f[flav],state,m_Npt[flav],m_exc.back()[i][flav]);
            get_hop_path(state,f[flav],m_Npt[flav],m_exc[i].back()[flav]);
        }
    }
}

void WaveFunction::Hop(size_t khop)
{
    State::Hop(m_exc[m_c_exc][khop]);
    m_c_exc=khop;
}

const vector<vector<hop_path_t> >& WaveFunction::GetHops() const
{
    return m_exc[m_c_exc];
}

const vector<hop_path_t>& WaveFunction::GetHop(size_t k) const
{
    return m_exc[m_c_exc][k];
}

void WaveFunction::Save()
{
    int rank=0, size=1;
#ifdef USEMPI
    MPI_Comm comm=m_filemanager->GetCalcComm();
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);
#endif//USEMPI
    if(rank==0){
        hid_t wfile=m_filemanager->WriteSimple("WaveFunction");
        for(size_t flav=0;flav<m_Nfl;++flav){
            vector<int> fock(m_exc.size()*m_Nfs[flav]);
            hsize_t dims[2]={m_exc.size(),m_Nfs[flav]};
            for(size_t e=0;e<m_exc.size();++e){
                uint_vec_t focke(m_Nfs[flav]);
                get_hop_state(GetHop(e)[flav],m_fock[flav],m_Npt[flav],focke);
                for(size_t s=0;s<m_Nfs[flav];++s)
                    fock[e*m_Nfs[flav]+s]=(focke[s]!=m_Npt[flav]);
            }
            ostringstream sn;
            sn<<"states_"<<flav;
            H5LTmake_dataset_int(wfile,sn.str().c_str(),2,dims,fock.data());
            // open dataset to add attributes
            hid_t dataset=H5Dopen(wfile,sn.str().c_str(),H5P_DEFAULT);
            // create attributes to store quantum numbers.
            hsize_t dim=m_Nfs[flav];
            hid_t ds;
            ds=H5Screate_simple(1,&dim,NULL);
            map<string,size_t> qn;
            map<string,hid_t> ats;
            map<string,vector<unsigned int> > qns;
            map<string,size_t>::iterator qniter;
            this->quantum_numbers(0,flav,qn);
            for(qniter=qn.begin();qniter!=qn.end();qniter++){
                ats[qniter->first]=H5Acreate(dataset,qniter->first.c_str(),H5T_STD_U32BE,ds,H5P_DEFAULT,H5P_DEFAULT);
                qns[qniter->first]=vector<unsigned int>(m_Nfs[flav]);
            }
            H5Sclose(ds);
            for(size_t f=0;f<m_Nfs[flav];++f){
                this->quantum_numbers(f,flav,qn);
                for(qniter=qn.begin();qniter!=qn.end();qniter++){
                    qns[qniter->first][f]=qniter->second;
                }
            }
            for(map<string,hid_t>::iterator atiter=ats.begin();atiter!=ats.end();atiter++){
                H5Awrite(atiter->second,H5T_NATIVE_UINT,&qns[atiter->first][0]);
                H5Aclose(atiter->second);
            }
            H5Dclose(dataset);
        }
        H5Fclose(wfile);
    }
}
