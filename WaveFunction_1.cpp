#include "WaveFunction_1.h"
#include "FileManager.h"
#include <hdf5_hl.h>
#include <hdf5.h>
#include <cstring>
#ifdef USEMPI
#include <mpi.h>
#endif

using namespace std;

WaveFunction_1::WaveFunction_1()
{}

void WaveFunction_1::build_base(const uint_vec_t& Npt, const uint_vec_t& Nfs)
{
    State_1::build_base(Npt,Nfs);
    m_cache=vector<vector<complex<double> > >(Npt.size());
    for(size_t fl=0;fl<m_Nfl;++fl)
        m_cache[fl]=vector<complex<double> >(size_t(pow(m_Nfs[fl],2)));
}

WaveFunction_1::~WaveFunction_1()
{}

void WaveFunction_1::init_matrices()
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

void WaveFunction_1::AddState(const vector<uint_vec_t>& f)
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

const vector<vector<hop_path_t> >& WaveFunction_1::GetHops() const
{
    return m_exc[m_c_exc];
}

const vector<hop_path_t>& WaveFunction_1::GetHop(size_t k) const
{
    return m_exc[m_c_exc][k];
}

void WaveFunction_1::Save(FileManager* fm)
{
    int rank=0, size=1;
#ifdef USEMPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#endif//USEMPI
    if(size==1 || rank==1){ //not root process except if non-mpi
        hid_t wfile=fm->WriteSimple("WaveFunction");
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
        }
        H5Fclose(wfile);
    }
}
