#include <iostream>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <sstream>

using namespace std;

bool parse(int argc, char* argv[], int& Nsamp, vector<string>& infile, string& outfile);
bool open_files(const vector<string>& infile, const string& outfile, vector<hid_t>& fin, hid_t& fout);
void close_files(const vector<hid_t>& fin, hid_t fout);
void bunch(int Nsamp, const vector<int>& statsin, vector<int>& statsout, vector<vector<int> >& args);
void getstats(const vector<hid_t>& in, vector<int>& sout, vector<hid_t>& idout, vector<string>& pathout);
void copy_attributes(hid_t iid,hid_t oid);
herr_t copy_att_cb(hid_t loc, const char* attr_name, const H5A_info_t* ainfo, void* op_data);
void concatenate(const vector<hid_t> idin, const vector<string>& pathin, const vector<int>& statsout, const vector<vector<int> > args, hid_t fout);

struct attr_op_data {
    hid_t oid;
};

int main(int argc, char *argv[]){
    string outfile;
    vector<string> infile;
    vector<hid_t> fin;
    hid_t fout(-1);
    int Nsamp=0;
    bool success=parse(argc,argv,Nsamp,infile,outfile);
    if(!success) return 1;
    success=open_files(infile,outfile,fin,fout);
    if(!success) return 2;
    vector<int> statsin,statsout;
    vector<vector<int> > args;
    vector<hid_t> idin;
    vector<string> pathin;
    getstats(fin,statsin,idin,pathin);
    if(Nsamp==0) Nsamp=statsin.size();
    bunch(Nsamp,statsin,statsout,args);
    copy_attributes(fin[0],fout);
    concatenate(idin,pathin,statsout,args,fout);
    close_files(fin,fout);
    return 0;
}

bool parse(int argc, char* argv[], int& Nsamp, vector<string>& infile, string& outfile)
{
    if(argc<3){
        cout<<"Usage:"<<endl
            <<"concatenate [-N <Nsamp>] infile1 infile2 ... outfile"<<endl;
        return false;
    }
    infile.clear();
    Nsamp=0;
    int i=1;
    while(i<argc){
        if(string(argv[i])=="-N"){
            if(argc<5){
                cerr<<"Not enough input arguments"<<endl;
                return false;
            }
            ++i;
            try {
                Nsamp=stoi(string(argv[i]));
            } catch(const invalid_argument& err) {
                cerr<<"Invalid argument on \"-N\": "<<err.what()<<endl;
                return false;
            }
        } else if(i==(argc-1)){
            outfile=string(argv[i]);
        } else {
            infile.push_back(argv[i]);
        }
        ++i;
    }
    return true;
}

bool open_files(const vector<string>& infile, const string& outfile, vector<hid_t>& fin, hid_t& fout)
{
    fin.resize(infile.size());
    for(size_t i=0;i<infile.size();++i){
        fin[i]=H5Fopen(infile[i].c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
        if(fin[i]<0){
            cerr<<"Could not open file \""<<infile[i]<<"\""<<endl;
            return false;
        }
    }
    fout=H5Fcreate(outfile.c_str(),H5F_ACC_EXCL,H5P_DEFAULT,H5P_DEFAULT);
    if(fout<0){
        cerr<<"Could not create file \""<<outfile<<"\""<<endl;
        return false;
    }
    return true;
}

void close_files(const vector<hid_t>& fin,hid_t fout)
{
    for(size_t i=0;i<fin.size();++i){
        H5Fclose(fin[i]);
    }
    H5Fclose(fout);
}

void bunch(int Nsamp, const vector<int>& statsin, vector<int>& statsout, vector<vector<int> >& args)
{
    vector<int> sin(statsin);
    vector<int> idx(statsin.size());
    for(size_t i=0;i<sin.size();++i) idx[i]=i;
    sort(idx.begin(),idx.end(),[&](int i,int j){return sin[i]<sin[j];});
    statsout.resize(Nsamp);
    args.resize(Nsamp);
    int n=0;
    while(sin.size()>size_t(2*Nsamp+1)){
        args[n].push_back(idx.front());
        args[n].push_back(idx.back());
        statsout[n]+=sin[idx.front()]+sin[idx.back()];
        idx.assign(idx.begin()+1,idx.end()-1);
        sin.assign(sin.begin()+1,sin.end()-1);
        n+=1;
        if(n==Nsamp)
            n=0;
    }
    while(sin.size()){
        int bmin=min_element(statsout.begin(),statsout.end())-statsout.begin();
        args[bmin].push_back(idx.back());
        statsout[bmin]+=sin.back();
        sin.pop_back();
        idx.pop_back();
    }
}

void getstats(const vector<hid_t>& in, vector<int>& sout, vector<hid_t>& idout, vector<string>& pathout)
{
    sout.clear();
    idout.clear();
    pathout.clear();
    for(size_t f=0;f<in.size();++f){
        int g=1;
        string gn="/rank-1";
        while(H5Lexists(in[f],gn.c_str(),H5P_DEFAULT)){
            int d=0;
            string dn=gn+"/data-0";
            while(H5Lexists(in[f],dn.c_str(),H5P_DEFAULT)){
                idout.push_back(in[f]);
                pathout.push_back(dn);
                //get stat attribute
                double adat;
                H5LTget_attribute_double(in[f],dn.c_str(),"statistics",&adat);
                sout.push_back(adat);
                ++d;
                dn=gn+"/data-"+to_string(d);
            }
            ++g;
            gn="/rank-"+to_string(g);
        }
    }
}

void copy_attributes(hid_t iid,hid_t oid)
{
    struct attr_op_data op;
    op.oid=oid;
    H5Aiterate2(iid,H5_INDEX_NAME,H5_ITER_NATIVE,NULL,copy_att_cb,(void*)&op);
}

herr_t copy_att_cb(hid_t loc, const char* attr_name, const H5A_info_t* ainfo, void* op_data)
{
    hid_t attr=H5Aopen(loc,attr_name,H5P_DEFAULT);
    hid_t type=H5Aget_type(attr);
    hid_t space=H5Aget_space(attr);
    if(attr<0 || type<0 || space<0){
        cerr<<"Error opening attribute \""<<attr_name<<"\"."<<endl;
        return -1;
    }
    hid_t oattr=H5Acreate2(((struct attr_op_data*)op_data)->oid,attr_name,type,space,H5P_DEFAULT,H5P_DEFAULT);
    if(oattr<0){
        cerr<<"Error creating attribute \""<<attr_name<<"\"."<<endl;
        return -1;
    }
    hsize_t size=H5Aget_storage_size(attr);
    char* buf=new char[size];
    herr_t status;
    status=H5Aread(attr,type,buf);
    if(status<0){
        cerr<<"Error reading attribute \""<<attr_name<<"\"."<<endl;
        return -1;
    }
    status=H5Awrite(oattr,type,buf);
    if(status<0){
        cerr<<"Error writing attribute \""<<attr_name<<"\"."<<endl;
        return -1;
    }
    delete [] buf;
    H5Aclose(attr);
    H5Aclose(oattr);
    H5Tclose(type);
    H5Sclose(space);
    return 0;
}

void concatenate(const vector<hid_t> idin, const vector<string>& pathin, const vector<int>& statsout, const vector<vector<int > > args, hid_t fout)
{
    hid_t g=H5Gcreate2(fout,"/rank-1",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    for(size_t s=0;s<args.size();++s){
        string pathout=string("/rank-1/data-")+to_string(s);
        hid_t din,dout,space,type;
        hsize_t dims[2];
        din=H5Dopen2(idin[args[s][0]],pathin[args[s][0]].c_str(),H5P_DEFAULT);
        space=H5Dget_space(din);
        H5Sget_simple_extent_dims(space,dims,NULL);
        vector<double> buf(dims[0]*dims[1]),conc(dims[0]*dims[1]);
        type=H5Dget_type(din);
        H5Dread(din,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf.data());
        conc=buf;
        H5Dclose(din);
        for(size_t si=1;si<args[s].size();++si){
            din=H5Dopen2(idin[args[s][si]],pathin[args[s][si]].c_str(),H5P_DEFAULT);
            H5Dread(din,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,buf.data());
            for(size_t d=0;d<buf.size();++d) conc[d]+=buf[d];
            H5Dclose(din);
        }
        for(size_t d=0;d<conc.size();++d) conc[d]/=args[s].size();
        dout=H5Dcreate2(fout,pathout.c_str(),type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dout,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,conc.data());
        H5Dclose(dout);
        H5LTset_attribute_int(fout,pathout.c_str(),"statistics",&statsout[s],1);
        H5Sclose(space);
        H5Tclose(type);
    }
    H5Gclose(g);
}
