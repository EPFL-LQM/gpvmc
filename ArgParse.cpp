#include "ArgParse.h"
#include <stdexcept>

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
