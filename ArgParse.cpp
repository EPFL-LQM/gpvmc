#include "ArgParse.h"

using namespace std;

ArgParse::ArgParse(int argc, char* argv[])
{
    for(int a=1;a<argc;++a){
        if(argv[a][0]=='-'){
            string in(argv[a]),trig,val;
            trig=in.substr(1,in.find('='));
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
                if(in[0]=='#') continue;
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
    return stod(m_streams.at(str));
}

int ArgParse::i(const string& str) const
{
    return stoi(m_streams.at(str));
}

bool ArgParse::b(const string& str) const
{
    return stoi(m_streams.at(str));
}
 string ArgParse::s(const string& str) const
{
    return m_streams.at(str);
}
