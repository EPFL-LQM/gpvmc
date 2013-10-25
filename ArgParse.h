#ifndef _ARGPARSE_H
#define _ARGPARSE_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>

class ArgParse {
    private:
        std::map<std::string,std::string> m_streams;
    public:
        ArgParse(int argc, char* argv[]);
        double d(const std::string& str) const;
        int i(const std::string& str) const;
        bool b(const std::string& str) const;
        std::string s(const std::string& str) const;
        bool HasArg(const std::string& str) const;
        void SetupParams(std::map<std::string,bool>& map_bool,
                         std::map<std::string,size_t>& map_size_t,
                         std::map<std::string,int>& map_int,
                         std::map<std::string,double>& map_double,
                         std::map<std::string,std::string>& map_string);
};

#endif//_ARGPARSE_H
