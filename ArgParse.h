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
};

#endif//_ARGPARSE_H
