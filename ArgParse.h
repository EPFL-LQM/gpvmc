#ifndef _ARGPARSE_H
#define _ARGPARSE_H

#include <iostream>
#include <string>
#include <map>
#include <vector>

int ArgParse(int argc, char* argv[],
             std::map<std::string,double>& double_opts,
             std::map<std::string,int>& int_opts,
             std::map<std::string,size_t>& uint_opts,
             std::map<std::string,bool>& bool_opts,
             std::map<std::string,std::string>& string_opts);

#endif//_ARGPARSE_H
