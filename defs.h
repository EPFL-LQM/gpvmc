#include <vector>
#include <utility>
#ifndef _DEFS_H
#define _DEFS_H
typedef std::pair<std::size_t,std::size_t> hop_t;
typedef std::vector<hop_t> hop_path_t;
typedef std::vector<std::size_t> uint_vec_t;
bool uint_vec_t_comp(const uint_vec_t& a,const uint_vec_t& b);
#endif
