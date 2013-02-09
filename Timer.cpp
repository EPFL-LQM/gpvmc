#include "Timer.h"
#include <sstream>

using namespace std;

Timer Timer::s_ti;

Timer::tim::~tim()
{
    map<string,tim*>::iterator it=m_timtree.begin();
    while(it!=m_timtree.end()){
        delete it->second;
        it++;
    }
}

void Timer::tic(string timer)
{
    size_t p=timer.find("/");
    tim* t=&s_ti.m_timers[timer.substr(0,p)];
    if(p!=string::npos)
        t=t->get_timer(timer.substr(p+1));
#ifdef USEPARA
    t->m_t0=omp_get_wtime();
#else
    t->m_t0=clock();
#endif
}

double Timer::toc(string timer)
{
    size_t p=timer.find("/");
    tim* t=&s_ti.m_timers[timer.substr(0,p)];
    if(p!=string::npos)
        t=t->get_timer(timer.substr(p+1));
#ifdef USEPARA
    t->m_t+=omp_get_wtime()-t->m_t0;
#else
    t->m_t+=(clock()-t->m_t0)/double(CLOCKS_PER_SEC);
#endif
    return t->m_t;
}

const Timer& Timer::timer()
{
    return s_ti;
}

double Timer::timer(std::string ti)
{
    size_t p=ti.find("/");
    if(p==string::npos)
        return s_ti.m_timers[ti].m_t;
    else
        return s_ti.m_timers[ti.substr(0,p)].get_timer(ti.substr(p+1))->m_t;
}

void Timer::reset(std::string ti)
{
    s_ti.m_timers[ti].m_t=0;
}

Timer::tim* Timer::tim::get_timer(string timer)
{
    size_t p=timer.find("/");
    if(p==string::npos){
        if(m_timtree.find(timer)==m_timtree.end())
            m_timtree[timer]=new tim;
        return m_timtree[timer];
    } else {
        if(m_timtree.find(timer.substr(0,p))==m_timtree.end())
            m_timtree[timer.substr(0,p)]=new tim;
        return m_timtree[timer.substr(0,p)]->get_timer(timer.substr(p+1));
    }
}

string Timer::report()
{
    ostringstream strout;
    map<string,tim>::iterator it=s_ti.m_timers.begin();
    while(it!=s_ti.m_timers.end())
    {
        strout<<it->first<<": "<<it->second.report()<<endl;
        ++it;
    }
    return strout.str();
}

string Timer::tim::report()
{
    ostringstream strout;
    strout<<"total: "<<m_t<<", self: "<<get_self_time()<<endl;
    map<string,tim*>::iterator it=m_timtree.begin();
    while(it!=m_timtree.end()){
        strout<<it->first<<": "<<it->second->report();
        it++;
    }
    return strout.str();
}

double Timer::tim::get_self_time()
{
    double out=m_t;
    map<string,tim*>::iterator it=m_timtree.begin();
    while(it!=m_timtree.end()){
        out-=it->second->m_t;
        it++;
    }
    return out;
}
