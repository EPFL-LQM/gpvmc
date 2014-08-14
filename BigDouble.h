#ifndef _BIGDOUBLE_H
#define _BIGDOUBLE_H
#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;

//! A class to store complex numbers and their scale separatly to avoid under- or over-flow.
class BigDouble
{
    public:
        double m_c;
        int    m_e;

        void bound()
        {
#ifdef CRAY
            if(isnan(m_c)){
#ifdef EXCEPT
                throw(std::logic_error("BigDouble::bound: NaN encountered"));
#else
                cerr<<"BigComplex::bound: NaN or Inf encountered"<<endl;
                abort();
#endif//EXCEPT
            }
            if(isinf(m_c)/* || (abs(m_c)<1e4 && abs(m_c)>1e-4)*/) return;
#else
            if(std::isnan(m_c)){
#ifdef EXCEPT
                throw(std::logic_error("BigDouble::bound: NaN encountered"));
#else
                cerr<<"BigDouble::bound: NaN encountered"<<endl;
                abort();
#endif//EXCPET
            }
            if(std::isinf(m_c)/* || (abs(m_c)<1e4 && abs(m_c)>1e-4)*/) return;
#endif//CRAY
            while(std::abs(m_c)>10.0){
                m_c/=10.0;
                m_e++;
            }
            if(std::abs(m_c)>0){
                while(std::abs(m_c)<1.0){
                    m_c*=10.0;
                    m_e--;
                }
            } else m_e=0;
        }

        BigDouble()
            :m_c(0), m_e(0)
        {}

        BigDouble(const double& c)
            :m_c(c), m_e(0)
        {
            bound();
        }
        BigDouble(const double& c, const int& e)
            :m_c(c), m_e(e)
        {
            bound();
        }
        BigDouble(const BigDouble& bc)
            :m_c(bc.m_c),m_e(bc.m_e)
        {}

        BigDouble& operator=(const BigDouble& bd)
        {
            if(this==&bd) return *this;
            m_c=bd.m_c;
            m_e=bd.m_e;
            return *this;
        }

        BigDouble& operator=(const double& d)
        {
            m_c=d;
            m_e=0;
            bound();
            return *this;
        }

        BigDouble& operator=(const int& d)
        {
            m_c=d;
            m_e=0;
            bound();
            return *this;
        }

        int exp() const {return m_e;}

        const BigDouble& operator+() const {return *this;}
        const BigDouble operator-() const {return BigDouble(-m_c,m_e);}

        operator double() const {return m_c*pow(10.0,m_e);}

        BigDouble& operator+=(const BigDouble& bd)
        {
            int mexp=std::max(m_e, bd.m_e);
            m_c=m_c*pow(10.0,m_e-mexp)+bd.m_c*pow(10.0,bd.m_e-mexp);
            m_e=mexp;
            bound();
            return *this;
        }

        BigDouble& operator-=(const BigDouble& bd)
        {
            int mexp=std::max(m_e, bd.m_e);
            m_c=m_c*pow(10.0,m_e-mexp)-bd.m_c*pow(10.0,bd.m_e-mexp);
            m_e=mexp;
            bound();
            return *this;
        }

        BigDouble& operator*=(const BigDouble& bd)
        {
            m_c*=bd.m_c;
            m_e+=bd.m_e;
            bound();
            return *this;
        }

        BigDouble& operator/=(const BigDouble& bd)
        {
            m_c/=bd.m_c;
            m_e-=bd.m_e;
            bound();
            return *this;
        }

        friend const BigDouble operator+(const BigDouble& left, const BigDouble& right);
        friend const BigDouble operator-(const BigDouble& left, const BigDouble& right);
        friend const BigDouble operator*(const BigDouble& left, const BigDouble& right);
        friend const BigDouble operator/(const BigDouble& left, const BigDouble& right);

        friend const BigDouble operator+(const BigDouble& left, const double& right);
        friend const BigDouble operator-(const BigDouble& left, const double& right);
        friend const BigDouble operator*(const BigDouble& left, const double& right);
        friend const BigDouble operator/(const BigDouble& left, const double& right);

        friend const BigDouble operator+(const double& left, const BigDouble& right);
        friend const BigDouble operator-(const double& left, const BigDouble& right);
        friend const BigDouble operator*(const double& left, const BigDouble& right);
        friend const BigDouble operator/(const double& left, const BigDouble& right);

        friend bool operator==(const BigDouble& left, const BigDouble& right);
        friend bool operator!=(const BigDouble& left, const BigDouble& right);
        friend bool operator<(const BigDouble& left, const BigDouble& right);
        friend bool operator>(const BigDouble& left, const BigDouble& right);

        friend bool operator==(const double& left, const BigDouble& right);
        friend bool operator!=(const double& left, const BigDouble& right);
        friend bool operator<(const double& left, const BigDouble& right);
        friend bool operator>(const double& left, const BigDouble& right);

        friend bool operator==(const BigDouble& left, const double& right);
        friend bool operator!=(const BigDouble& left, const double& right);
        friend bool operator<(const BigDouble& left, const double& right);
        friend bool operator>(const BigDouble& left, const double& right);

        friend std::ostream& operator<<(std::ostream& os,const BigDouble& bc);
};

inline const BigDouble abs(const BigDouble& bd)
{
    return BigDouble(std::abs(bd.m_c),bd.m_e);
}

inline const BigDouble operator+(const BigDouble& left, const BigDouble& right)
{
    int mexp=std::max(left.m_e,right.m_e);
    return BigDouble(left.m_c*pow(10,left.m_e-mexp)+right.m_c*pow(10,right.m_e-mexp),mexp);
}

inline const BigDouble operator-(const BigDouble& left, const BigDouble& right)
{
    int mexp=std::max(left.m_e,right.m_e);
    return BigDouble(left.m_c*pow(10,left.m_e-mexp)-right.m_c*pow(10,right.m_e-mexp),mexp);
}

inline const BigDouble operator*(const BigDouble& left, const BigDouble& right)
{
    return BigDouble(left.m_c*right.m_c,left.m_e+right.m_e);
}

inline const BigDouble operator/(const BigDouble& left, const BigDouble& right)
{
    return BigDouble(left.m_c/right.m_c,left.m_e-right.m_e);
}

inline const BigDouble operator+(const BigDouble& left, const double& right)
{
    return BigDouble(left.m_c+right*pow(10,-left.m_e),left.m_e);
}

inline const BigDouble operator-(const BigDouble& left, const double& right)
{
    return BigDouble(left.m_c-right*pow(10,-left.m_e),-left.m_e);
}

inline const BigDouble operator*(const BigDouble& left, const double& right)
{
    return BigDouble(left.m_c*right,left.m_e);
}

inline const BigDouble operator/(const BigDouble& left, const double& right)
{
    return BigDouble(left.m_c/right,left.m_e);
}

inline const BigDouble operator+(const double& left, const BigDouble& right)
{
    return BigDouble(right.m_c+left*pow(10,-right.m_e),right.m_e);
}

inline const BigDouble operator-(const double& left, const BigDouble& right)
{
    return BigDouble(right.m_c-left*pow(10,-right.m_e),-right.m_e);
}

inline const BigDouble operator*(const double& left, const BigDouble& right)
{
    return BigDouble(right.m_c*left,right.m_e);
}

inline const BigDouble operator/(const double& left, const BigDouble& right)
{
    return BigDouble(right.m_c/left,right.m_e);
}

inline bool operator==(const BigDouble& left, const BigDouble& right)
{
    BigDouble l(left), r(right);
    l.bound(); r.bound();
    return l.m_c==r.m_c && l.m_e==r.m_e;
}

inline bool operator!=(const BigDouble& left, const BigDouble& right)
{
    BigDouble l(left), r(right);
    l.bound(); r.bound();
    return l.m_c!=r.m_c || l.m_e!=r.m_e;
}

inline bool operator<(const BigDouble& left, const BigDouble& right)
{
    BigDouble l(left), r(right);
    l.bound(); r.bound();
    return l.m_e<r.m_e || (l.m_e==r.m_e && l.m_c<r.m_c);
}

inline bool operator<=(const BigDouble& left, const BigDouble& right)
{
    BigDouble l(left), r(right);
    l.bound(); r.bound();
    return l.m_e<r.m_e || (l.m_e==r.m_e && l.m_c<=r.m_c);
}

inline bool operator>(const BigDouble& left, const BigDouble& right)
{
    BigDouble l(left), r(right);
    l.bound(); r.bound();
    return l.m_e>r.m_e || (l.m_e==r.m_e && l.m_c>r.m_c);
}

inline bool operator>=(const BigDouble& left, const BigDouble& right)
{
    BigDouble l(left), r(right);
    l.bound(); r.bound();
    return l.m_e>r.m_e || (l.m_e==r.m_e && l.m_c>=r.m_c);
}

inline bool operator==(const double& left, const BigDouble& right)
{
    return left==double(right);
}

inline bool operator!=(const double& left, const BigDouble& right)
{
    return left!=double(right);
}

inline bool operator<(const double& left, const BigDouble& right)
{
    return left<double(right);
}

inline bool operator>(const double& left, const BigDouble& right)
{
    return left>double(right);
}

inline bool operator<=(const double& left, const BigDouble& right)
{
    return left<=double(right);
}

inline bool operator>=(const double& left, const BigDouble& right)
{
    return left>=double(right);
}

inline bool operator==(const BigDouble& left, const double& right)
{
    return double(left)==right;
}

inline bool operator!=(const BigDouble& left, const double& right)
{
    return double(left)!=right;
}

inline bool operator<(const BigDouble& left, const double& right)
{
    return double(left)<right;
}

inline bool operator>(const BigDouble& left, const double& right)
{
    return double(left)>right;
}

inline bool operator<=(const BigDouble& left, const double& right)
{
    return double(left)<=right;
}

inline bool operator>=(const BigDouble& left, const double& right)
{
    return double(left)>=right;
}

inline std::ostream& operator<<(std::ostream& os,const BigDouble& bc)
{
    os<<bc.m_c<<"e"<<bc.m_e;
    return os;
}

#endif//_BIGDOUBLE_H
