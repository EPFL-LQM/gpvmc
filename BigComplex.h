#ifndef _BIGCOMPLEX_H
#define _BIGCOMPLEX_H
#include <complex>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "BigDouble.h"

//! A class to store complex numbers and their scale separatly to avoid under- or over-flow.
class BigComplex
{
    public:
        std::complex<double> m_c;
        int             m_e;

        void bound()
        {
            bool unbounded=false;
#ifdef CRAY
            unbounded=(isnan(real(m_c)) || isinf(real(m_c)) || isnan(imag(m_c)) || isinf(imag(m_c)));
#else
            unbounded=(std::isnan(real(m_c)) || std::isinf(real(m_c)) || std::isnan(imag(m_c)) || std::isinf(imag(m_c)));
#endif
            if(unbounded){
#ifdef EXCEPT
                throw(std::logic_error("BigComplex::bound: NaN or Inf encountered"));
#else
                std::cerr<<"BigComplex::bound: NaN or Inf encountered"<<std::endl;
                abort();
#endif
            }
            //if(norm(m_c) < 1e8 && norm(m_c)>1e-8) return;
            while(std::abs(m_c)>10.0){
                m_c/=10.0;
                m_e++;
            }
            if(std::abs(m_c)>0){
                while(std::abs(m_c)<1.0){
                    m_c*=10.0;
                    m_e--;
                }
            }

        }

        bool isnaninf()
        {
#ifdef CRAY
            return isnan(real(m_c)) || isinf(real(m_c)) || isnan(imag(m_c)) || isinf(imag(m_c));
#else
            return std::isnan(real(m_c)) || std::isinf(real(m_c)) || std::isnan(imag(m_c)) || std::isinf(imag(m_c));
#endif
        }

        BigComplex(const std::complex<double>& c, int e=0)
            :m_c(c), m_e(e)
        {
            bound();
        }

        BigComplex(const double& r=0, const double& i=0, const int& e=0)
            :m_c(std::complex<double>(r,i)),m_e(e)
        {
            bound();
        }

        BigComplex(const BigComplex& bc):m_c(bc.m_c),m_e(bc.m_e){}

        BigComplex& operator=(const BigComplex& bc) 
        {
            if(this==&bc) return *this;
            m_c=bc.m_c;
            m_e=bc.m_e;
            return *this;
        }

        BigComplex& operator=(const std::complex<double>& c)
        {
            m_c=c;
            m_e=0;
            bound();
            return *this;
        }

        BigComplex& operator=(const BigDouble& bd)
        {
            m_c=bd.m_c;
            m_e=bd.m_e;
            return *this;
        }

        BigComplex& operator=(const int& d)
        {
            m_c=d;
            m_e=0;
            bound();
            return *this;
        }

        BigComplex& operator=(const double& d)
        {
            m_c=d;
            m_e=0;
            bound();
            return *this;
        }

        int exp() const {return m_e;}

        operator std::complex<double>() const {return m_c*pow(10.0,m_e);}

        BigComplex& operator+=(const BigComplex& bd)
        {
            int mexp=std::max(m_e, bd.m_e);
            m_c=m_c*pow(10.0,m_e-mexp)+bd.m_c*pow(10.0,bd.m_e-mexp);
            m_e=mexp;
            bound();
            return *this;
        }

        BigComplex& operator-=(const BigComplex& bd)
        {
            int mexp=std::max(m_e, bd.m_e);
            m_c=m_c*pow(10.0,m_e-mexp)-bd.m_c*pow(10.0,bd.m_e-mexp);
            m_e=mexp;
            bound();
            return *this;
        }

        BigComplex& operator*=(const BigComplex& bd)
        {
            m_c*=bd.m_c;
            m_e+=bd.m_e;
            bound();
            return *this;
        }

        BigComplex& operator/=(const BigComplex& bd)
        {
            m_c/=bd.m_c;
            m_e-=bd.m_e;
            bound();
            return *this;
        }

        BigComplex& operator+=(const BigDouble& bd)
        {
            int mexp=std::max(m_e, bd.m_e);
            m_c=m_c*pow(10.0,m_e-mexp)+bd.m_c*pow(10.0,bd.m_e-mexp);
            m_e=mexp;
            bound();
            return *this;
        }

        BigComplex& operator-=(const BigDouble& bd)
        {
            int mexp=std::max(m_e, bd.m_e);
            m_c=m_c*pow(10.0,m_e-mexp)-bd.m_c*pow(10.0,bd.m_e-mexp);
            m_e=mexp;
            bound();
            return *this;
        }

        BigComplex& operator*=(const BigDouble& bd)
        {
            m_c*=bd.m_c;
            m_e+=bd.m_e;
            bound();
            return *this;
        }

        BigComplex& operator/=(const BigDouble& bd)
        {
            m_c/=bd.m_c;
            m_e-=bd.m_e;
            bound();
            return *this;
        }

        BigComplex& operator+=(const double& d)
        {
            m_c+=d;
            bound();
            return *this;
        }

        BigComplex& operator-=(const double& d)
        {
            m_c-=d;
            bound();
            return *this;
        }

        BigComplex& operator*=(const double& d)
        {
            m_c*=d;
            bound();
            return *this;
        }

        BigComplex& operator/=(const double& d)
        {
            m_c/=d;
            bound();
            return *this;
        }

        friend const BigComplex operator+(const BigComplex& left, const BigComplex& right);
        friend const BigComplex operator-(const BigComplex& left, const BigComplex& right);
        friend const BigComplex operator*(const BigComplex& left, const BigComplex& right);
        friend const BigComplex operator/(const BigComplex& left, const BigComplex& right);

        friend const BigComplex operator+(const BigComplex& left, const BigDouble& right);
        friend const BigComplex operator-(const BigComplex& left, const BigDouble& right);
        friend const BigComplex operator*(const BigComplex& left, const BigDouble& right);
        friend const BigComplex operator/(const BigComplex& left, const BigDouble& right);

        friend const BigComplex operator+(const BigComplex& left, const double& right);
        friend const BigComplex operator-(const BigComplex& left, const double& right);
        friend const BigComplex operator*(const BigComplex& left, const double& right);
        friend const BigComplex operator/(const BigComplex& left, const double& right);

        friend const BigComplex operator+(const double& left, const BigComplex& right);
        friend const BigComplex operator-(const double& left, const BigComplex& right);
        friend const BigComplex operator*(const double& left, const BigComplex& right);
        friend const BigComplex operator/(const double& left, const BigComplex& right);

        friend bool operator==(const BigComplex& left, const BigComplex& right);
        friend bool operator==(const std::complex<double>& left, const BigComplex& right);
        friend bool operator==(const BigComplex& left, const std::complex<double>& right);
        friend bool operator==(const BigComplex& left, const double& right);
        friend bool operator==(const double& left, const BigComplex& right);

        friend std::ostream& operator<<(std::ostream& os,const BigComplex& bc);
};

inline const BigDouble norm(const BigComplex& bc)
{
    return BigDouble(norm(bc.m_c),2*bc.m_e);
}

inline const BigDouble abs(const BigComplex& bc)
{
    return BigDouble(std::abs(bc.m_c),bc.m_e);
}

inline const BigComplex conj(const BigComplex& bc)
{
    return BigComplex(conj(bc.m_c),bc.m_e);
}

inline const BigComplex operator+(const BigComplex& left, const BigComplex& right)
{
    int mexp=std::max(left.m_e,right.m_e);
    return BigComplex(left.m_c*pow(10,left.m_e-mexp)+right.m_c*pow(10,right.m_e-mexp),mexp);
}

inline const BigComplex operator-(const BigComplex& left, const BigComplex& right)
{
    int mexp=std::max(left.m_e,right.m_e);
    return BigComplex(left.m_c*pow(10,left.m_e-mexp)-right.m_c*pow(10,right.m_e-mexp),mexp);
}

inline const BigComplex operator*(const BigComplex& left, const BigComplex& right)
{
    return BigComplex(left.m_c*right.m_c,left.m_e+right.m_e);
}

inline const BigComplex operator/(const BigComplex& left, const BigComplex& right)
{
    return BigComplex(left.m_c/right.m_c,left.m_e-right.m_e);
}

inline const BigComplex operator+(const BigComplex& left, const BigDouble& right)
{
    int mexp=std::max(left.m_e,right.m_e);
    return BigComplex(left.m_c*pow(10,left.m_e-mexp)+right.m_c*pow(10,right.m_e-mexp),mexp);
}

inline const BigComplex operator-(const BigComplex& left, const BigDouble& right)
{
    int mexp=std::max(left.m_e,right.m_e);
    return BigComplex(left.m_c*pow(10,left.m_e-mexp)-right.m_c*pow(10,right.m_e-mexp),mexp);
}

inline const BigComplex operator*(const BigComplex& left, const BigDouble& right)
{
    return BigComplex(left.m_c*right.m_c,left.m_e+right.m_e);
}

inline const BigComplex operator/(const BigComplex& left, const BigDouble& right)
{
    return BigComplex(left.m_c/right.m_c,left.m_e-right.m_e);
}

inline const BigComplex operator+(const BigComplex& left, const double& right)
{
    return BigComplex(left.m_c+right*pow(10,-left.m_e),left.m_e);
}

inline const BigComplex operator-(const BigComplex& left, const double& right)
{
    return BigComplex(left.m_c-right*pow(10,-left.m_e),-left.m_e);
}

inline const BigComplex operator*(const BigComplex& left, const double& right)
{
    return BigComplex(left.m_c*right,left.m_e);
}

inline const BigComplex operator/(const BigComplex& left, const double& right)
{
    return BigComplex(left.m_c/right,left.m_e);
}

inline const BigComplex operator+(const double& left, const BigComplex& right)
{
    return BigComplex(right.m_c+left*pow(10,-right.m_e),right.m_e);
}

inline const BigComplex operator-(const double& left, const BigComplex& right)
{
    return BigComplex(left*pow(10,-right.m_e)-right.m_c,-right.m_e);
}

inline const BigComplex operator*(const double& left, const BigComplex& right)
{
    return BigComplex(right.m_c*left,right.m_e);
}

inline const BigComplex operator/(const double& left, const BigComplex& right)
{
    return BigComplex(left/right.m_c,-right.m_e);
}

inline bool operator==(const BigComplex& left, const BigComplex& right)
{
    return left.m_c==right.m_c && left.m_e==right.m_e;
}

inline bool operator==(const std::complex<double>& left, const BigComplex& right)
{
    return left==std::complex<double>(right);
}

inline bool operator==(const BigComplex& left, const std::complex<double>& right)
{
    return std::complex<double>(left)==right;
}

inline bool operator==(const BigComplex& left, const double& right)
{
    return std::complex<double>(left)==right;
}

inline bool operator==(const double& left, const BigComplex& right)
{
    return left==std::complex<double>(right);
}

inline std::ostream& operator<<(std::ostream& os,const BigComplex& bc)
{
    os<<bc.m_c<<"e"<<bc.m_e;
    return os;
}

#endif//_BIGCOMPLEX_H
