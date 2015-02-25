//
//  ExampleFunctors.h
//  OdeSolver
//
//  Beispiel-Gleichungssysteme zur Zeitintegration


#ifndef OdeSolver_ExampleFunctors_h
#define OdeSolver_ExampleFunctors_h

#include "Definitions.h"
#include <cmath>

// x'''(t) +x''(t)-4x'(t)-4x(t) = exp(-t)+sin(t)
template <class T>
class third_order_system{
public:
    void operator()(const std::vector<T>& x, std::vector<T>& dxdt, const T t ){
        //x'(t) = RHS(x,t)
        dxdt[0] = x[1];
        dxdt[1] = x[2];
        dxdt[2] = -x[2]+4*x[1]+4*x[0]+exp(-1*t)+sin(t);
    }
};

template <class T>
class harm_osc {
    double m_gam;
public:
    harm_osc( double gam ) : m_gam(gam) { }
    
    void operator() ( const std::vector<T> &x , std::vector<T> &dxdt , const T /* t */ )
    {
        dxdt[0] = x[1];
        dxdt[1] = -x[0] - m_gam*x[1];
    }
};

template <class T>
class first_order_system{
public:
    void operator()(const std::vector<T> &x, std::vector<T> &dxdt, const T t){
        dxdt[0] = exp(t);
    }
};

template <class T>
class first_order_system_complex{
public:
    void operator()(const std::vector<std::complex<T> > &x, std::vector<std::complex<T> > &dxdt, const T /* t */){
        dxdt[0] = x[0];
    }
};

template <class T>
class second_order_system_complex{
public:
    void operator()(const std::vector<std::complex<T> > &x, std::vector<std::complex<T> > &dxdt, const T /* t */){
        dxdt[0] = x[1];
        dxdt[1] = -x[0];
    }
};

template <class T>
class third_order_system_valarray{
public:
    void operator()(const std::valarray<T>& x, std::valarray<T>& dxdt, const T t ){
        dxdt[0] = x[1];
        dxdt[1] = x[2];
        dxdt[2] = -x[2]+4*x[1]+4*x[0]+exp(-1*t)+sin(t);
    }
};


#endif
