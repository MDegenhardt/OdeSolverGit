//
//  NonlinearFunctors.h
//  OdeSolver
//
//  Beispiel-Functoren fuer das Newton-Verfahren (jeweils Funktionen und Ableitungen)

#ifndef OdeSolver_NonlinearFunctors_h
#define OdeSolver_NonlinearFunctors_h

#include "Definitions.h"
#include <cmath>

template<class T, class U>
class nonlinear_system1{
public:
    void operator()(const std::vector<T>& x, std::vector<U>& fx){
        //f(x)=1-2/x^2
        fx[0] = 1-2/(x[0]*x[0]);
        //f(x)=x^2-612
        fx[1] = (x[1]*x[1]) -612;
        //f(x)=x^2-1
        fx[2] = (x[2]*x[2]) -1;
    }
};
//template<class T, class U>
//class nonlinear_system1derivate{
//public:
//    void operator()(const std::vector<T>& x, std::vector<U>& fx){
//        //f'(x)=4/x^3
//        fx[0] = 4/(x[0]*x[0]*x[0]);
//        //f'(x)=2x
//        fx[1] = 2*x[1];
//        //f'(x)=2x
//        fx[2] = 2*x[2];
//    }
//};

template<class T, class U>
class nonlinear_function1{
public:
    void operator()(const std::vector<T>& x, std::vector<U>& fx){
        //f(x)=1-x^2
        fx[0] = 2-(x[0]*x[0]);
    }
};

//class nonlinear_function1derivate{
//public:
//    void operator()(const state_type_double& x, state_type_double& fx){
//        //f'(x)=-2x
//        fx[0] = -2.0*x[0];
//    }
//};

template <class T>
class complex_nonlinear_system1{
public:
    void operator()(const std::vector<std::complex<T>>& x, std::vector<std::complex<T>>& fx){
        //f(x)=1-2/x^2
        fx[0] = cdouble(1)-cdouble(2)/(x[0]*x[0]);
        //f(x)=x^2-612
        fx[1] = (x[1]*x[1]) -cdouble(612);
        //f(x)=x^2+1
        fx[2] = (x[2]*x[2]) +cdouble(1);
    }
};

template <class T>
class complex_nonlinear_system1derivate{
public:
    void operator()(const std::vector<std::complex<T>>& x, std::vector<std::complex<T>>& fx){
        //f'(x)=4/x^3
        fx[0] = cdouble(4)/(x[0]*x[0]*x[0]);
        //f'(x)=2x
        fx[1] = cdouble(2)*x[1];
        //f'(x)=2x
        fx[2] = cdouble(2)*x[2];
    }
};

template <class T>
class implicit_system1{
public:
    void operator()(const std::vector<T>& x, std::vector<T>& dx, T t){
        //x'(t) = RHS(x,t)
        //x'(t) = exp(t)
        dx[0] = exp(t);
        //x'(t) = 20t+t^2
        dx[1] = 20.0*t+t*t;
    }
};

template <class T>
class implicit_system1derivate{
public:
    void operator()(const std::vector<T>& x, std::vector<T>& dx, T t){
        //x''(t) = exp(t)
        dx[0] = exp(t);
        //x''(t)=20+2t
        dx[1] = 20.0+2*t;

    }
};

template<class T, class U>
class implicit_system1_ad{
public:
    void operator()(const std::vector<T>& x, std::vector<U>& dx, double t){
        //x'(t) = RHS(x,t)
        //x'(t) = exp(t)
        dx[0] = exp(t);
        //x'(t) = 20t+t^2
        dx[1] = 20.0*t+t*t;
    }
};


#endif
