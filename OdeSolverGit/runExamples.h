//
//  runExamples.h
//  OdeSolver
//
//  Funktion um sonstige Beispiele zusammenzufassen, kann in main.cpp aufgerufen werden

#ifndef OdeSolver_runExamples_h
#define OdeSolver_runExamples_h

#include "RungeKutta4.h"
#include "Euler.h"
#include "Dopri5.h"
#include "Fehlberg78.h"
#include "RungeKutta14.h"
#include "implicit_euler.h"
#include "implicit_euler_ad.h"

void runExamples(){
    
    RungeKutta4 rk4;
    Euler e;
    Dopri5 dp5;
    Fehlberg78 fb78;

    //RK4 Table
    ButcherTable but1;
    but1.a21 = 0.5, but1.a31 = 0.0, but1.a32=0.5, but1.a43=1.0;
    but1.c1 = 0.0, but1.c2=0.5, but1.c3=0.5, but1.c4=1.0;
    but1.b1=1.0/6.0, but1.b2 = 1.0/3.0, but1.b3 = 1.0/3.0,but1.b4 = 1.0/6.0;

    RungeKutta14 rk14(but1);

    
    //System 3.Ordnung
    state_type_double x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    
    third_order_system<double> tos;
    size_t steps = integrate_const(fb78, tos, x, 0.0, 5.0, 0.1, false, false);
    
    //System 2.Ordnung
    state_type_double x2(2);
    x2[0] = 1.0;
    x2[1] = 0.0;
    
    harm_osc<double> hos(0.15);
    size_t steps2 = integrate_const(e, hos, x2, 0.0, 10.0, 0.1, false, false);
    
    //System 1.Ordnung
    state_type_double d0 = {1};
    first_order_system<double> fos;
    e.stepperName = "Euler_exp(t)0_5";
    size_t steps3 = integrate_const(e, fos, d0, 0.0, 5.01, 0.5, false, false);
    e.stepperName = "Euler";
    
    d0 = {1};
    rk4.stepperName ="RK4_exp(t)0_5";
    steps3 = integrate_const(rk4, fos, d0, 0.0, 5.01, 0.5, false, false);
    rk4.stepperName = "RungeKutta4";
    
    //System 2.Ordnung, Complex
    std::complex<double> c1(1.0, 2.0);
    std::complex<double> c2(3.0, 4.0);
    state_type_complex c(2);
    c[0] = c1;
    c[1] = c2;
    
    second_order_system_complex<double> fosc;
    size_t stepsC = integrate_const(e, fosc, c, 0.0, 10.0, 0.1, false, false);
    
    //Newton complex 3 Gleichungen
    state_type_complex nc(3);
    cdouble nc1(2.0, 0.0);
    cdouble nc2(10.0, 0.0);
    cdouble nc3(1.0,1.0);
    nc[0] = nc1;
    nc[1] = nc2;
    nc[2] = nc3;
    complex_nonlinear_system1<double> cfx;
    complex_nonlinear_system1derivate<double> cfdxdt;
    size_t maxIterations = newton_raphson(cfx, cfdxdt, nc, false);
    
    //Newton 3 Gleichungen
    state_type_double n(3);
    n[0] = 1;
    n[1] = 10;
    n[2] = 1;
    nonlinear_system1<RallNo<double>, RallNo<double> > fx;
    size_t  cmaxIterations = newton_raphson_ad(fx, n, false);
    
    
    //Newton 1 Gleichung
    state_type_double n1(1);
    n1[0] = 1.0;
    nonlinear_function1<RallNo<double>, RallNo<double> > fx1;
    size_t maxIterations1 = newton_raphson_ad(fx1, n1, false, 1000, 1e-8);
    
    //Valarray
    const double vals[] = {1.0, 2.0, 3.0};
    val_array va(vals, 3);
    
    third_order_system_valarray<double> tosv;
    size_t stepsVa = integrate_const(e, tosv, va, 0.0, 10.0, 0.1, false, false);
    
    
    //System 2.Ordnung, Implizit
    state_type_double xImp2(2);
    xImp2[0] = 1.0;
    xImp2[1] = 1.0;
    
    implicit_system1<double> fxI;
    implicit_system1derivate<double> fxdxI;
    
    implicit_system1_ad<RallNo<double>, RallNo<double> > fxIad;
    
    size_t stepsImp = implicit_euler(fxI, fxdxI, xImp2, 0.0, 5.0, 0.001, false, false);
    
    xImp2[0] = 1.0;
    xImp2[1] = 1.0;
    
    stepsImp = implicit_euler_ad(fxIad, xImp2, 0.0, 5.0, 0.001, false, false);
    
    
    
}

#endif
