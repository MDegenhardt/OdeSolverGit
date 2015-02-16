//
//  RungeKutta4.h
//  OdeSolver
//
//  Klassisches Runge-Kutta-Verfahren 4.Ordnung zur Ausfuehrung des Integrationsschrittes

#ifndef OdeSolver_RungeKutta4_h
#define OdeSolver_RungeKutta4_h
#include <iostream>

class RungeKutta4{
public:
    
    //Konstruktor
    RungeKutta4(){
        stepperName = "RungeKutta4";
    }
    RungeKutta4(std::string solverName){
        stepperName = solverName;
    }
    
    //Copy-Konstruktor
    RungeKutta4(const RungeKutta4& otherRK4){
        stepperName = otherRK4.stepperName;
    }
    
    //Destruktor
    ~RungeKutta4(){
    }
    
    //Name des Steppers
    std::string stepperName;
    
    //fuehrt den RK4-Schritt aus => 4 Funktionsauswertungen
    //bekommt Functor sys, aktuellen Zustandswert x_current, Speichervariable fuer neuen Zustandswert x_new, aktuelle Zeit t, die Schrittweite dt
    template<class System, typename State, typename Value = double, typename Deriv = State, typename Time = Value>
    void doStep(System sys, State x_current, Deriv &x_new, Time t, Time dt){
        
        const double a21=0.5, a31=0.0, a41=0.0, a32=0.5, a42=0.0, a43=1.0,
        b1=1.0/6.0, b2=1.0/3.0, b3=1.0/3.0, b4=1.0/6.0, c1=0.0, c2=0.5, c3=0.5, c4=1.0;
        
        size_t state_vec_size = x_current.size();
        //speichern des naechsten x-Werte-Vektors
        State x_temp(state_vec_size);
        
        //RK4 Koeffizienten
        State k1(state_vec_size), k2(state_vec_size), k3(state_vec_size), k4(state_vec_size);
        
        //k1 berechnen
        //Funktionsauswertung bei x_current, t+c1*dt
        sys(x_current, k1, t+c1*dt);
        
        //k2 berechnen
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*a21*k1[i];
        }
        //Funktionsauswertung bei x_temp, t+c2*dt
        sys(x_temp, k2, t+c2*dt);
        
        
        //k3
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a31*k1[i]+a32*k2[i]);
        }
        //Funktionsauswertung bei x_temp,t+c3*dt
        sys(x_temp, k3, t+c3*dt);
        
        
        //k4
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a41*k1[i]+a42*k2[i]+ a43*k3[i]);
        }
        //Funktionsauswertung bei x_temp,t+c4*dt
        sys(x_temp, k4, t+c4*dt);
        
        
        //RK4-Schritt
        for (int i=0; i<x_current.size(); i++) {
            x_new[i] = x_current[i] + dt*(b1*k1[i]+b2*k2[i]+b3*k3[i]+b4*k4[i]);
        }
        
    }
    
};

#endif