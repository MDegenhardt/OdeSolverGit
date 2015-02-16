//
//  Euler.h
//  OdeSolver
//
//  Explizites Euler-Verfahren 1.Ordnung zur Ausfuehrung des Integrationsschrittes

#ifndef OdeSolver_Euler_h
#define OdeSolver_Euler_h
#include <iostream>


class Euler{
public:
    
    //Konstruktor
    Euler(){
        stepperName = "Euler";
    }
    Euler(std::string solverName){
        stepperName = solverName;
    }
    
    //Copy-Konstruktor
    Euler(const Euler& otherEuler){
        stepperName = otherEuler.stepperName;
    }
    
    //Destruktor
    ~Euler(){
    }
    
    //Name des Steppers
    std::string stepperName;
    
    //fuehrt den Euler-Schritt aus => 1 Funktionsauswertung
    //bekommt Functor sys, aktuellen Zustandswert x_current, Speichervariable fuer neuen Zustandswert x_new, aktuelle Zeit t, die Schrittweite dt
    template<class System, typename State, typename Value = double, typename Deriv = State, typename Time = Value>
    void doStep(System sys, State x_current, Deriv &x_new, Time t, Time dt){
        
        State vectorEvaluated(x_current.size());
        //Funktionsauswertung
        sys(x_current, vectorEvaluated, t);
        
        for (int i=0; i<x_current.size(); i++) {
            x_new[i] = x_current[i] + vectorEvaluated[i]*dt;
        }


    }
    
    
    
};
#endif
