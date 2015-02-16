//
//  RungeKutta14.h
//  OdeSolver
//
//  Allgemeines explizites Runge-Kuuta-Verfahren zur Ausfuehrung des Integrationsschrittes, kann bis zu 14 Funktionsauswertungen durchfuehren


#ifndef OdeSolver_RungeKutta14_h
#define OdeSolver_RungeKutta14_h

#include "ButcherTable.h"

class RungeKutta14{
public:
    
    //Konstruktor
    RungeKutta14(ButcherTable bTable){
        stepperName = "RK14";
        bt = bTable;
    }
    RungeKutta14(ButcherTable bTable, std::string solverName){
        stepperName = solverName;
        bt = bTable;
    }
    
    //Copy-Konstruktor
    RungeKutta14(const RungeKutta14& otherRK14){
        stepperName = otherRK14.stepperName;
        bt = otherRK14.bt;
    }
    
    //Destruktor
    ~RungeKutta14(){
    }
    
    //Name des Steppers
    std::string stepperName;
    //Butcher-Table des RK14 Verfahrens
    ButcherTable bt;
    
    //fuehrt den RK14-Schritt aus
    //bekommt Functor sys, aktuellen Zustandswert x_current, Speichervariable fuer neuen Zustandswert x_new, aktuelle Zeit t, die Schrittweite dt
    template<class System, typename State, typename Value = double, typename Deriv = State, typename Time = Value>
    void doStep(System sys, State x_current, Deriv &x_new, Time t, Time dt){
        
        
        size_t state_vec_size = x_current.size();
        //speichern des naechsten x-Werte-Vektors
        State x_temp(state_vec_size);
        
        //RK14 Koeffizienten
        State k1(state_vec_size), k2(state_vec_size), k3(state_vec_size), k4(state_vec_size), k5(state_vec_size),
        k6(state_vec_size), k7(state_vec_size), k8(state_vec_size), k9(state_vec_size), k10(state_vec_size),
        k11(state_vec_size), k12(state_vec_size), k13(state_vec_size), k14(state_vec_size);
        
        //k1 berechnen
        //Funktionsauswertung bei x_current, t+c1*dt
        sys(x_current, k1, t+bt.c1*dt);
        
        //k2 berechnen
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*bt.a21*k1[i];
        }
        //Funktionsauswertung bei x_temp, t+c2*dt
        sys(x_temp, k2, t+bt.c2*dt);
        
        //k3
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a31*k1[i]+bt.a32*k2[i]);
        }
        //Funktionsauswertung bei x_temp,t+c3*dt
        sys(x_temp, k3, t+bt.c3*dt);
        
        //k4
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a41*k1[i]+bt.a42*k2[i]+bt.a43*k3[i]);
        }
        //Funktionsauswertung bei x_temp,t+c4*dt
        sys(x_temp, k4, t+bt.c4*dt);
        
        //k5
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a51*k1[i]+bt.a52*k2[i]+bt.a53*k3[i]+bt.a54*k4[i]);
        }
        //Funktionsauswertung bei x_temp,t+c5*dt
        sys(x_temp, k5, t+bt.c5*dt);
        
        //k6
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a61*k1[i]+bt.a62*k2[i]+bt.a63*k3[i]+bt.a64*k4[i]+bt.a65*k5[i]);
        }
        //Funktionsauswertung bei x_temp,t+c6*dt
        sys(x_temp, k6, t+bt.c6*dt);
        
        //k7
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a71*k1[i]+bt.a72*k2[i]+bt.a73*k3[i]+bt.a74*k4[i]+bt.a75*k5[i]+bt.a76*k6[i]);
        }
        //Funktionsauswertung bei x_temp,t+c7*dt
        sys(x_temp, k7, t+bt.c7*dt);
        
        //k8
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a81*k1[i]+bt.a82*k2[i]+bt.a83*k3[i]+bt.a84*k4[i]+bt.a85*k5[i]+bt.a86*k6[i]+bt.a87*k7[i]);
        }
        //Funktionsauswertung bei x_temp,t+c8*dt
        sys(x_temp, k8, t+bt.c8*dt);
        
        //k9
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a91*k1[i]+bt.a92*k2[i]+bt.a93*k3[i]+bt.a94*k4[i]+bt.a95*k5[i]+bt.a96*k6[i]+bt.a97*k7[i]+bt.a98*k8[i]);
        }
        //Funktionsauswertung bei x_temp,t+c9*dt
        sys(x_temp, k9, t+bt.c9*dt);
        
        //k10
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a101*k1[i]+bt.a102*k2[i]+bt.a103*k3[i]+bt.a104*k4[i]+bt.a105*k5[i]+bt.a106*k6[i]+bt.a107*k7[i]+bt.a108*k8[i]+bt.a109*k9[i]);
        }
        //Funktionsauswertung bei x_temp,t+c10*dt
        sys(x_temp, k10, t+bt.c10*dt);
        
        //k11
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a111*k1[i]+bt.a112*k2[i]+bt.a113*k3[i]+bt.a114*k4[i]+bt.a115*k5[i]+bt.a116*k6[i]+bt.a117*k7[i]+bt.a118*k8[i]+bt.a119*k9[i]+bt.a1110*k10[i]);
        }
        //Funktionsauswertung bei x_temp,t+c11*dt
        sys(x_temp, k11, t+bt.c11*dt);
        
        //k12
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a121*k1[i]+bt.a122*k2[i]+bt.a123*k3[i]+bt.a124*k4[i]+bt.a125*k5[i]+bt.a126*k6[i]+bt.a127*k7[i]+bt.a128*k8[i]+bt.a129*k9[i]+bt.a1210*k10[i]+bt.a1211*k11[i]);
        }
        //Funktionsauswertung bei x_temp,t+c12*dt
        sys(x_temp, k12, t+bt.c12*dt);
        
        //k13
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a131*k1[i]+bt.a132*k2[i]+bt.a133*k3[i]+bt.a134*k4[i]+bt.a135*k5[i]+bt.a136*k6[i]+bt.a137*k7[i]+bt.a138*k8[i]+bt.a139*k9[i]+bt.a1310*k10[i]+bt.a1311*k11[i]+bt.a1312*k12[i]);
        }
        //Funktionsauswertung bei x_temp,t+c13*dt
        sys(x_temp, k13, t+bt.c13*dt);
        
        //k14
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(bt.a141*k1[i]+bt.a142*k2[i]+bt.a143*k3[i]+bt.a144*k4[i]+bt.a145*k5[i]+bt.a146*k6[i]+bt.a147*k7[i]+bt.a148*k8[i]+bt.a149*k9[i]+bt.a1410*k10[i]+bt.a1411*k11[i]+bt.a1412*k12[i]+bt.a1413*k13[i]);
        }
        //Funktionsauswertung bei x_temp,t+c14*dt
        sys(x_temp, k14, t+bt.c14*dt);
        
        
        //RK14-Schritt
        for (int i=0; i<x_current.size(); i++) {
            x_new[i] = x_current[i] + dt*(bt.b1*k1[i]+bt.b2*k2[i]+bt.b3*k3[i]+bt.b4*k4[i]+bt.b5*k5[i]+bt.b6*k6[i]+bt.b7*k7[i]+bt.b8*k8[i]+bt.b9*k9[i]+bt.b10*k10[i]+bt.b11*k11[i]+bt.b12*k12[i]+bt.b13*k13[i]+bt.b14*k14[i]);
        }
        
    }
    
};


#endif
