//
//  Dopri5.h
//  OdeSolver
//
//  Verfahren nach Dormand-Prince 5.Ordnung zur Ausfuehrung des Integrationsschrittes


#ifndef OdeSolver_Dopri5_h
#define OdeSolver_Dopri5_h


class Dopri5{
public:
    
    //Konstruktor
    Dopri5(){
        stepperName = "Dopri5";
    }
    Dopri5(std::string solverName){
        stepperName = solverName;
    }
    
    //Copy-Konstruktor
    Dopri5(const Dopri5& otherDopri5){
        stepperName = otherDopri5.stepperName;
    }
    
    //Destruktor
    ~Dopri5(){
    }
    
    //Name des Steppers
    std::string stepperName;
    
    //fuehrt den Dopri5-Schritt aus => 6 Funktionsauswertungen
    //bekommt Functor sys, aktuellen Zustandswert x_current, Speichervariable fuer neuen Zustandswert x_new, aktuelle Zeit t, die Schrittweite dt
    template<class System, typename State, typename Value = double, typename Deriv = State, typename Time = Value>
    void doStep(System sys, State x_current, Deriv &x_new, Time t, Time dt){
        
        
        const double
        a21=0.2,
        a31=3.0/40.0, a32=9.0/40.0,
        a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
        a51=19372.0/6561.0, a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
        a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
        
        b1=35.0/384.0, b2=0.0, b3=500.0/1113.0, b4=125.0/192.0, b5=-2187.0/6784.0, b6=11.0/84.0,
        c1=0.0, c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0, c6=1.0;
        
        size_t state_vec_size = x_current.size();
        //speichern des naechsten x-Werte-Vektors
        State x_temp(state_vec_size);
        
        //Dopri5 Koeffizienten
        State k1(state_vec_size), k2(state_vec_size), k3(state_vec_size), k4(state_vec_size), k5(state_vec_size), k6(state_vec_size);
        
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
            x_temp[i] = x_current[i] + dt*(a41*k1[i]+a42*k2[i]+a43*k3[i]);
        }
        //Funktionsauswertung bei x_temp,t+c4*dt
        sys(x_temp, k4, t+c4*dt);
        
        //k5
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
        }
        //Funktionsauswertung bei x_temp,t+c5*dt
        sys(x_temp, k5, t+c5*dt);
        
        //k6
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a61*k1[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
        }
        //Funktionsauswertung bei x_temp,t+c6*dt
        sys(x_temp, k6, t+c6*dt);
        
        
        //Dopri5-Schritt
        for (int i=0; i<x_current.size(); i++) {
            x_new[i] = x_current[i] + dt*(b1*k1[i]+b2*k2[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]);
        }
        
    }
    
};



#endif
