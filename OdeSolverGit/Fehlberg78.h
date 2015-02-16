//
//  Fehlberg78.h
//  OdeSolver
//
//  Verfahren nach Fehlberg 7. bzw. 8. Ordnung zur Ausfuehrung des Integrationsschrittes


#ifndef OdeSolver_Fehlberg78_h
#define OdeSolver_Fehlberg78_h

class Fehlberg78{
public:
    
    //Konstruktor
    Fehlberg78(){
        stepperName = "Fehlberg78";
    }
    Fehlberg78(std::string solverName){
        stepperName = solverName;
    }
    
    //Copy-Konstruktor
    Fehlberg78(const Fehlberg78& otherFB78){
        stepperName = otherFB78.stepperName;
    }
    
    //Destruktor
    ~Fehlberg78(){
    }
    
    //Name des Steppers
    std::string stepperName;
    
    //fuehrt den Fehlberg78-Schritt aus => 13 Funktionsauswertungen
    //bekommt Functor sys, aktuellen Zustandswert x_current, Speichervariable fuer neuen Zustandswert x_new, aktuelle Zeit t, die Schrittweite dt
    template<class System, typename State, typename Value = double, typename Deriv = State, typename Time = Value>
    void doStep(System sys, State x_current, Deriv &x_new, Time t, Time dt){
        
        
        const double
        a21=2.0/27.0,
        a31=1.0/36.0, a32=1.0/12.0,
        a41=1.0/24.0, a42=0.0, a43=1.0/8.0,
        a51=5.0/12.0, a52=0.0, a53=-25.0/16.0, a54=25.0/16.0,
        a61=1.0/20.0, a62=0.0, a63=0.0, a64=1.0/4.0, a65=1.0/5.0,
        a71=-25.0/108.0, a72=0.0, a73=0.0, a74=125.0/108.0, a75=-65.0/27.0, a76=125.0/54.0,
        a81=31.0/300.0, a82=0.0, a83=0.0, a84=0.0, a85=61.0/225.0, a86=-2.0/9.0, a87=13.0/900.0,
        a91=2.0, a92=0.0, a93=0.0, a94=-53.0/6.0, a95=704.0/45.00, a96=-107.0/9.0, a97=67.0/90.0, a98=3.0,
        a101=-91.0/108.0, a102=0.0, a103=0.0, a104=23.0/108.0, a105=-976.0/135.0, a106=311.0/54.0, a107=-19.0/60.0, a108=17.0/6.0, a109=-1.0/12.0,
        a111=2383.0/4100.0, a112=0.0, a113=0.0, a114=-341.0/164.0, a115=4496.0/1025.0, a116=-301.0/82.0, a117=2133.0/4100.0, a118=45.0/82.0, a119=45.0/164.0, a1110=18.0/41.0,
        a121=3.0/205.0, a122=0.0, a123=0.0, a124=0.0, a125=0.0, a126=-6.0/41.0, a127=-3.0/205.0, a128=-3.0/41.0, a129=3.0/41.0, a1210=6.0/41.0, a1211=0.0,
        a131=-1777.0/4100.0, a132=0.0, a133=0.0, a134=-341.0/164.0, a135=4496.0/1025.0, a136=-289.0/82.0, a137=2193.0/4100.0, a138=51.0/82.0, a139=33.0/164.0, a1310=12.0/41.0, a1311=0.0, a1312=1.0,
        
        //7.Ordnung
//        b1=41.0/840.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=34.0/105.0, b7=9.0/35.0, b8=9.0/35.0, b9=9.0/280.0, b10=9.0/280.0, b11=41.0/840.0, b12=0.0, b13=0.0,
        //8.Ordnung
        b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=34.0/105.0, b7=9.0/35.0, b8=9.0/35.0, b9=9.0/280.0, b10=9.0/280.0, b11=0.0, b12=41.0/840.0, b13=41.0/840.0,
        
        c1=0.0, c2=2.0/27.0, c3=1.0/9.0, c4=1.0/6.0, c5=5.0/12.0, c6=0.5, c7=5.0/6.0, c8=1.0/6.0, c9=2.0/3.0, c10=1.0/3.0, c11=1.0, c12=0.0, c13=1.0;
        
        size_t state_vec_size = x_current.size();
        //speichern des naechsten x-Werte-Vektors
        State x_temp(state_vec_size);
        
        //Fehlberg Koeffizienten
        State k1(state_vec_size), k2(state_vec_size), k3(state_vec_size), k4(state_vec_size), k5(state_vec_size),
        k6(state_vec_size), k7(state_vec_size), k8(state_vec_size), k9(state_vec_size), k10(state_vec_size),
        k11(state_vec_size), k12(state_vec_size), k13(state_vec_size);
        
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
        
        //k7
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a71*k1[i]+a72*k2[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
        }
        //Funktionsauswertung bei x_temp,t+c7*dt
        sys(x_temp, k7, t+c7*dt);
        
        //k8
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a81*k1[i]+a82*k2[i]+a83*k3[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
        }
        //Funktionsauswertung bei x_temp,t+c8*dt
        sys(x_temp, k8, t+c8*dt);
        
        //k9
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a91*k1[i]+a92*k2[i]+a93*k3[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]);
        }
        //Funktionsauswertung bei x_temp,t+c9*dt
        sys(x_temp, k9, t+c9*dt);
        
        //k10
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a101*k1[i]+a102*k2[i]+a103*k3[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i]+a109*k9[i]);
        }
        //Funktionsauswertung bei x_temp,t+c10*dt
        sys(x_temp, k10, t+c10*dt);
        
        //k11
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a111*k1[i]+a112*k2[i]+a113*k3[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);
        }
        //Funktionsauswertung bei x_temp,t+c11*dt
        sys(x_temp, k11, t+c11*dt);
        
        //k12
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a121*k1[i]+a122*k2[i]+a123*k3[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]+a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k11[i]);
        }
        //Funktionsauswertung bei x_temp,t+c12*dt
        sys(x_temp, k12, t+c12*dt);
        
        //k13
        //naechsten x-Wert-Vektor x_temp berechnen
        for (int i=0; i<state_vec_size; i++) {
            x_temp[i] = x_current[i] + dt*(a131*k1[i]+a132*k2[i]+a133*k3[i]+a134*k4[i]+a135*k5[i]+a136*k6[i]+a137*k7[i]+a138*k8[i]+a139*k9[i]+a1310*k10[i]+a1311*k11[i]+a1312*k12[i]);
        }
        //Funktionsauswertung bei x_temp,t+c13*dt
        sys(x_temp, k13, t+c13*dt);
        
        
        //Fehlberg78-Schritt
        for (int i=0; i<x_current.size(); i++) {
            x_new[i] = x_current[i] + dt*(b1*k1[i]+b2*k2[i]+b3*k3[i]+b4*k4[i]+b5*k5[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k11[i]+b12*k12[i]+b13*k13[i]);
        }

        
    }
    
};


#endif
