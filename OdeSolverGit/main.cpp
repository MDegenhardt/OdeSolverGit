//
//  main.cpp
//  OdeSolver
//


#include "OdeSolver.h"


int main(int argc, const char * argv[])
{

    // *** Strahl im Gauss-Feld ***
    //Start-Bedingungen: [px0, py0, x0, y0]
    state_type_double rayState(4);
    rayState[0] = 1.0;
    rayState[1] = 0.0;
    rayState[2] = -1.0;
    rayState[3] = 1.0;
    
    //Namen fuer Datei definieren
    std::stringstream nameEuler;
    nameEuler << "GaussEulerP" << rayState[0] << "_" << rayState[1]
    << "X" << rayState[2]<< "_" << rayState[3];
    std::stringstream nameRK4K;
    nameRK4K << "GaussRK4P" << rayState[0] << "_" << rayState[1]
    << "X" << rayState[2]<< "_" << rayState[3] << "K";
    std::stringstream nameRK4A;
    nameRK4A << "GaussRK4P" << rayState[0] << "_" << rayState[1]
    << "X" << rayState[2]<< "_" << rayState[3] << "A";
    
    //Functor fuer Gauss-Feld erstellen
    fourth_order_system_gauss<double> gaussField;
    //Stepper initialisieren
    Euler eg(nameEuler.str());
    RungeKutta4 rk4gK(nameRK4K.str());
    RungeKutta4 rk4gA(nameRK4A.str());
    // Strahl durch Integration berechnen
    size_t stepsGaussK = integrate_const(rk4gK, gaussField, rayState, 0.0, 7.1, 0.5, false, false);
    
    //Anfangsbedingungen zuruecksetzen
    rayState[0] = 1.0;
    rayState[1] = 0.0;
    rayState[2] = -1.0;
    rayState[3] = 1.0;
    // Strahl durch Integration berechnen
    size_t stepsGaussA = integrate_adaptive(rk4gA, gaussField, rayState, 0.0, 7.1, 0.5, false, false);
    
    // *** 3D Strahl im Gauss-Feld ***
    //Start-Bedingungen: [px0, py0, pz0, x0, y0, z0]
    state_type_double rayState3d(6);
    rayState3d[0] = 0.0;
    rayState3d[1] = 0.0;
    rayState3d[2] = 1.0;
    rayState3d[3] = 1.0;
    rayState3d[4] = 1.0;
    rayState3d[5] = -1.0;
    
    std::stringstream nameRK4K3d;
    nameRK4K3d << "GaussRK43d_P" << rayState3d[0] << "_" << rayState3d[1] << "_" << rayState3d[2] << "X" << rayState3d[3] << "_" << rayState3d[4] << "_" << rayState3d[5] << "K";

    //Functor fuer Gauss-Feld erstellen
    fourth_order_system_gauss_3d<double> gaussField3d;
    //Stepper initialisieren
    RungeKutta4 rk4gK3d(nameRK4K3d.str());
    // Strahl durch Integration berechnen
    size_t stepsGaussK3d = integrate_const(rk4gK3d, gaussField3d, rayState3d, 0.0, 5.0, 0.5, false, false);
    
    
    // *** Strahl im Gauss-Feld mit Schnittpunktberechnung ***
    
    // 1. adaptives Vorgehen
    //Start-Bedingungen: [px0, py0, x0, y0]
    state_type_double rayRDState(4);
    rayRDState[0] = 1.0;
    rayRDState[1] = 0.0;
    rayRDState[2] = -1.0;
    rayRDState[3] = 1.0;
    
    //Radius muss selben Typ wie State-Type-Elemente haben
    Raindrop rd1(1.5);
    //Name fuer Datei definieren
    std::stringstream nameRaindrop;
    nameRaindrop << "RK4RaindropP" << rayRDState[0] << "_" << rayRDState[1]
    << "X" << rayRDState[2]<< "_" << rayRDState[3] << "_rek";
    //Stepper initialisieren
    RungeKutta4 rk4rd(nameRaindrop.str());
    
    //Berechnung des Schnittpunktes mit dem Rand durch rekursives Vorgehen mit Halbierung der Schrittweite
//    size_t stepsRdA = rd1.calcBoundaryIntersectionAdaptive(rk4rd, gaussField, rayRDState, 0.5, false);
    
    // 2. Interpolation+Newton
    //Start-Bedingungen: [px0, py0, x0, y0]
    rayRDState[0] = 1.0;
    rayRDState[1] = 0.0;
    rayRDState[2] = -1.0;
    rayRDState[3] = 1.0;
    
    //Name fuer Datei definieren
    std::stringstream nameRaindropHermite;
    nameRaindropHermite << "RK4RaindropP" << rayRDState[0] << "_" << rayRDState[1]
    << "X" << rayRDState[2]<< "_" << rayRDState[3] << "_hermite";
    
    //Stepper initialisieren
    RungeKutta4 rk4rdh(nameRaindropHermite.str());
    
    //Berechnung des Schnittpunktes mit dem Rand durch Interpolation mit Hermite-Spline und Newton-Verfahren mit AD
//    size_t stepsRdH = rd1.calcBoundaryIntersectionHermite(rk4rdh, gaussField, rayRDState, 1.0, false);

    
    //sonstige Beispiele aufrufen
    runExamples();
    

    
    runGaussPicture();
    
    
    
    
    
    return 0;
}

