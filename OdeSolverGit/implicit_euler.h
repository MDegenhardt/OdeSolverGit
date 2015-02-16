//
//  implicit_euler.h
//  OdeSolver
//
//  Implizites Euler-Verfahren 1.Ordnung zur Ausfuehrung des Integrationsschrittes


#ifndef OdeSolver_implicit_euler_h
#define OdeSolver_implicit_euler_h

#include "Definitions.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>
    

template<class System, class System2, class State, class Time>
size_t implicit_euler(System sys, System2 deriv, State &x, Time startTime, Time endTime, Time timeStep,  bool consoleOutput = false, bool fileOutput = true){

    Time tCurrent = startTime;
    Time tNext = startTime; // = startTime vector initialisieren (Werte egal)
    State xCurrent = x;
    State xNext = x; //initialisieren
    
    size_t iterationCount = 0;
    size_t state_vec_size = x.size();
    
    //Output
    std::cout.precision(6);
    std::cout.setf( std::ios::showpoint | std::ios::scientific);
    //Consolen-Ausgabe Startwerte
    if (consoleOutput) {
        std::cout << "T(" << iterationCount << "): " << tCurrent << ", X(" << iterationCount << "): ";
        
        for (size_t i=0; i<state_vec_size; i++) {
            std::cout << xCurrent[i] << ", ";
        }
        std::cout << "\n";
    }
    //File-Output
    std::stringstream filename;
    filename << "SolutionImplicitEuler.dat";
    
    //testen, ob Datei schon existiert
    int fileCount = 1;
    while (file_exists(filename.str())) {
        filename.str("");
        filename << "SolutionImplicitEuler(" << fileCount << ").dat";
        fileCount++;
    }
    std::string fname = filename.str();
    std::ofstream write_file(fname);
    
    const int PRECISION = 7;
    const int WIDTH     = 20;
    if (!fileOutput) {
        //file wieder loeschen
        write_file.close();
        std::remove ( fname.c_str() );
    }
    if (fileOutput){
        write_file.precision(PRECISION);
        write_file.setf(std::ios::scientific | std::ios::showpoint);
        
        //Startwerte schreiben
        write_file << tCurrent << std::setw(WIDTH);
        
        for (size_t i=0; i<state_vec_size; i++) {
            write_file << xCurrent[i] << std::setw(WIDTH);
        }
        
        write_file << "\n";
    }
    // end Output
    
    do{
        iterationCount++;
        //naechsten Zeitpunkt berechnen
        tNext = tCurrent + timeStep;
        
        // impliziter Euler-Schritt
        //naechsten Wert durch Anwendung des Newton-Verfahrens berechnen
        double eps = 1e-6;
        size_t maxIterations = 50;
        
        for (int i=0; i<state_vec_size; i++) {
            
            double error = 0.0;
            size_t iterationsNewton = 0;
            
            State functionEvaluated(state_vec_size);
            State derivateEvaluated(state_vec_size);
            
            //Startwert fuer Newton
            xNext[i] = xCurrent[i];
            
            do{
                iterationsNewton++;
                
                //F(x) auswerten an xNext und  tNext
                sys(xNext, functionEvaluated, tNext);
                //F'(x) auswerten
                deriv(xNext, derivateEvaluated, tNext);
                
                //Newton-Schritt mit Formel des impliziten Euler-Verfahrens
                //x_{n+1} = x_{n+1} - F(x_{n+1})/F'(x_{n+1}), F(x_{n+1}) = x_{n+1} - x_n - dt*f(t_{n+1}, x_{n+1}), F'(x_{n+1}) = 1 - dt*f'(t_{n+1}, x_{n+1})
                xNext[i]    = xNext[i] - (xNext[i]-xCurrent[i]-timeStep*functionEvaluated[i])/(1.0 - timeStep*derivateEvaluated[i]);
                error = std::abs(xCurrent[i]-xNext[i]);
                
            } while ( error > eps && (iterationsNewton < maxIterations));
            
        }
        // Ende Euler-Schritt

        //Zeitpunkt und Werte fuer naechste Iteration weitersetzen
        tCurrent = tNext;
        xCurrent = xNext;
        
        //Output
        //Consolen-Ausgabe
        if (consoleOutput) {
            std::cout << "T(" << iterationCount << "): " << tCurrent << ", X(" << iterationCount << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << xCurrent[i] << ", ";
            }
            
            std::cout << "\n";
            
        }
        //File-Output
        if (fileOutput){
            write_file << tCurrent << std::setw(WIDTH);
            
            for (size_t i=0; i<state_vec_size; i++) {
                write_file << xCurrent[i] << std::setw(WIDTH);
            }
            
            write_file << "\n";
            
        }
        //end Output
        
        
    } while (tNext < endTime-timeStep);
    x = xNext;
    
    return iterationCount;
    
}


#endif
