//
//  integrate_const.h
//  OdeSolver
//
//  fuehrt Integration des ODE-Systems mit konstanter Schrittweite aus

#ifndef OdeSolver_integrate_const_h
#define OdeSolver_integrate_const_h

#include "Definitions.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>



template<class Stepper, class System, class State, class Time>
size_t integrate_const(Stepper stepper, System sys, State &x, Time startTime, Time endTime, Time timeStep,  bool consoleOutput = false, bool fileOutput = true){
    
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
    filename << "Solution" << stepper.stepperName << ".dat";
    
    //testen, ob Datei schon existiert
    int fileCount = 1;
    while (file_exists(filename.str())) {
        filename.str("");
        filename << "Solution" << stepper.stepperName << "(" << fileCount << ").dat";
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
        
        //Schritt ausfuehren
        //bekommt Functor sys, aktuellen Zustandswert xCurrent, Speichervariable fuer neuen Zustandswert xNext, aktuelle Zeit tCurrent, die Schrittweite timeStep
        stepper.doStep(sys, xCurrent, xNext, tCurrent, timeStep );
        
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
