//
//  integrate_adaptive.h
//  OdeSolver
//
//  fuehrt Integration des ODE-Systems mit variabler Schrittweite aus


#ifndef OdeSolver_integrate_adaptive_h
#define OdeSolver_integrate_adaptive_h

#include "Definitions.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdio>



template<class Stepper, class System, class State, class Time>
size_t integrate_adaptive(Stepper stepper, System sys, State &x, Time startTime, Time endTime, Time timeStep,  bool consoleOutput = false, bool fileOutput = true){
    
    //*** Initialisierung ***
    
    Time tCurrent = startTime;
    Time tNext = startTime; // = startTime vector initialisieren (Werte egal)
    Time calculatedTimeStep = timeStep; //initialisieren
    
    State xCurrent = x;
    State xNextFull = x; //initialisieren
    State xNextHalf = x; //initialisieren
    
    size_t iterationCountAll = 0; //zaehlt alle Iterationen
    size_t iterationCountAcc = 0; //zahelt nur Iterationen mit akzeptierter SW
    size_t state_vec_size = x.size();
    
    //Output
    std::cout.precision(6);
    std::cout.setf( std::ios::showpoint | std::ios::scientific);
    //Consolen-Ausgabe Startwerte
    if (consoleOutput) {
        std::cout << "T(" << iterationCountAll << "): " << tCurrent << ", X(" << iterationCountAll << "): ";
        
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
    
    // *** Werte einstellen ***
    
    const int dbout = 0; //debug-output, Infos zu den Iterationen ein-/ausstellen
    const double SCALE = 2.0; //Aenderung der SW um max. diesen Faktor
    const double EPS = 5*1e-1; //absolute Genauigkeit
    const Time timeStepMIN = 1e-8; //minimal zul. SW
    const Time timeStepMAX = 5; //maximal zul. SW
    double trunc_error; //lokaler Diskretisierungsfehler
    
    // *** Iterationen durchfuehren ***
    
    do{
        iterationCountAll++;
        iterationCountAcc++;
        
        //*** Schritte ausfuehren ***
        
        //vollen Schritt ausfuehren, Wert in xNextFull speichern
        //bekommt Functor sys, aktuellen Zustandswert xCurrent, Speichervariable fuer neuen Zustandswert xNextFull, aktuelle Zeit tCurrent, die Schrittweite timeStep
        stepper.doStep(sys, xCurrent, xNextFull, tCurrent, timeStep);
        
        // 2 Halbschritte ausfuehren, Wert in xNextHalf speichern
        stepper.doStep(sys, xCurrent, xNextHalf, tCurrent, timeStep/2.0 );
        stepper.doStep(sys, xCurrent, xNextHalf, tCurrent, timeStep/2.0 );

        
        if (dbout) {
            std::cout << "T(" << iterationCountAll << "): " << tCurrent << ", XCurrent(" << iterationCountAll << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << xCurrent[i] << ", ";
            }
            std::cout << "\n";
            
            std::cout << "T(" << iterationCountAll << "): " << tCurrent << ", XFull(" << iterationCountAll << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << xNextFull[i] << ", ";
            }
            
            std::cout << "\n";
            
            std::cout << "T(" << iterationCountAll << "): " << tCurrent << ", X2Half(" << iterationCountAll << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << xNextHalf[i] << ", ";
            }
            
            std::cout << "\n";
            
        }
        
        // *** lokalen Diskretisierungsfehler berechnen mit Euklid-Norm ***
        double sum = 0.0;
        for (int i=0; i<state_vec_size; i++) {
            sum += pow(fabs(xNextFull[i] - xNextHalf[i]),2.0);
        }
        trunc_error = sqrt(sum);
        
        //Division durch 0 verhindern
        if (trunc_error == 0.0) {
            trunc_error = 1e-15;
        }
        
        if(dbout) std::cout << "trunc_error: " << trunc_error << "\n";
        
        // *** neuen timeStep berechnen ***
        
        //berechnete optimale SW
        calculatedTimeStep = timeStep * pow(fabs(EPS/ trunc_error), 0.2);
        
        if(dbout) std::cout << "Calc. TimeStep: " << calculatedTimeStep << "\n";
        
        //berechnete SW > als max. zulaessige SW in diesem Schritt => urspruengliche SW mit Skalierungsfaktor multiplizieren
        if (calculatedTimeStep / timeStep > SCALE) {
            timeStep = timeStep * SCALE;
            if(dbout) std::cout << "SW groesser: " << (timeStep * SCALE) << "\n";
        }
        //berechnete SW < als min. zulaessige SW in diesem Schritt => urspruengliche SW durch Skalierungsfaktor dividieren
        else if(calculatedTimeStep / timeStep < 1.0/SCALE){
            timeStep = timeStep / SCALE;
             if(dbout) std::cout << "SW kleiner: " << (timeStep / SCALE) << "\n";
            
        }
        //berechnete SW ist OK => annehmen
        else{
            timeStep = calculatedTimeStep;
            if(dbout) std::cout << "SW akzeptiert: " << timeStep << "\n";
        }
        
        //max. zulaessige ges. SW beschraenken
        if (fabs(timeStep) > timeStepMAX) {
            timeStep = timeStepMAX * timeStep /fabs(timeStep);
           if(dbout)  std::cout << "SW beschraenkt(MAX): " << timeStep << "\n";
        }
        
        //min. zulaessige ges. SW beschraenken
        if (fabs(timeStep) < timeStepMIN) {
            timeStep = timeStepMIN * timeStep /fabs(timeStep);
        if(dbout) std::cout << "SW beschraenkt(MIN): " << timeStep << "\n";
        }
        
        if(dbout) std::cout << "Verwendete SW: " << timeStep << "\n";
        
        // *** lokalen Diskretisierungsfehler testen ***
        
        //Diskretisierungsfehler klein genug
        if (fabs(trunc_error - EPS) < 1e-8)  {
            
            //naechsten Zeitpunkt berechnen
            tNext = tCurrent + timeStep;
            
            //Zeitpunkt und Werte fuer naechste Iteration weitersetzen
            tCurrent = tNext;
            xCurrent = xNextFull;
            
            if(dbout) std::cout << "Diskr. Fehler OK: " << trunc_error << "\n";
            
            //Output
            //Consolen-Ausgabe
            if (consoleOutput) {
                std::cout << "T(" << iterationCountAcc << "): " << tCurrent << ", X(" << iterationCountAcc << "): ";
                
                for (size_t i=0; i<state_vec_size; i++) {
                    std::cout << xCurrent[i] << ", ";
                }
                
                std::cout << "\n";
                
                if(dbout){
                std::cout << "T(" << iterationCountAcc << "): " << tCurrent << ", X2Half(" << iterationCountAcc << "): ";
                
                for (size_t i=0; i<state_vec_size; i++) {
                    std::cout << xNextHalf[i] << ", ";
                }
                
                std::cout << "\n";
                }
                
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
            
        }
        
        //Diskretisierungsfehler zu gross => Schritt wiederholen
        else{
            iterationCountAcc--;
            //naechster Zeitpunkt = aktueller Zeitpunkt (nochmal)
            tNext = tCurrent;
            tCurrent = tNext;
            xCurrent = xCurrent;
            if(dbout) std::cout << "Diskr. Fehler zu gross: " << trunc_error << "\n";
        }
        
        if(dbout) std::cout << "\n";
        
        
    } while (tNext < endTime-timeStep);
    x = xNextFull;
    
    return 0;
    
}

#endif
