//
//  newton_raphson_ad.h
//  OdeSolver
//
//  Newton-Raphson Verfahren mit automatischer Differentation fuer Systeme von Gleichungen mit einer Variable x
//

#ifndef OdeSolver_newton_raphson_ad_h
#define OdeSolver_newton_raphson_ad_h

#include <fstream>
#include <iomanip>
#include "autoderiv.h"



//bekommt 1 Funktor (Vektor der Funktionen) und Vektor der Startwerte (State == state_type_double)
template<class System, class State>
size_t newton_raphson_ad(System Sysx, State& x, bool consoleOutput = false, size_t maxIterations = 1000, double epsilon = 1e-8){
    
    bool fileOutput = false;
    
    size_t state_vec_size = x.size();
    size_t maxIterationCount = 0;
    
    //Iterationen fuer jede Funktion
    size_t iterationCount[state_vec_size];
    
    bool allCalculated = false;
    //ist die gewuenschte Genauigkeit (epsilon) der einzelnen Funktionen erreicht?
    bool isCalculated[state_vec_size];
    for (int i=0; i<state_vec_size; i++) {
        isCalculated[i] = false;
        iterationCount[i] = 0;
    }
    
    //Aktuelle x-Werte
    State xCurrent = x;
    //neu berechnete x-Werte (state nur initialisieren)
    State xNext    = x;
    //Unterschied zur vorigen Iteration
    State error    = x;
    
    //Vektor zum Speichern der Funktionswerte
    State functionEvaluated(state_vec_size);
    //Vektor zum Speichern der Ableitungswerte
    State derivateEvaluated(state_vec_size);
    
    //Vektor von RallNo zum Speichern des Vektors von Funktionswert und Ableitungswert
    state_type_rall_double fEvaluated(state_vec_size);
    //Vektor von RallNo zum Speichern des Vektors der Auswertungsstelle (x-Wert)
    state_type_rall_double xCurrentRall(state_vec_size);
    
    //Console-Output
    if (consoleOutput) {
        std::cout << "X(0): ";
        for (size_t i=0; i<state_vec_size; i++) {
            std::cout << x[i] << ", ";
        }
        std::cout << "\n";
    }
    
    //File-Output
    const int PRECISION = 15;
    const int WIDTH     = 25;
    
    std::stringstream filename;
    filename << "SolutionNewtonRaphson" << ".dat";
    //testen, ob Datei schon existiert
    int fileCount = 1;
    while (file_exists(filename.str())) {
        filename.str("");
        filename << "SolutionNewtonRaphson(" << fileCount << ").dat";
        fileCount++;
    }
    
    std::string fname = filename.str();
    std::ofstream write_file(fname);
    
    if (!fileOutput) {
        //file wieder loeschen
        write_file.close();
        std::remove ( fname.c_str() );
    }
    
    if (fileOutput) {
        std::ofstream write_file(filename.str());
        write_file.precision(PRECISION);
        write_file.setf(std::ios::scientific | std::ios::showpoint);
        
        //Startwerte schreiben
        for (size_t i=0; i<state_vec_size; i++) {
            write_file << xCurrent[i] << std::setw(WIDTH);
        }
        write_file << "\n";
        // end Output
    }

    
    do{
        maxIterationCount++;
        
        xCurrent = xNext;
        // Vektor der x-Werte in Vektor von RallNos umwandeln und jeweils ableiten
        for (int i=0; i<state_vec_size; i++) {
            xCurrentRall[i] = xCurrent[i];
            xCurrentRall[i].derive();
        }
        //f(x) und f'(x) auswerten an xCurrent, speichern in fEvaluated
        Sysx(xCurrentRall, fEvaluated);
        
        // Funktionswert und Ableitung an der Stelle speichern
        for (int i=0; i<state_vec_size; i++) {
            functionEvaluated[i] = (fEvaluated[i]).val();
            derivateEvaluated[i] = (fEvaluated[i]).der();
        }
        
        //Value-Output
        bool valueOutput = false;
        if (valueOutput) {
            std::cout << "FUN(" << maxIterationCount << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << functionEvaluated[i] << ", ";
            }
            std::cout << "\n";
            
            std::cout << "DER(" << maxIterationCount << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << derivateEvaluated[i] << ", ";
            }
            std::cout << "\n";
        }
        
        for (int i=0; i<state_vec_size; i++) {
            
            if (isCalculated[i] == false) {
                xNext[i]    = xCurrent[i] - functionEvaluated[i]/derivateEvaluated[i];
                
                error[i] = std::abs(xCurrent[i]-xNext[i]);
                (iterationCount[i])++;
            }
            
            if (std::abs(error[i]) < epsilon) {
                isCalculated[i] = true;
            }
        }
        
        //Console-Output
        if (consoleOutput) {
            std::cout << "X(" << maxIterationCount << "): ";
            
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << xNext[i] << ", ";
            }
            std::cout << "\n";
        }
        
        if (fileOutput) {
            //File-Output
            
            for (size_t i=0; i<state_vec_size; i++) {
                write_file << xNext[i] << std::setw(WIDTH);
            }
            write_file << "\n";
            //end Output
        }

        
        //alle Werte berechnet(nur true im Vektor) ?
        int calcCount = 0;
        for (int i=0; i<state_vec_size; i++) {
            if (isCalculated[i] == true) {
                calcCount++;
            }
            if (calcCount == state_vec_size) {
                allCalculated = true;
            }
        }
        
    } while ( !allCalculated && (maxIterationCount < maxIterations));
    
    x= xNext;
    
    //Console-Output
    if (consoleOutput) {
        std::cout << "XFinal: ";
        typename State::const_iterator it;
        for (it=x.begin(); it!=x.end(); it++) {
            std::cout << *it << ", ";
        }
        std::cout << "\n";
        std::cout << "Iterationen: ";
        for (int i=0; i<state_vec_size; i++) {
            std::cout << iterationCount[i] << ", ";
        }
        std::cout << "\n";
        std::cout << "Errors: ";
        for (int i=0; i<state_vec_size; i++) {
            std::cout << error[i] << ", ";
        }
        std::cout << "\n";
    }
    if (fileOutput) {
        write_file.close();
    }
    
    
    
    return maxIterationCount;
    
}


#endif
