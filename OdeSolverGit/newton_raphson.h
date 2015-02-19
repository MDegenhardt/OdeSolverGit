//
//  newton_raphson.h
//  OdeSolver
//
//  Newton-Raphson Verfahren fuer Systeme von Gleichungen mit einer Variable x
//

#ifndef OdeSolver_newton_raphson_h
#define OdeSolver_newton_raphson_h

#include <fstream>
#include <iomanip>



//bekommt 2 Funktoren (Vektor der Funktionen, Vektor der Ableitungen) und Vektor der Startwerte
template<class System, class System2, class State>
size_t newton_raphson(System Sysx, System2 Sysdxdt, State& x, bool consoleOutput = false, size_t maxIterations = 1000, double epsilon = 1e-8){
    
    size_t state_vec_size = x.size(); 
    size_t maxIterationCount = 0;

    //Iterationen fuer jede Funktion
    size_t* iterationCount = new size_t[state_vec_size];
    
    bool allCalculated = false;
    //ist die gewuenschte Genauigkeit (epsilon) der einzelnen Funktionen erreicht?
    bool* isCalculated = new bool[state_vec_size];
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
    
    State functionEvaluated(state_vec_size);
    State derivateEvaluated(state_vec_size);

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

//    std::stringstream filename;
//    filename << "SolutionNewtonRaphson" << ".dat";
//    //testen, ob Datei schon existiert
//    int fileCount = 1;
//    while (file_exists(filename.str())) {
//        filename.str("");
//        filename << "SolutionNewtonRaphson(" << fileCount << ").dat";
//        fileCount++;
//    }
//    
//    std::ofstream write_file(filename.str());
//    write_file.precision(PRECISION);
//    write_file.setf(std::ios::scientific | std::ios::showpoint);
//    
//    //Startwerte schreiben
//    for (size_t i=0; i<state_vec_size; i++) {
//        write_file << xCurrent[i] << std::setw(WIDTH);
//    }
//    write_file << "\n";
//    // end Output
    
    do{
        maxIterationCount++;
        
        xCurrent = xNext;
        
        //f(x) auswerten
        Sysx(xCurrent, functionEvaluated);
        //f'(x) auswerten
        Sysdxdt(xCurrent, derivateEvaluated);
        
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
//        //File-Output
//        
//        for (size_t i=0; i<state_vec_size; i++) {
//            write_file << xNext[i] << std::setw(WIDTH);
//        }
//        write_file << "\n";
//        //end Output
  
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

//    write_file.close();

	delete[] iterationCount;
	delete[] isCalculated;
    
    return maxIterationCount;
    
}



#endif
