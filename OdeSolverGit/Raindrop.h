//
//  Raindrop.h
//  OdeSolver
//
//  Die Klasse stellt einen Regentropfen, bzw. ein kreisfoermiges Gebiet dar, das beim angegebenen Radius eine Unstetigkeit aufweist
//  Es existieren 2 verschiedene Funktionen, um den Schnittpunkt eines Lichtstrahls mit dem Kreisrand zu bestimmen:
//  1. Interpolation der beiden Punkte um die Unstetigkeit durch einen kubischen Hermite-Spline und anschliessende Loesung mit dem Newton-Verfahren
//  2. Rekursives Vorgehen mit Halbierung des Schrittweite bei Ueberschreiten der Unstetigkeit


#ifndef __OdeSolver__Raindrop__
#define __OdeSolver__Raindrop__

#include <iostream>
#include "integrate_const.h"
#include "Definitions.h"
#include "IntersectionFunctors.h"

class Raindrop{
    
public:
    Raindrop(){
        mRadius = 1.0;
    }
    template <class ScalarType>
    Raindrop(ScalarType radius){
        mRadius = radius;
    }
    Raindrop(const Raindrop& otherRaindrop ){
        mRadius = otherRaindrop.mRadius;
    }
    ~Raindrop(){
        
    }
    
    double mRadius;
    size_t iterationCountAdaptive;
    
    //Funktion zur Interpolation der beiden Punkte um die Unstetigkeit durch einen kubischen Hermite-Spline und anschliessende Loesung mit dem Newton-Verfahren
    // return: Iterationen des Newton-Verfahrens
    template<class Stepper, class System, class State, class Time>
    size_t calcBoundaryIntersectionHermite(Stepper stepper, System sys, State &xCurrent, Time timeStep, bool consoleOutput = false){
        
        
        // *** 1. State des letzten Punktes innerhalb des Kreises und des ersten Punktes ausserhalb des Kreises bestimmen ***
        
        // xCurrent == [px, py, x, y]
        double &px_ = xCurrent[0]; double &py_ = xCurrent[1]; double &x_ = xCurrent[2]; double &y_ = xCurrent[3];
        
        size_t state_vec_size = xCurrent.size();
        size_t iterationCount = 0;
        
        //F wird in-place aktualisiert (Abstand zum Kreisrand => wenn neg. dann innerhalb des Kreises)
        double F =  x_*x_+y_*y_ - mRadius*mRadius;
        
        //Start bei t = 0
        Time tCurrent = 0.0;
        //Werte der vorigen Iteration
        Time tPrevious = tCurrent;
        State xPrevious = xCurrent; //initialisieren
        
        //solange der aktuelle Punkt im Kreis ist, wird ein Iterationsschritt ausgefuehrt
        while (F < 0) {
            iterationCount++;
            
            //Werte der vorigen Iteration speichern
            xPrevious = xCurrent;
            tPrevious = tCurrent;
            //Dateiausgabe von integrate_const verhindern
            integrate_const(stepper, sys, xCurrent, tCurrent, timeStep, timeStep, false, false);
            tCurrent = tCurrent+timeStep;
            F =  x_*x_+y_*y_ - mRadius*mRadius;
            
        }
        
        //States abspeichern
        State xInCircle = xPrevious;
        Time  tInCircle = tPrevious;
        State xOutOfCircle = xCurrent;
        Time  tOutOfCircle = tCurrent;
        
        //File-Output IntersectionPoint
        const int PRECISION = 7;
        const int WIDTH     = 20;
        std::stringstream filenameIP;
        filenameIP << "Solution" << stepper.stepperName << "IntersectionPoint.dat";
        
        //testen, ob Datei schon existiert
        int fileCountIP = 1;
        while (file_exists(filenameIP.str())) {
            filenameIP.str("");
            filenameIP << "Solution" << stepper.stepperName << "IntersectionPoint(" << fileCountIP << ").dat";
            fileCountIP++;
        }
        std::ofstream write_fileIP(filenameIP.str());
        write_fileIP.precision(PRECISION);
        write_fileIP.setf(std::ios::scientific | std::ios::showpoint);
        
        
        // *** 2. Funktion F(s) definieren ***

        // die beiden Punkte und den Radius benutzen, um die Funktion F(s) zu definieren
        intersectionEquation<RallNo<double>, RallNo<double>>  Fs(xInCircle, xOutOfCircle, mRadius);
        
        // *** 3. Newton-Verfahren mit AD anwenden ***
        
        // Startpunkt des Newton-Verfahrens waehlen
        state_type_double startvec(1);
        startvec[0] = 0.5;
        // Newton-Verfahren zur Loesung der nichtlinearen Gleiung F(s) = x(s)^2+y(s)^2-r^2 anwenden
//        size_t iterationsNewton = newton_raphson(Fs, Fsds, startvec,false, 1000, 1e-6);
        
        size_t iterationsNewton = newton_raphson_ad(Fs, startvec, false, 1000, 1e-6);
        
        
        if (consoleOutput) {
            std::cout << "SP bei Bahnparameter s = " << startvec[0] << "\n";
        }
        

        
        
        // *** 4. Bestimmung des Schnittpunktes durch Einsetzen in die Spline-Gleichung ***
        
        //x(s), y(s)
        hermite_systemXY xy_s(xInCircle, xOutOfCircle);
        state_type_double intersectionPoint(2);
        //berechneten Bahnparamter des Schnittpunktes in die Spline-Gleichung einsetzen
        xy_s(startvec, intersectionPoint);
        
        if (consoleOutput) {
            std::cout << "xSP:" << intersectionPoint[0] << "\n";
            std::cout << "ySP:" << intersectionPoint[1] << "\n";
        }

        
        //File-Output
        write_fileIP << startvec[0] << std::setw(WIDTH);
        write_fileIP << intersectionPoint[0] << std::setw(WIDTH) << intersectionPoint[1] << std::setw(WIDTH);
        write_fileIP << "\n";
        write_fileIP.close();
        //end Output
        
        
        // *** 5. Hermite-Spline Punkte berechen fuer Grafik ***
        
        //File-Output
        std::stringstream filenameHS;
        filenameHS << "Solution" << stepper.stepperName << "Spline.dat";
        //testen, ob Datei schon existiert
        int fileCountHS = 1;
        while (file_exists(filenameHS.str())) {
            filenameHS.str("");
            filenameHS << "Solution" << stepper.stepperName << "Spline(" << fileCountHS << ").dat";
            fileCountHS++;
        }
        std::ofstream write_fileHS(filenameHS.str());
        write_fileHS.precision(PRECISION);
        write_fileHS.setf(std::ios::scientific | std::ios::showpoint);
        

        // x == [px, py, x, y]
        double px_In = xInCircle[0]; double py_In = xInCircle[1]; double x_In = xInCircle[2]; double y_In = xInCircle[3];
        double px_Out = xOutOfCircle[0]; double py_Out = xOutOfCircle[1]; double x_Out = xOutOfCircle[2]; double y_Out = xOutOfCircle[3];
        
        const int SEGMENTS = 5;
        
        //Intervall auf s = [0,1] normieren
        for (int i = 0; i<SEGMENTS+1; i++) {
            double s = (double)i/SEGMENTS;

            //blending-Functions der Hermite Kurve berechnen
            // b(s)
            double b0 = 2*s*s*s - 3*s*s + 1;
            double b1 = -2*s*s*s + 3*s*s;
            double b2 = s*s*s - 2*s*s + s;
            double b3 = s*s*s - s*s;
            
            //Koordinaten x(s) und y(s) der Kurve bei s berechnen
            double xHermite = b0*x_In + b1*x_Out + b2*px_In + b3*px_Out;
            double yHermite = b0*y_In + b1*y_Out + b2*py_In + b3*py_Out;
            
            //Console-Output
            if (consoleOutput) {
                std::cout << "s(" << i << "):" << s << "\n";
                std::cout << "x(" << i << "):" << xHermite << "\n";
                std::cout << "y(" << i << "):" << yHermite << "\n";
            }
            
            //File-Output
            write_fileHS << s << std::setw(WIDTH);
            write_fileHS << xHermite << std::setw(WIDTH) << yHermite << std::setw(WIDTH);
            write_fileHS << "\n";
            
            //end Output

        }
        write_fileHS.close();
        
        return iterationsNewton;
    
    }
    

    //Funktion zur rekursiven Berechnung des Schnittpunktes
    //Funktion initialisiert die Startwerte und die Datei-, Consolenausgabe
    //ruft die Funktion calcBoundaryIntersectionRek rekursiv auf mit xCurrent aus Vorgabe und tCurrent = 0 als Startwerte
    template<class Stepper, class System, class State, class Time>
    size_t calcBoundaryIntersectionAdaptive(Stepper stepper, System sys, State &xCurrent, Time timeStep, bool consoleOutput = false){
        
        size_t state_vec_size = xCurrent.size();
        iterationCountAdaptive = 0;
        
        //Start bei t = 0
        Time tCurrent = 0.0;
        
        //Console-Output Startwerte
        if (consoleOutput) {
            std::cout << "T(" << iterationCountAdaptive << "): " << tCurrent << ", X(" << iterationCountAdaptive << "): ";
            for (size_t i=0; i<state_vec_size; i++) {
                std::cout << xCurrent[i] << ", ";
            }
            std::cout << "\n";
        }
        
        //File-Output initialisieren
        const int PRECISION = 7;
        const int WIDTH     = 20;
        std::stringstream filename;
        filename << "Solution" << stepper.stepperName << ".dat";
        
        //testen, ob Datei schon existiert
        int fileCount = 1;
        while (file_exists(filename.str())) {
            filename.str("");
            filename << "Solution" << stepper.stepperName << "(" << fileCount << ").dat";
            fileCount++;
        }
        
        std::ofstream write_file(filename.str());
        write_file.precision(PRECISION);
        write_file.setf(std::ios::scientific | std::ios::showpoint);
        
        //Startwerte schreiben
        write_file << tCurrent << std::setw(WIDTH);
        for (size_t i=0; i<state_vec_size; i++) {
            write_file << xCurrent[i] << std::setw(WIDTH);
        }
        write_file << "\n";
        // end Output
        
        //die rekursive Funktion wird aufgerufen (hier nicht rekursiv)
        calcBoundaryIntersectionAdaptiveRek(stepper, sys, xCurrent, tCurrent, timeStep, iterationCountAdaptive, write_file, WIDTH, consoleOutput);
        
        return iterationCountAdaptive;
        
    }
    
    //rekursive Funktion zur naeherungsweisen Berechnung des Schnittpunktes mit dem Rand
    //ruft sich rekursiv, ausgehend vom jeweils letzten gueltigen Punkt, mit halber Schrittweite auf
    template<class Stepper, class System, class State, class Time>
    void calcBoundaryIntersectionAdaptiveRek(Stepper stepper, System sys, State &xCurrent, Time startTime, Time timeStep, size_t iterationCount, std::ofstream &write_file, const int WIDTH,  bool consoleOutput = false){
        
        size_t state_vec_size = xCurrent.size();
        // xCurrent == [px, py, x, y]
        double &px_ = xCurrent[0]; double &py_ = xCurrent[1]; double &x_ = xCurrent[2]; double &y_ = xCurrent[3];
        //F wird in-place aktualisiert (Abstand zum Kreisrand => wenn neg. dann innerhalb des Kreises)
        double F =  x_*x_+y_*y_ - mRadius*mRadius;
        
        Time tCurrent = startTime;
        Time tPrevious = tCurrent; //initialisieren
        State xPrevious = xCurrent;
        
        //solange der aktuelle Punkt im Kreis ist, wird noch ein Iterationsschritt ausgefuehrt
        while (F < 0) {
            iterationCount++;
            
            //Werte der vorigen Iteration speichern
            xPrevious = xCurrent;
            tPrevious = tCurrent;
            //Dateiausgabe von integrate_const verhindern
            integrate_const(stepper, sys, xCurrent, tCurrent, timeStep, timeStep, false, false);
            tCurrent = tCurrent+timeStep;
            F =  x_*x_+y_*y_ - mRadius*mRadius;
            
            //Wenn Wert ausserhalb des Kreises, nicht ausgeben
            if (F < 0) {
                //Consolen-Ausgabe
                if (consoleOutput) {
                    std::cout << "F(" << iterationCount << "): " << F << "\n";
                    
                    std::cout << "T(" << iterationCount << "): " << tCurrent << ", X(" << iterationCount << "): ";
                    for (size_t i=0; i<state_vec_size; i++) {
                        std::cout << xCurrent[i] << ", ";
                    }
                    std::cout << "\n";

                }
                //File-Output
                write_file << tCurrent << std::setw(WIDTH);
                for (size_t i=0; i<state_vec_size; i++) {
                    write_file << xCurrent[i] << std::setw(WIDTH);
                }
                write_file << "\n";
                //end Output
            } else {
                //Iterationsschritt nicht zaehlen
                iterationCount--;
            }
            
        }
        
        //Vom vorletzten Punkt rekursiv mit halber Schrittweite starten
        if (F > 10e-6 ) {
            //die rekursive Funktion wird rekursiv aufgerufen
            calcBoundaryIntersectionAdaptiveRek(stepper, sys, xPrevious, tPrevious, 0.5*timeStep, iterationCount, write_file, WIDTH, consoleOutput);
        }
        
    }
    
};


#endif /* defined(__OdeSolver__Raindrop__) */
