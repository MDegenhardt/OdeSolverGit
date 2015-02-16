//
//  PicturePlane.h
//  OdeSolver
//
// fuehrt Zeitintegration eines Strahls von der FCP bis zur NCP durch mittels Hermite-Spline

#ifndef OdeSolver_PicturePlane_h
#define OdeSolver_PicturePlane_h

#include <iostream>
#include "integrate_const.h"
#include "Definitions.h"
#include "IntersectionFunctors.h"


class PicturePlane{
    
public:
    
    PicturePlane(double z_ncp_, double z_fcp_){
        z_ncp = z_ncp_;
        z_fcp = z_fcp_;
    }
    
    PicturePlane(const PicturePlane& otherPlane){
        z_ncp = otherPlane.z_ncp;
        z_fcp = otherPlane.z_fcp;
    }
    
    ~PicturePlane(){
        
    }
    
    // z-Koordinate der Near Clipping Plane
    double z_ncp;
    //    // z-Koordinate der Far Clipping Plane
    double z_fcp;
    
    //Funktion zur Interpolation der beiden Punkte um die near clipping plane (NCP) durch einen kubischen Hermite-Spline und anschliessende Loesung mit dem Newton-Verfahren
    // return: Iterationen des Newton-Verfahrens
    template<class Stepper, class System, class State, class Time>
    size_t calcNCPIntersectionHermite3d(Stepper stepper, System sys, State& xCurrent, Time timeStep, State& intersectionPoint, bool consoleOutput = false){
        
        //        consoleOutput = true;
        
        // *** 1. State des letzten Punktes innerhalb des Sichtvolumens (zwischen FCP und NCP) und des ersten Punktes ausserhalb des Sichtvolumens bestimmen ***
        
        // xCurrent == [px0, py0, pz0, x0, y0, z0]
        double &z_ = xCurrent[5];
        
        size_t state_vec_size = xCurrent.size();
        size_t iterationCount = 0;
        
        //F wird in-place aktualisiert (Abstand zur NCP => wenn neg. dann innerhalb des Sichtvolumens)
        double F =  z_+z_ncp;
        
        //Start bei t = 0
        Time tCurrent = 0.0;
        //Werte der vorigen Iteration
        Time tPrevious = tCurrent;
        State xPrevious = xCurrent; //initialisieren
        
        //solange der aktuelle Punkt im Sichtvolumen ist, wird ein Iterationsschritt ausgefuehrt
        // nach n Iterationen abbrechen, da der Strahl die NCP wahrscheinlich nicht erreicht -> Punkt x=0, y=0, z=z_ncp zurueckgeben
        int n = 1000;
        bool rayReturned = false;
        while (F < 0 ) {
            iterationCount++;
            if (iterationCount > n) {
                std::cout << "Ray returned \n";
                rayReturned = true;
                break;
            }
            
            //Werte der vorigen Iteration speichern
            xPrevious = xCurrent;
            tPrevious = tCurrent;
            //Dateiausgabe von integrate_const verhindern
            integrate_const(stepper, sys, xCurrent, tCurrent, timeStep, timeStep, false, false);
            tCurrent = tCurrent+timeStep;
            F =  z_+z_ncp;
            
        }
        
        
        //States abspeichern
        State xInArea = xPrevious;
        Time  tInArea = tPrevious;
        State xOutOfArea = xCurrent;
        Time  tOutOfArea = tCurrent;
        
        //        //File-Output IntersectionPoint
        const int PRECISION = 7;
        const int WIDTH     = 20;
        size_t iterationsNewton = 0;
        
        //wenn der Strahl die NCP erreicht hat
        if (!rayReturned) {
            
            // *** 2. Funktion F(s) definieren ***
            // die beiden Punkte und den Radius benutzen, um die Funktion F(s) zu definieren
            intersectionEquationPlane3d<RallNo<double>, RallNo<double>>  Fs(xInArea, xOutOfArea, z_ncp);
            
            // *** 3. Newton-Verfahren mit AD anwenden ***
            // Startpunkt des Newton-Verfahrens waehlen
            State startvec(1);
            startvec[0] = 0.5;
            // Newton-Verfahren zur Loesung der nichtlinearen Gleiung F(s) = z(s)+z_ncp anwenden
            iterationsNewton = newton_raphson_ad(Fs, startvec, false, 1000, 1e-6);
            
            if (consoleOutput) {
                std::cout << std::scientific << "SP bei Bahnparameter s = " << startvec[0] << "\n";
            }
            
            // *** 4. Bestimmung des Schnittpunktes durch Einsetzen in die Spline-Gleichung ***
            // x(s), y(s), z(s)
            hermite_systemXYZ xyz_s(xInArea, xOutOfArea);
            //berechneten Bahnparameter des Schnittpunktes in die Spline-Gleichung einsetzen
            xyz_s(startvec, intersectionPoint);
            
        }
        // Strahl hat NCP nicht erreicht -> Pixel in den Ursprung setzen
        else {
            intersectionPoint[0] = 0.0;
            intersectionPoint[1] = 0.0;
            intersectionPoint[2] = z_ncp;
        }
        if (consoleOutput) {
            std::cout << "xSP:" << intersectionPoint[0] << "\n";
            std::cout << "ySP:" << intersectionPoint[1] << "\n";
            std::cout << "zSP:" << intersectionPoint[2] << "\n";
        }
        
        return iterationsNewton;
        
    }
    
    
    
};


#endif
