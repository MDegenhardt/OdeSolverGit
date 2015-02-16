//
//  Definitions.h
//  OdeSolver
//
//  Einige Definitionen, die in anderen Dateien benoetigt werden


#ifndef OdeSolver_Definitions_h
#define OdeSolver_Definitions_h

#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>
#include "autoderiv.h"


typedef std::vector<double> state_type_double;
typedef std::vector<RallNo<double>> state_type_rall_double;
typedef std::vector<std::complex<double>> state_type_complex;
typedef std::complex<double> cdouble;
typedef std::valarray<double> val_array;

//testen ob eine Datei schon existiert
bool file_exists(const std::string& name)
{
    std::ifstream file(name);
    if(!file)    //Datei nicht gefunden -> file ist 0
        return false;    //Datei nicht gefunden
    else         //Datei gefunden -> file ist ungleich 0
        return true;     //Datei gefunden
}


#endif
