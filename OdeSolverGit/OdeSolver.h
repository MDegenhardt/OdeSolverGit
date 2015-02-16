//
//  OdeSolver.h
//  OdeSolver
//
//  Bindet alle sonstigen header-files ein, die fuer den OdeSolver benoetigt werden

#ifndef OdeSolver_OdeSolver_h
#define OdeSolver_OdeSolver_h

#include <iostream>
#include <cmath>
#include <complex>
#include <valarray>


#include "Definitions.h"
#include "integrate_const.h"
#include "integrate_adaptive.h"
#include "newton_raphson.h"
#include "newton_raphson_ad.h"
#include "autoderiv.h"
#include "ExampleFunctors.h"
#include "NonlinearFunctors.h"
#include "GaussFieldFunctor.h"
#include "Raindrop.h"
#include "runExamples.h"
#include "PicturePlane.h"
#include "GaussPicture.h"
#include "BMPLoader.h"

#include "Euler.h"
#include "RungeKutta4.h"
#include "Dopri5.h"
#include "Fehlberg78.h"
#include "ButcherTable.h"
#include "RungeKutta14.h"
#include "implicit_euler.h"
#include "implicit_euler_ad.h"

#endif
