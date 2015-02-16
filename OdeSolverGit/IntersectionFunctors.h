//
//  IntersectionFunctors.h
//  OdeSolver
//
//  definiert Functoren, die zur Schnittpunktberechnung mit dem Hermite-Spline benoetigt werden

#ifndef OdeSolver_IntersectionFunctors_h
#define OdeSolver_IntersectionFunctors_h

// Gleichungssystem, das den Hermite-Spline darstellt (als Functor)
// x(s) = ...
// y(s) = ...
// z(s) = ...
class hermite_systemXYZ{
    state_type_double xInArea;
    state_type_double xOutOfArea;
public:
    
    hermite_systemXYZ( state_type_double xI, state_type_double xO ){
        xInArea = xI;
        xOutOfArea = xO;
    }
    
    void operator()(const state_type_double& svec, state_type_double& fs){
        
        double s = svec[0];
        
        // x == [px0, py0, pz0, x0, y0, z0]
        double px_In = xInArea[0]; double py_In = xInArea[1]; double pz_In = xInArea[2]; double x_In = xInArea[3]; double y_In = xInArea[4]; double z_In = xInArea[5];
        double px_Out = xOutOfArea[0]; double py_Out = xOutOfArea[1]; double pz_Out = xOutOfArea[2]; double x_Out = xOutOfArea[3]; double y_Out = xOutOfArea[4]; double z_Out = xOutOfArea[5];
        
        double b0 = 2*s*s*s - 3*s*s + 1;
        double b1 = -2*s*s*s + 3*s*s;
        double b2 = s*s*s - 2*s*s + s;
        double b3 = s*s*s - s*s;
        
        //x(s)
        fs[0] = b0*x_In + b1*x_Out + b2*px_In + b3*px_Out;
        //y(s)
        fs[1] = b0*y_In + b1*y_Out + b2*py_In + b3*py_Out;
        //z(s)
        fs[2] = b0*z_In + b1*z_Out + b2*pz_In + b3*pz_Out;
        
    }
};

// Gleichungssystem, das den Hermite-Spline darstellt (als Functor)
// x(s) = ...
// y(s) = ...
class hermite_systemXY{
    state_type_double xInArea;
    state_type_double xOutOfArea;
public:
    
    hermite_systemXY( state_type_double xI, state_type_double xO ){
        xInArea = xI;
        xOutOfArea = xO;
    }
    
    void operator()(const state_type_double& svec, state_type_double& fs){
        
        double s = svec[0];
        
        // x == [px, py, x, y]
        double px_In = xInArea[0]; double py_In = xInArea[1]; double x_In = xInArea[2]; double y_In = xInArea[3];
        double px_Out = xOutOfArea[0]; double py_Out = xOutOfArea[1]; double x_Out = xOutOfArea[2]; double y_Out = xOutOfArea[3];
        
        double b0 = 2*s*s*s - 3*s*s + 1;
        double b1 = -2*s*s*s + 3*s*s;
        double b2 = s*s*s - 2*s*s + s;
        double b3 = s*s*s - s*s;
        
        //x(s)
        fs[0] = b0*x_In + b1*x_Out + b2*px_In + b3*px_Out;
        //y(s)
        fs[1] = b0*y_In + b1*y_Out + b2*py_In + b3*py_Out;
        
    }
};

//Functor, der die Gleichung F(s) = x(s)^2+y(s)^2-r^2 darstellt
// Konstruktor erwartet State-Type des letzten inneren Punktes im Kreis und des ersten aeusseren Punktes und den Radius des Kreises
template<class T, class U>
class intersectionEquation{
    state_type_double xInCircle;
    state_type_double xOutOfCircle;
    double r;
public:
    
    intersectionEquation( state_type_double xI, state_type_double xO, double radius ){
        xInCircle = xI;
        xOutOfCircle = xO;
        r = radius;
    }
    
    void operator()(const std::vector<T>& svec, std::vector<U>& result){
        
        RallNo<double> s = svec[0];
        
        // x == [px, py, x, y]
        double px_In = xInCircle[0]; double py_In = xInCircle[1]; double x_In = xInCircle[2]; double y_In = xInCircle[3];
        double px_Out = xOutOfCircle[0]; double py_Out = xOutOfCircle[1]; double x_Out = xOutOfCircle[2]; double y_Out = xOutOfCircle[3];
        
        // blending-functions
        // blending-functions
        RallNo<double> b0 = 2*s*s*s - 3*s*s + 1;
        RallNo<double> b1 = -2*s*s*s + 3*s*s;
        RallNo<double> b2 = s*s*s - 2*s*s + s;
        RallNo<double> b3 = s*s*s - s*s;
        
        //F(s) = x(s)^2 + y(s)^2 -r^2
        
        //x(s)
        result[0] = (b0*x_In + b1*x_Out + b2*px_In + b3*px_Out)*
        (b0*x_In + b1*x_Out + b2*px_In + b3*px_Out)+
        //y(s)
        (b0*y_In + b1*y_Out + b2*py_In + b3*py_Out)*
        (b0*y_In + b1*y_Out + b2*py_In + b3*py_Out)-
        r*r;
        
    }
};

////Functor, der die Ableitung von F(s), also F'(s) = 2*x(s)*x'(s)+ 2*y(s)*y'(s) darstellt
//// Konstruktor erwartet State-Type des letzten inneren Punktes im Kreis und des ersten aeusseren Punktes
//class intersectionEquationDerivate{
//    state_type_double xInCircle;
//    state_type_double xOutOfCircle;
//public:
//    
//    intersectionEquationDerivate( state_type_double xI, state_type_double xO ){
//        xInCircle = xI;
//        xOutOfCircle = xO;
//    }
//    
//    void operator()(const state_type_double& svec, state_type_double& dFds){
//        
//        double s = svec[0];
//        
//        // x == [px, py, x, y]
//        double px_In = xInCircle[0]; double py_In = xInCircle[1]; double x_In = xInCircle[2]; double y_In = xInCircle[3];
//        double px_Out = xOutOfCircle[0]; double py_Out = xOutOfCircle[1]; double x_Out = xOutOfCircle[2]; double y_Out = xOutOfCircle[3];
//        
//        double b0 = 2*s*s*s - 3*s*s + 1;
//        double b1 = -2*s*s*s + 3*s*s;
//        double b2 = s*s*s - 2*s*s + s;
//        double b3 = s*s*s - s*s;
//        double b0ds = 6*s*s-6*s;
//        double b1ds = -6*s*s+6*s;
//        double b2ds = 3*s*s-4*s+1;
//        double b3ds = 3*s*s-2*s;
//        
//        //x(s)
//        dFds[0] = 2*(b0*x_In + b1*x_Out + b2*px_In + b3*px_Out)*
//        //dxds(s)
//        (b0ds*x_In + b1ds*x_Out + b2ds*px_In + b3ds*px_Out)+
//        //y(s)
//        2*(b0*y_In + b1*y_Out + b2*py_In + b3*py_Out)*
//        //dyds(s)
//        (b0ds*y_In + b1ds*y_Out + b2ds*py_In + b3ds*py_Out);
//        
//    }
//};


//Functor, der die Gleichung F(s) = z(s)+z_ncp darstellt
// Konstruktor erwartet State-Type des letzten inneren Punktes im Sichtfeld und des ersten aeusseren Punktes und den Abstand zur NCP
template<class T, class U>
class intersectionEquationPlane3d{
    state_type_double xInPlane;
    state_type_double xOutOfPlane;
    double z_ncp;
public:
    
    intersectionEquationPlane3d( state_type_double xI, state_type_double xO, double zncp_ ){
        xInPlane = xI;
        xOutOfPlane = xO;
        z_ncp = zncp_;
    }
    
    void operator()(const std::vector<T>& svec, std::vector<U>& result){
        
        RallNo<double> s = svec[0];
        
        // x == [px0, py0, pz0, x0, y0, z0]
        double pz_In = xInPlane[2]; double z_In = xInPlane[5];
        double pz_Out = xOutOfPlane[2]; double z_Out = xOutOfPlane[5];
        
        // blending-functions
        RallNo<double> b0 = 2*s*s*s - 3*s*s + 1;
        RallNo<double> b1 = -2*s*s*s + 3*s*s;
        RallNo<double> b2 = s*s*s - 2*s*s + s;
        RallNo<double> b3 = s*s*s - s*s;
        
        //F(s) = z(s) + z_ncp
        result[0] = (b0*z_In + b1*z_Out + b2*pz_In + b3*pz_Out)+ z_ncp;
        
    }
};

//Functor, der die Gleichung F(s) = z(s)+z_ncp darstellt
// Konstruktor erwartet State-Type des letzten inneren Punktes im Sichtfeld und des ersten aeusseren Punktes und den Abstand zur NCP
template<class T, class U>
class intersectionEquationPlane2d{
    state_type_double xInPlane;
    state_type_double xOutOfPlane;
    double z_ncp;
public:
    
    intersectionEquationPlane2d( state_type_double xI, state_type_double xO, double zncp_ ){
        xInPlane = xI;
        xOutOfPlane = xO;
        z_ncp = zncp_;
    }
    
    void operator()(const std::vector<T>& svec, std::vector<U>& result){
        
        RallNo<double> s = svec[0];
        
        // x == [pz, px, z, x]
        double pz_In = xInPlane[0]; double px_In = xInPlane[1]; double z_In = xInPlane[2]; double x_In = xInPlane[3];
        double pz_Out = xOutOfPlane[0]; double px_Out = xOutOfPlane[1]; double z_Out = xOutOfPlane[2]; double x_Out = xOutOfPlane[3];
        
        // blending-functions
        RallNo<double> b0 = 2*s*s*s - 3*s*s + 1;
        RallNo<double> b1 = -2*s*s*s + 3*s*s;
        RallNo<double> b2 = s*s*s - 2*s*s + s;
        RallNo<double> b3 = s*s*s - s*s;
        
        //F(s) = z(s) + z_ncp
        result[0] = (b0*z_In + b1*z_Out + b2*pz_In + b3*pz_Out)+ z_ncp;
        
    }
};





#endif
