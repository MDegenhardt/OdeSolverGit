//
//  GaussFieldFunctor.h
//  OdeSolver
//
//Gauss-Feld Functor


#ifndef OdeSolver_GaussFieldFunctor_h
#define OdeSolver_GaussFieldFunctor_h

#include "Definitions.h"
#include <cmath>

//x == [Fx, Fy, px, py]
//x == [0.5 grad(n(x,y)), p(x,y)]
template <class T>
class fourth_order_system_gauss{
public:
    void operator()(const std::vector<T>& x, std::vector<T>& dxdt, const T t ){
        T px_ = x[0]; T py_ = x[1]; T x_ = x[2]; T y_ = x[3];
        dxdt[0] = -1.0* x_ *exp(-1.0*(x_*x_+y_*y_));
        dxdt[1] = -1.0* y_ *exp(-1.0*(x_*x_+y_*y_));
        dxdt[2] = px_;
        dxdt[3] = py_;
    }
};

//x == [Fx, Fy, Fz, px, py, pz]
//x == [0.5 grad(n(x,y,z)), p(x,y,z)]
template <class T>
class fourth_order_system_gauss_3d{
public:
    void operator()(const std::vector<T>& x, std::vector<T>& dxdt, const T t ){
        T px_ = x[0]; T py_ = x[1]; T pz_ = x[2]; T x_ = x[3]; T y_ = x[4]; T z_ = x[5];
        dxdt[0] = -1.0* x_ *exp(-1.0*(x_*x_+y_*y_+z_*z_));
        dxdt[1] = -1.0* y_ *exp(-1.0*(x_*x_+y_*y_+z_*z_));
        dxdt[2] = -1.0* z_ *exp(-1.0*(x_*x_+y_*y_+z_*z_));
        dxdt[3] = px_;
        dxdt[4] = py_;
        dxdt[5] = pz_;
    }
};


//x == [Fx, Fy, Fz, px, py, pz]
//x == [0.5 grad(n(x,y,z)), p(x,y,z)]
template <class T>
class fourth_order_system_gauss_offset_2D{
public:
    //Offset des Gauss-Feldes
    T delta_x;
    T delta_y;
    T delta_z;
    //Verstaerkung von xz
    T delta_a;
    //Verstaerkung in der Breite
    T delta_b;
    
    fourth_order_system_gauss_offset_2D(T deltaX, T deltaY, T deltaZ, T deltaA, T deltaB){
        delta_x = deltaX;
        delta_y = deltaY;
        delta_z = deltaZ;
        delta_a = deltaA;
        delta_b = deltaB;
    }
    
    void operator()(const std::vector<T>& x, std::vector<T>& dxdt, const T t ){
        T px_ = x[0]; T py_ = x[1]; T pz_ = x[2]; T x_ = x[3]; T y_ = x[4]; T z_ = x[5];
        
        T x_off = x_ - delta_x;
        T y_off = y_ - delta_y;
        T z_off = z_ - delta_z;
        //y-Werte werden in 2D nicht beachtet
        T exponent = exp(-0.5/delta_b*(x_off*x_off + z_off*z_off));
        dxdt[0] = -1.0*delta_a*x_off*exponent;
//        dxdt[1] = dxdt[1];
        dxdt[2] = -1.0*delta_a*z_off*exponent;
        dxdt[3] = px_;
//        dxdt[4] = dxdt[4];
        dxdt[5] = pz_;
    }
};

//x == [Fx, Fy, Fz, px, py, pz]
//x == [0.5 grad(n(x,y,z)), p(x,y,z)]
template <class T>
class fourth_order_system_gauss_offset_3D{
public:
    //Offset des Gauss-Feldes
    T delta_x;
    T delta_y;
    T delta_z;
    //Verstaerkung von xyz
    T delta_a;
    //Verstaerkung in der Breite
    T delta_b;
    
    fourth_order_system_gauss_offset_3D(T deltaX, T deltaY, T deltaZ, T deltaA, T deltaB){
        delta_x = deltaX;
        delta_y = deltaY;
        delta_z = deltaZ;
        delta_a = deltaA;
        delta_b = deltaB;
    }
    
    void operator()(const std::vector<T>& x, std::vector<T>& dxdt, const T t ){
        T px_ = x[0]; T py_ = x[1]; T pz_ = x[2]; T x_ = x[3]; T y_ = x[4]; T z_ = x[5];
        
        T x_off = x_ - delta_x;
        T y_off = y_ - delta_y;
        T z_off = z_ - delta_z;
        
        T exponent = exp(-0.5/delta_b*(x_off*x_off + y_off*y_off + z_off*z_off));
        dxdt[0] = -1.0*delta_a*x_off*exponent;
        dxdt[1] = -1.0*delta_a*y_off*exponent;
        dxdt[2] = -1.0*delta_a*z_off*exponent;
        dxdt[3] = px_;
        dxdt[4] = py_;
        dxdt[5] = pz_;
    }
};

//x == [Fx, Fy, Fz, px, py, pz]
//x == [0.5 grad(n(x,y,z)), p(x,y,z)]
template <class T>
class fourth_order_system_sphere_offset_3D{
public:
    //Offset des Gauss-Feldes
    T delta_x;
    T delta_y;
    T delta_z;
    
    //Verstaerkung aller Werte
    T delta_a;
    //Verstaerkung in der Breite
    T delta_b;
    
    fourth_order_system_sphere_offset_3D(T deltaX, T deltaY, T deltaZ, T deltaA, T deltaB){
        delta_x = deltaX;
        delta_y = deltaY;
        delta_z = deltaZ;
        delta_a = deltaA;
        delta_b = deltaB;
    }
    
    void operator()(const std::vector<T>& x, std::vector<T>& dxdt, const T t ){
        T px_ = x[0]; T py_ = x[1]; T pz_ = x[2]; T x_ = x[3]; T y_ = x[4]; T z_ = x[5];
        
        T x_off = x_ - delta_x;
        T y_off = y_ - delta_y;
        T z_off = z_ - delta_z;
        
        T denom = (x_off*x_off + y_off*y_off + z_off*z_off)*(x_off*x_off + y_off*y_off + z_off*z_off);
        dxdt[0] = -1.0*delta_a*5*x_off/denom*5*delta_b;
        dxdt[1] = -1.0*delta_a*5*y_off/denom*5*delta_b;
        dxdt[2] = -1.0*delta_a*5*z_off/denom*5*delta_b;
        dxdt[3] = px_;
        dxdt[4] = py_;
        dxdt[5] = pz_;
    }
};




#endif
