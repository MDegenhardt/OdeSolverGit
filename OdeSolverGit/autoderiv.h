//
//  autoderiv.h
//  AutoDeriv
//


#ifndef AutoDeriv_autoderiv_h
#define AutoDeriv_autoderiv_h

#include <cmath>

template <typename T>
class RallNo{
    
private:
    T m_val;
    T m_der;
    bool m_containsX;
    
public:
    
    //    Konstruktor
    RallNo(){
        m_containsX = false;
    }
    //    Copy-Konstuktor
    RallNo(const RallNo<T> &otherRall){
        m_val = otherRall.m_val;
        m_containsX = otherRall.m_containsX;
        //        wenn die Ableitung nicht konstant -> kopieren
        if (m_containsX) {
            m_der = otherRall.m_der;
        }
    }
    //    Konstruktor mit Value Eingabe (konst. Zahl -> RallNo)
    template <class ConstType>
    RallNo(const ConstType &cVal){
        m_val = cVal;
        m_containsX = false;
    }
    
    //    = Operator macht aus Zahl eine RallNo
    template <class ConstType>
    RallNo<T>& operator=(const ConstType &cVal){
        m_val = cVal;
        m_containsX = false;
        return *this;
    }
    
    //    = Operator weist Werte einer RallNo einer anderen RallNo zu
    RallNo<T>& operator=(const RallNo<T> &otherRallNo){
        if (this==&otherRallNo) {
            return *this;
        }
        m_val = otherRallNo.m_val;
        m_containsX = otherRallNo.m_containsX;
        if (m_containsX) {
            m_der = otherRallNo.m_der;
        }
        return *this;
    }
    //    gibt Wert der Funktion zurueck
    const T& val() const{
        return m_val;
    }
    
    T& val() {
        return m_val;
    }
    //    gibt Wert der Ableitung zurueck
    const T& der() const{
        return m_der;
    }
    //    setzt den Wert der Ableitung
    void setDer(T der){
        m_der = der;
    }
    //    Ableitung = 0 wenn Fkt = const
    T& der()
    {
        if (m_containsX) return m_der;
        static T z = 0.0;
        return z;
    }
    
    //    Ableitung von x -> d/dx = 1
    T& derive(){
        m_der = 1.0;
        m_containsX = true;
        return m_der;
    }
    
    //    ist eine RallNo von x abhaengig?
    bool containsX() const{
        return m_containsX;
    }
    
    //    Abhaengigkeit einer RallNo von X setzen
    void setContainsX(){
        m_containsX = true;
    }
    
};

// *** ADDITION ***

// a=const, b=f(x)
template<typename T, typename ConstType>
RallNo<T> add1(const ConstType &a, const RallNo<T> &b){
    
    RallNo<T> c(a+b.val());
    if (!b.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(b.der());
    return c;
}

// a=f(x), b=const
template<typename T, typename ConstType>
RallNo<T> add2(const RallNo<T> &a, const ConstType &b){
    
    RallNo<T> c(a.val()+b);
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(a.der());
    return c;
}



template<typename T, typename ConstType>
RallNo<T> operator+(const ConstType &a, const RallNo<T> &b){
    return add1(a,b);
}

template<typename T, typename ConstType>
RallNo<T> operator+(const RallNo<T> &a, const ConstType &b){
    return add2(a,b);
}


template<typename T>
RallNo<T> add3(const RallNo<T> &a, const RallNo<T> &b){
    RallNo<T> c(a.val()+b.val());
    //c ist abh. von x
    c.setContainsX();
    //Ableitungen addieren
    c.setDer(a.der()+b.der());
    return c;
}


template<typename T>
RallNo<T> operator+(const RallNo<T> &a, const RallNo<T> &b){
    switch ((a.containsX()?1:0)|(b.containsX()?2:0)) {
            // a!=f(x), b!=f(x)
        case 0:
            return RallNo<T>(a.val()+b.val());
            // a=f(x), b!=f(x)
        case 1:
            return add2(a,b.val());
            // a!=f(x), b=f(x)
        case 2:
            return add1(a.val(),b);
    }
    // a=f(x), b=f(x)
    return add3(a,b);
}

// *** SUBTRAKTION ***

// a=const, b=f(x)
template<typename T, typename ConstType>
RallNo<T> sub1(const ConstType &a, const RallNo<T> &b){
    
    RallNo<T> c(a-b.val());
    if (!b.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(-b.der());
    return c;
}

// a=f(x), b=const
template<typename T, typename ConstType>
RallNo<T> sub2(const RallNo<T> &a, const ConstType &b){
    
    RallNo<T> c(a.val()-b);
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(a.der());
    return c;
}



template<typename T, typename ConstType>
RallNo<T> operator-(const ConstType &a, const RallNo<T> &b){
    return sub1(a,b);
}

template<typename T, typename ConstType>
RallNo<T> operator-(const RallNo<T> &a, const ConstType &b){
    return sub2(a,b);
}


template<typename T>
RallNo<T> sub3(const RallNo<T> &a, const RallNo<T> &b){
    RallNo<T> c(a.val()-b.val());
    //c ist abh. von x
    c.setContainsX();
    //Ableitungen subtrahieren
    c.setDer(a.der()-b.der());
    return c;
}


template<typename T>
RallNo<T> operator-(const RallNo<T> &a, const RallNo<T> &b){
    switch ((a.containsX()?1:0)|(b.containsX()?2:0)) {
            // a!=f(x), b!=f(x)
        case 0:
            return RallNo<T>(a.val()-b.val());
            // a=f(x), b!=f(x)
        case 1:
            return sub2(a,b.val());
            // a!=f(x), b=f(x)
        case 2:
            return sub1(a.val(),b);
    }
    // a=f(x), b=f(x)
    return sub3(a,b);
}


// *** MULTIPLIKATION ***

// a=const, b=f(x)
template<typename T, typename ConstType>
RallNo<T> mul1(const ConstType &a, const RallNo<T> &b){
    
    RallNo<T> c(a*b.val());
    if (!b.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(b.der()*a);
    return c;
}

// a=f(x), b=const
template<typename T, typename ConstType>
RallNo<T> mul2(const RallNo<T> &a, const ConstType &b){
    
    RallNo<T> c(a.val()*b);
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(a.der()*b);
    return c;
}

template<typename T, typename ConstType>
RallNo<T> operator*(const ConstType &a, const RallNo<T> &b){
    return mul1(a,b);
}

template<typename T, typename ConstType>
RallNo<T> operator*(const RallNo<T> &a, const ConstType &b){
    return mul2(a,b);
}


template<typename T>
RallNo<T> mul3(const RallNo<T> &a, const RallNo<T> &b){
    RallNo<T> c(a.val()*b.val());
    //c ist abh. von x
    c.setContainsX();
    //Ableitung nach Kettenregel berechnen
    c.setDer(a.der()*b.val()+b.der()*a.val());
    return c;
}

template<typename T>
RallNo<T> operator*(const RallNo<T> &a, const RallNo<T> &b){
    switch ((a.containsX()?1:0)|(b.containsX()?2:0)) {
            // a!=f(x), b!=f(x)
        case 0:
            return RallNo<T>(a.val()*b.val());
            // a=f(x), b!=f(x)
        case 1:
            return mul2(a,b.val());
            // a!=f(x), b=f(x)
        case 2:
            return mul1(a.val(),b);
    }
    // a=f(x), b=f(x)
    return mul3(a,b);
}


// *** DIVISION ***

// a=const, b=f(x)
template<typename T, typename ConstType>
RallNo<T> div1(const ConstType &a, const RallNo<T> &b){
    
    RallNo<T> c(a/b.val());
    if (!b.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c'= -a/b^2*b'
    c.setDer(-a/( b.val()*b.val() ) *b.der());
    return c;
}

// a=f(x), b=const
template<typename T, typename ConstType>
RallNo<T> div2(const RallNo<T> &a, const ConstType &b){
    
    RallNo<T> c(a.val()/b);
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = a'/b
    c.setDer(a.der()/b);
    return c;
}

template<typename T, typename ConstType>
RallNo<T> operator/(const ConstType &a, const RallNo<T> &b){
    return div1(a,b);
}

template<typename T, typename ConstType>
RallNo<T> operator/(const RallNo<T> &a, const ConstType &b){
    return div2(a,b);
}


template<typename T>
RallNo<T> div3(const RallNo<T> &a, const RallNo<T> &b){
    RallNo<T> c(a.val()/b.val());
    //c ist abh. von x
    c.setContainsX();
    //Ableitung nach Quotientenregel berechnen
    //c=a/b, c' = (a'b-ab')/b^2 = (a'-cb')/b
    c.setDer( ( a.der()-c.val()*b.der() ) /b.val());
    return c;
}

template<typename T>
RallNo<T> operator/(const RallNo<T> &a, const RallNo<T> &b){
    switch ((a.containsX()?1:0)|(b.containsX()?2:0)) {
            // a!=f(x), b!=f(x)
        case 0:
            return RallNo<T>(a.val()/b.val());
            // a=f(x), b!=f(x)
        case 1:
            return div2(a,b.val());
            // a!=f(x), b=f(x)
        case 2:
            return div1(a.val(),b);
    }
    // a=f(x), b=f(x)
    return div3(a,b);
}

// *** UNAERES + ***

template<typename T>
RallNo<T> operator+(const RallNo<T> &a){
    
    RallNo<T> c(a.val());
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(a.der());
    return c;
}

// *** UNAERES - ***

template<typename T>
RallNo<T> operator-(const RallNo<T> &a){
    
    RallNo<T> c(-a.val());
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(-a.der());
    return c;
}

// *** EXP(X) ***

template<typename T>
RallNo<T> exp(const RallNo<T> &a){
    
    RallNo<T> c(exp(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    c.setDer(a.der()*c.val() );
    return c;
}

// *** LN(X) ***

template<typename T>
RallNo<T> log(const RallNo<T> &a){
    
    RallNo<T> c(log(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = 1/a(x)*a'(x)
    c.setDer(1/a.val()*a.der() );
    return c;
}

// *** LOG10(X) ***

template<typename T>
RallNo<T> log10(const RallNo<T> &a){
    
    RallNo<T> c(log10(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = 1/a(x)*a'(x)
    c.setDer(1/a.val()*a.der() );
    return c;
}

// *** SQR(X) ***

template<typename T>
RallNo<T> sqr(const RallNo<T> &a){
    
    RallNo<T> c(a.val()*a.val());
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = 2*a(x)*a'(x)
    c.setDer(2*a.val()*a.der() );
    return c;
}

template<typename T = double, typename ConstType >
RallNo<T> sqr(ConstType &a){
    
    RallNo<T> c(a*a);
    return c;
}

// *** SQRT(X) ***

template<typename T>
RallNo<T> sqrt(const RallNo<T> &a){
    
    RallNo<T> c(sqrt(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = 1/(2*sqrt(a(x)))*a'(x)
    c.setDer(1/( 2*sqrt(a.val()) )*a.der() );
    return c;
}


// *** POW(X,Y) ***


// a=const, b=f(x)
template<typename T, typename ConstType>
RallNo<T> pow1(const ConstType &a, const RallNo<T> &b){
    
    RallNo<T> c(pow(a, b.val()));
    if (!b.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = a^b*ln(a)*b'
    c.setDer(c.val()*log(a)*b.der());
    return c;
}

// a=f(x), b=const
template<typename T, typename ConstType>
RallNo<T> pow2(const RallNo<T> &a, const ConstType &b){
    
    RallNo<T> c(pow(a.val(),b));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = b*a^(b-1)*a'
    c.setDer(b*pow(a.val(), b-1)*a.der());
    return c;
}

template<typename T, typename ConstType>
RallNo<T> pow(const ConstType &a, const RallNo<T> &b){
    return pow1(a, b);
}

template<typename T, typename ConstType>
RallNo<T> pow(const RallNo<T> &a, const ConstType &b){
    return pow2(a, b);
}


template<typename T>
RallNo<T> pow3(const RallNo<T> &a, const RallNo<T> &b){
    RallNo<T> c(pow(a.val(), b.val()));
    //c' = b*a^(b-1)*a' + a^b*ln(a)*b' = g*a' + h*b'
    T g(b.val()*pow(a.val(), b.val()-1));
    T h(c.val()*log(a.val()));
    c.setContainsX();
    c.setDer(g*a.der()+h*b.der());
    return c;
    
}

template<typename T>
RallNo<T> pow(const RallNo<T> &a, const RallNo<T> &b){
    
    switch ((a.containsX()?1:0)|(b.containsX()?2:0)) {
            // a!=f(x), b!=f(x)
        case 0:
            return RallNo<T>(pow(a.val(),b.val()));
            // a=f(x), b!=f(x)
        case 1:
            return pow2(a,b.val());
            // a!=f(x), b=f(x)
        case 2:
            return pow1(a.val(),b);
    }
    // a=f(x), b=f(x)
    return pow3(a,b);
    
}


// *** SIN(X) ***

template<typename T>
RallNo<T> sin(const RallNo<T> &a){
    
    RallNo<T> c(sin(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = cos(a(x))*a'(x)
    c.setDer(cos(a.val())*a.der() );
    return c;
}

// *** COS(X) ***

template<typename T>
RallNo<T> cos(const RallNo<T> &a){
    
    RallNo<T> c(cos(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = -sin(a(x))*a'(x)
    c.setDer(-sin(a.val())*a.der() );
    return c;
}

// *** TAN(X) ***

template<typename T>
RallNo<T> tan(const RallNo<T> &a){
    
    RallNo<T> c(tan(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = 1/cos^2(a(x))*a'(x)
    c.setDer(1/(cos(a.val())*cos(a.val()))*a.der() );
    return c;
}

// *** COT(X) ***

template<typename T>
RallNo<T> cot(const RallNo<T> &a){
    
    RallNo<T> c(1/tan(a.val()));
    if (!a.containsX()) {
        return c;
    }
    //c ist abh. von x
    c.setContainsX();
    //c' = -1/sin^2(a(x))*a'(x)
    c.setDer(-1/(sin(a.val())*sin(a.val()))*a.der() );
    return c;
}





#endif
