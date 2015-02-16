//
//  ButcherTable.h
//  OdeSolver
//
//  Butcher-Tableau zur Definition beliebiger expliziter Runge-Kutta-Verfahren mit bis zu 14 Funktionsauswertungen

#ifndef OdeSolver_ButcherTable_h
#define OdeSolver_ButcherTable_h

struct ButcherTable{
    
    double
    a21=0.0,
    a31=0.0, a32=0.0,
    a41=0.0, a42=0.0, a43=0.0,
    a51=0.0, a52=0.0, a53=0.0, a54=0.0,
    a61=0.0, a62=0.0, a63=0.0, a64=0.0, a65=0.0,
    a71=0.0, a72=0.0, a73=0.0, a74=0.0, a75=0.0, a76=0.0,
    a81=0.0, a82=0.0, a83=0.0, a84=0.0, a85=0.0, a86=0.0, a87=0.0,
    a91=0.0, a92=0.0, a93=0.0, a94=0.0, a95=0.0, a96=0.0, a97=0.0, a98=0.0,
    a101=0.0, a102=0.0, a103=0.0, a104=0.0, a105=0.0, a106=0.0, a107=0.0, a108=0.0, a109=0.0,
    a111=0.0, a112=0.0, a113=0.0, a114=0.0, a115=0.0, a116=0.0, a117=0.0, a118=0.0, a119=0.0, a1110=0.0,
    a121=0.0, a122=0.0, a123=0.0, a124=0.0, a125=0.0, a126=0.0, a127=0.0, a128=0.0, a129=0.0, a1210=0.0, a1211=0.0,
    a131=0.0, a132=0.0, a133=0.0, a134=0.0, a135=0.0, a136=0.0, a137=0.0, a138=0.0, a139=0.0, a1310=0.0, a1311=0.0, a1312=0.0,
    a141=0.0, a142=0.0, a143=0.0, a144=0.0, a145=0.0, a146=0.0, a147=0.0, a148=0.0, a149=0.0, a1410=0.0, a1411=0.0, a1412=0.0, a1413=0.0,
    
    b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0, b9=0.0, b10=0.0, b11=0.0, b12=0.0, b13=0.0, b14=0.0,
    c1=0.0, c2=0.0, c3=0.0, c4=0.0, c5=0.0, c6=0.0, c7=0.0, c8=0.0, c9=0.0, c10=0.0, c11=0.0, c12=0.0, c13=0.0, c14=0.0;
    
};


#endif
