//
//  BMPLoader.h
//  OdeSolver
//
//  Klasse zum einlesen, bearbeien und schreiben eines BMP-Bildes
//  bearbeiten: schwarz-weiss-Filter, Positionsberechnung der Pixel durch Strahlen im inhomogenen Feld

#ifndef OdeSolver_BMPLoader_h
#define OdeSolver_BMPLoader_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <vector>
#include "OdeSolver.h"

// boost ublas-matrix
#include <boost/numeric/ublas/matrix.hpp>


#pragma pack(2)
//size == 3byte
typedef struct
{
    unsigned char  rgbBlue;          /* Blue value */
    unsigned char  rgbGreen;         /* Green value */
    unsigned char  rgbRed;           /* Red value */
} RGB;

//size == 14byte
typedef struct                       /**** BMP file header structure ****/
{
    unsigned short bfType;           /* Magic number for file */
    unsigned int   bfSize;           /* Size of file */
    unsigned short bfReserved1;      /* Reserved */
    unsigned short bfReserved2;      /* ... */
    unsigned int   bfOffBits;        /* Offset to bitmap data */
} BITMAPFILEHEADER;

//size == 40byte
typedef struct                       /**** BMP file info structure ****/
{
    unsigned int   biSize;           /* Size of info header */
    int            biWidth;          /* Width of image */
    int            biHeight;         /* Height of image */
    unsigned short biPlanes;         /* Number of color planes */
    unsigned short biBitCount;       /* Number of bits per pixel */
    unsigned int   biCompression;    /* Type of compression to use */
    unsigned int   biSizeImage;      /* Size of image data */
    int            biXPelsPerMeter;  /* X pixels per meter */
    int            biYPelsPerMeter;  /* Y pixels per meter */
    unsigned int   biClrUsed;        /* Number of colors used */
    unsigned int   biClrImportant;   /* Number of important colors */
} BITMAPINFOHEADER;

class BMPLoader{
    
private:
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
public:
    std::vector<RGB> pixels;
    
public:
    
    // Bild einlesen
    std::vector<RGB> loadBMP(const char* filename) {
        
        std::ifstream input;
        input.open(filename, std::fstream::in | std::fstream::binary);
        
        if (!input.is_open()) {
            std::cout << "loadBMP: can not open " << filename << "\n";
            return pixels;
        }
        
        input.read((char*)&bfh, sizeof(bfh));
        input.read((char*)&bih, sizeof(bih));
        
        std::cout << "Picture: " << filename << " loaded.\n";
        std::cout << "Width: " << bih.biWidth << ", Height: " << bih.biHeight << ".\n";
        std::cout << "#Pixels: " << bih.biWidth*bih.biHeight << "\n\n";
        
        RGB pix;
        int pixInLine;
        
        for (int y=0; y<bih.biHeight; y++) {
            // wievieltes pixel in einer Zeile (fuer padding)
            pixInLine = 0;
            
            for (int x=0; x<bih.biWidth; x++) {
                
                input.read((char*)&pix, sizeof(pix));
                //                 printf( "PixelR %d: %3d %3d %3d\n", i+1, pix.rgbRed, pix.rgbGreen, pix.rgbBlue );
                
                pixInLine += sizeof(RGB);
                pixels.push_back(pix);
                
            }
            // padding
            if (pixInLine % 4 != 0) {
                pixInLine = 4- (pixInLine%4);
                input.read((char*)&pix, pixInLine);
            }
            
        }
        
        input.close();
        return pixels;
        
    }
    
    // Bild schreiben
    void writeBMP(std::vector<RGB> pxVec, const char* flname){
        
        std::stringstream filename;
        filename << flname << ".bmp";
        
        //testen, ob Datei schon existiert
        int fileCount = 1;
        while (file_exists(filename.str())) {
            filename.str("");
            filename << flname << "(" << fileCount << ").bmp";
            fileCount++;
        }
        std::string fname = filename.str();
        std::ofstream output(fname);
        
        output.write((char*)&bfh, sizeof(bfh));
        output.write((char*)&bih, sizeof(bih));
        
        RGB pix;
        int pixInLine;
        int k = 0;
        for (int y=0; y<bih.biHeight; y++) {
            // wievieltes pixel in einer Zeile (fuer padding)
            pixInLine = 0;
            
            for (int x=0; x<bih.biWidth; x++) {
                
                pix = pxVec[k];
                
                output.write((char*)&pix, sizeof(pix));
                printf( "PixelWrite %d: %3d %3d %3d\n", k+1, pix.rgbRed, pix.rgbGreen, pix.rgbBlue );
                
                pixInLine += sizeof(RGB);
                k++;
                
            }
            // padding
            if (pixInLine % 4 != 0) {
                pixInLine = 4- (pixInLine%4);
                output.write((char*)&pix, pixInLine);
            }
            
        }
        
        
    }
    
    // Filter fuer schwarz-weiss
    std::vector<RGB> makeBlackAndWhite(std::vector<RGB> pixVec){
        
        RGB pix, newPix;
        int k = 0;
        double sum = 0;
        for (int j=0; j<bih.biHeight; j++) {
            
            for (int i=0; i<bih.biWidth; i++) {
                
                pix = pixVec[k];
                sum = pix.rgbRed + pix.rgbGreen + pix.rgbBlue;
                if (sum >= 3.0*255.0/2.0) {
                    newPix.rgbRed = 255;
                    newPix.rgbGreen = 255;
                    newPix.rgbBlue = 255;
                } else {
                    newPix.rgbRed = 0;
                    newPix.rgbGreen = 0;
                    newPix.rgbBlue = 0;
                }
                pixVec[k] = newPix;
                printf( "PixelBW %d: %3d %3d %3d\n", i+1, pix.rgbRed, pix.rgbGreen, pix.rgbBlue );
                
                k++;
                
            }
            
        }
        
        return pixVec;
 
    }
    
    // Positionsberechnung der Pixel durch Strahlen im inhomogenen Feld
    std::vector<RGB> makeGaussPicture(std::vector<RGB> pxVec, std::string &newFileNm){
        
        // *** Werte einstellen ***
        
        //Offset des Feldes vom Ursprung
        double xOffsetGauss = 0.5*bih.biWidth; //0.5*width
        double yOffsetGauss = 0.5*bih.biHeight;
        double zOffsetGauss = -10.0; //-10
        
        //alle Werte verstaerken
        double aFactorGauss = 60.0; //50 //60
        //breiter machen
        double bFactorGauss = 90.0; //40 //60
        
        //z-Startpunkt der Strahlen
        double zFarClippingPlane = -14.0; //-30 //-14
        //z-Ankunfspunkt der Strahlen
        double zNearClippingPlane = -0.5; //-0.5
        
        //Impuls
        double xStartImpulse = 0.0; //0
        double yStartImpulse = 0.0; //0
        double zStartImpulse = 60.0; //18 /33 //80
        
        //Schrittweite fuer die Zeitintegration
        double stepWidth = 0.05;
        
        
        
        //Dateiname
        std::stringstream filenameSS;
        filenameSS << newFileNm << "__x_" << xOffsetGauss << "__y_" << yOffsetGauss << "__z_" << zOffsetGauss << "__aF_" << aFactorGauss << "__bF_" << bFactorGauss << "__zFCP_" << zFarClippingPlane << "__zNCP_" << zNearClippingPlane << "__px_" << xStartImpulse << "__py_" << yStartImpulse << "__pz_" << zStartImpulse;
        newFileNm = filenameSS.str();
        
        
        PicturePlane p1(zNearClippingPlane, zFarClippingPlane);
        RungeKutta4 rk4p("PicturePlane");
        //Functoren fuer die Felder erstellen
        fourth_order_system_gauss_offset_2D<double> gauss_off_2D(xOffsetGauss, yOffsetGauss, zOffsetGauss, aFactorGauss, bFactorGauss);
        fourth_order_system_gauss_offset_3D<double> gauss_off_3D(xOffsetGauss, yOffsetGauss, zOffsetGauss, aFactorGauss, bFactorGauss);
        fourth_order_system_sphere_offset_3D<double> sphere_off_3D(xOffsetGauss, yOffsetGauss, zOffsetGauss, aFactorGauss, bFactorGauss);
        
        //Pixel zuerst in Matrix uebertragen
        boost::numeric::ublas::matrix<RGB> xyPixPlaneORG(bih.biWidth, bih.biHeight); 
        int k = 0;
        for (int y=0; y<bih.biHeight; y++) {
            for (int x=0; x<bih.biWidth; x++) {
                xyPixPlaneORG(x,y) = pxVec[k];
                k++;
                //                printf( "PxPlane %d: %3d %3d %3d\n", k, xyPixPlane(x,y).rgbRed, xyPixPlane(x,y).rgbGreen, xyPixPlane(x,y).rgbBlue );
            }
        }
        
        //neue Matrix fuer transformierte Pixel
        boost::numeric::ublas::matrix<RGB> xyPixPlaneNEW(bih.biWidth, bih.biHeight);
        
        //Hintergrund auf schwarz setzen
        RGB blackPix;
        blackPix.rgbBlue = 0;
        blackPix.rgbGreen = 0;
        blackPix.rgbRed = 0;
        int l = 0;
        for (int y=0; y<bih.biHeight; y++) {
            for (int x=0; x<bih.biWidth; x++) {
                xyPixPlaneNEW(x,y) = blackPix;
                //                printf( "PxPlaneNEW %d: %3d %3d %3d\n", l+1, xyPixPlaneNEW(x,y).rgbRed, xyPixPlaneNEW(x,y).rgbGreen, xyPixPlaneNEW(x,y).rgbBlue );
                l++;
            }
        }
        
        //neuen Ort(x,y) der Pixel nach Transformation durch Strahl berechnen
        state_type_double rayPState3d(6);
        state_type_double rayPState2d(6);
        state_type_double intersectionPoint3d(3);
        state_type_double intersectionPoint2d(3);
        std::cout.precision(2);
        std::cout << "Progress: " << std::fixed << 0.00 << "%" << "\n";
        
        for (int y=0; y<bih.biHeight; y++) {
            
            for (int x=0; x<bih.biWidth; x++) {

                //Start-Bedingungen: [px0, py0, pz0, x0, y0, z0]
                rayPState3d[0] = xStartImpulse;
                rayPState3d[1] = yStartImpulse;
                rayPState3d[2] = zStartImpulse;
                
                // alle x-Werte und y-Werte durchgehen
                rayPState3d[3] = (double)x;
                rayPState3d[4] = (double)y;
                // z-Start wie FCP
                rayPState3d[5] = zFarClippingPlane;

                //neuen x und y-Wert berechnen
                
                //Feld auswaehlen
//                 size_t stepsP = p1.calcNCPIntersectionHermite3d(rk4p, gauss_off_2D, rayPState3d, stepWidth, intersectionPoint3d, false);
               size_t stepsP = p1.calcNCPIntersectionHermite3d(rk4p, sphere_off_3D, rayPState3d, stepWidth, intersectionPoint3d, false);
//               size_t stepsP = p1.calcNCPIntersectionHermite3d(rk4p, gauss_off_3D, rayPState3d, stepWidth, intersectionPoint3d, false);
                
                
                double xSP = intersectionPoint3d[0];
                double ySP = intersectionPoint3d[1];
                double zSP = intersectionPoint3d[2];
                
                // maximale Breite begrenzen
                if (xSP > bih.biWidth){
                    xSP = bih.biWidth-1.0;
                }
                if (xSP < 0.0){
                    xSP = 1.0;
                }
                // maximale Hoehe begrenzen
                if (ySP > bih.biHeight){
                    ySP = bih.biHeight-1.0;
                }
                if (ySP < 0.0){
                    ySP = 1.0;
                }

                //gerundeten x-Wert und y-Wert als neue Koordinaten des Pixels setzen
                // wenn Wert schon vorhanden -> Mittelwert beider Pixel nehmen
                //Pixel schwarz ? -> ueberschreiben, sonst Mittelwert
                if ((xyPixPlaneNEW(x,y).rgbGreen == blackPix.rgbGreen)&& (xyPixPlaneNEW(x,y).rgbBlue == blackPix.rgbBlue) && (xyPixPlaneNEW(x,y).rgbRed == blackPix.rgbRed)) {
                    xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)) = xyPixPlaneORG(x,y);
                } else {
                    xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)).rgbRed = ((xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)).rgbRed + xyPixPlaneORG(x,y).rgbRed))/2.0;
                    xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)).rgbBlue = ((xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)).rgbBlue + xyPixPlaneORG(x,y).rgbBlue))/2.0;
                    xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)).rgbGreen = ((xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)).rgbGreen + xyPixPlaneORG(x,y).rgbGreen))/2.0;
                    
                }
                xyPixPlaneNEW(floor(xSP+0.5), floor(ySP+0.5)) = xyPixPlaneORG(x,y);
                std::cout.precision(0);
                //                std::cout << "(X,Y) = (" << std::fixed << (int)floor(xSP+0.5)+1 << "," << (int)floor(ySP+0.5)+1 << ")\n";
                
            }
            std::cout.precision(2);
            std::cout << "Progress: " << std::fixed << (double)(y+1)/(double)bih.biHeight*100 << "%" << "\n";
            
            
        }
        
        //pixel in Vektor schreiben
        std::vector<RGB> retVec(bih.biHeight*bih.biWidth);
        int m=0;
        for (int y=0; y<bih.biHeight; y++) {
            for (int x=0; x<bih.biWidth; x++) {
                retVec[m] = xyPixPlaneNEW(x,y);
                //                printf( "PixelRet %d: %3d %3d %3d\n", m, retVec[m].rgbRed, retVec[m].rgbGreen, retVec[m].rgbBlue );
                m++;
            }
        }
        
        return retVec;
        
    }
    
    
    
};


#endif
