#ifndef MYGEN_H
#define MYGEN_H

#include <TH1F.h>
#include<TRandom3.h>
class MyGen : public TRandom3{

    public:
        
        MyGen();
        double GenZ() ;
        double GenZUni() ;
        double GenX() ;
        double GenY() ;
        int GenMultFixed() ;
        int GenMultUni() ;
        int GenMultDistro();
        double GenTheta();
        double GenPhi();
        double GenThetap();
        void GenNoise(double *u);

    private:
        
        void MyOpenFile();
        TH1F* fheta;
        TH1F* fhmult;
        
    ClassDef(MyGen,1)
}; 

#endif