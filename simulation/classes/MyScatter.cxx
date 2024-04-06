#include<Riostream.h>
#include<TObject.h>
#include<TMath.h>
#include"MyScatter.h"

ClassImp(MyScatter)

MyScatter::MyScatter():TObject(),
ftheta(0),
fphi(0),
fthetap(0),
fphip(0){
    //DEFALUT CONSTRUCTOR
}

MyScatter::MyScatter(double theta, double phi, double thetap, double phip):TObject(),
ftheta(theta),
fphi(phi),
fthetap(thetap),
fphip(phip){
    //STANDARD CONSTRUCTOR
}

//sono gli angoli nuovi nel sdR del laboratorio
void MyScatter::SetPhi(double* u){
    if( u[1] > 0 ) { //siamo nel I o II quadrante
        if(u[1]/u[0] >0)   fphi= TMath::ATan( u[1]/u[0] );    //I quandrante
        else fphi= TMath::ATan(u[1]/u[0] )+ TMath::Pi(); //II quadrante
    }
    if( u[1] < 0 ) { //siamo nel III o IV quadrante
        if(u[1]/u[0] >0) fphi= TMath::ATan(u[1]/u[0] )+ TMath::Pi(); // III quadrante
        else fphi= TMath::ATan( u[1]/u[0] ) + TMath::Pi()*2; // IV quadrante
    }
}

void MyScatter::SetTheta(double* u){
    ftheta=TMath::ACos(u[2]); //theta è univocamente determinato, in quanto vale tra 0 e Pi
}

void MyScatter::SetPhip(double phip){
    fphip=phip;
}

void MyScatter::SetThetap(double thetap){
    fthetap=thetap;
}

void MyScatter::rotate(double* u){
    double mr[3][3];
    mr[0][0]=-TMath::Sin(fphi);
    mr[1][0]=TMath::Cos(fphi);
    mr[2][0]=0.;
    mr[0][1]=-TMath::Cos(fphi)*TMath::Cos(ftheta);
    mr[1][1]=-TMath::Cos(ftheta)*TMath::Sin(fphi);
    mr[2][1]=TMath::Sin(ftheta);
    mr[0][2]=TMath::Sin(ftheta)*TMath::Cos(fphi);
    mr[1][2]=TMath::Sin(ftheta)*TMath::Sin(fphi);
    mr[2][2]=TMath::Cos(ftheta);  

    double cdp[3];
    cdp[0]=TMath::Sin(fthetap)*TMath::Cos(fphip);
    cdp[1]=TMath::Sin(fthetap)*TMath::Sin(fphip);
    cdp[2]=TMath::Cos(fthetap);
    for(int i=0; i<3; i++){
        u[i]=0.;
        for(int j=0; j<3; j++) u[i] += mr[i][j]*cdp[j];
    }
}

void MyScatter::StartScatter(double theta, double phi, double thetap, double phip){
    ftheta= theta;
    fphi=phi;
    fthetap=thetap;
    fphip=phip;
}

void MyScatter::ClearScatter(){
    ftheta=0;   
    fphi=0;
    fthetap=0;
    fphip=0;
}

