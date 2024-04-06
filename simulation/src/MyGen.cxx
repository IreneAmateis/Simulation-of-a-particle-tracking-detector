#include<Riostream.h>
#include<TMath.h>
#include<TFile.h>
#include<TH1F.h>
#include<TRandom3.h>
#include"MyGen.h"

ClassImp(MyGen)

MyGen::MyGen(): TRandom3(),
fheta(NULL),
fhmult(NULL){
    MyOpenFile();
    }

void MyGen::MyOpenFile(){
    TFile F("heta2.root");
    fheta =(TH1F*)F.Get("heta2");
    fheta->SetDirectory(0); 
    fhmult =(TH1F*)F.Get("hmul");
    fhmult->SetDirectory(0);
    F.Close();
}

double MyGen::GenZ() {
    double sigma=5.3; //cm
    double u= Gaus(0, sigma);
    return u;
}

double MyGen::GenZUni(){
    double u=-30 + 60*Rndm();
    return u;
}

double MyGen::GenX() {
    double u= Gaus(0, 0.01);
    return u;
}

double MyGen::GenY() {
    double u= Gaus(0, 0.01);
    return u;
}

int MyGen::GenMultFixed() {
    int mult= 50; 
    return mult;
}

int MyGen::GenMultUni() {
    int mult=1 + 80*Rndm();
    return mult;
}

int MyGen::GenMultDistro() {
    int mult= fhmult->GetRandom();
    return mult;
}

double MyGen::GenPhi() {
    double u=2.*(TMath::Pi())*Rndm();
    return u;
}

double MyGen::GenTheta() {
    double eta;
    eta =fheta->GetRandom();
    return 2*TMath::ATan(TMath::Exp(-eta));
}

double MyGen::GenThetap(){
    double u=Gaus(0, 0.001*(TMath::Sqrt(2)) ); //in radianti
    return u;
}

void MyGen::GenNoise(double *u){
    u[0]=-13.5 + 27*Rndm(); //coordinata z
    u[1]=2.*(TMath::Pi())*Rndm(); //angolo phi
}
