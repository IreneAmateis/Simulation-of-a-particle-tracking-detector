#include "MyInt.h"
#include<Riostream.h>
#include<TMath.h>


ClassImp(MyInt)

    MyInt::MyInt(): TObject(),
fX(0.),
fY(0.),
fZ(0.),
fTheta(0.),
fPhi(0.){
    //DEFAULT CONSTRUCTOR
}

void MyInt::StartPoint(double x, double y, double z){
    fX=x;
    fY=y;
    fZ=z;
}

void MyInt::SetX(double R, double* u, bool ms ){
    if( ms ) fX= fX + u[0]*tInt(fX, fY, fZ, u, R );
    else fX= fX+ Calc1(fPhi, fTheta)*tInt( fX, fY, fZ, fPhi, fTheta, R);
}

void MyInt::SetY(double R, double* u, bool ms){
    if( ms ) fY= fY + u[1]*tInt(fX, fY, fZ, u, R );
    else fY= fY+ Calc2(fPhi, fTheta)*tInt( fX, fY, fZ, fPhi, fTheta, R);
}

void MyInt::SetZ(double R, double* u, bool ms){
    if ( ms ) {
        double result= tInt(fX, fY, fZ, u, R );
        fZ= fZ + u[2]*result; }
    else fZ= fZ+ Calc3(fTheta)*tInt( fX, fY, fZ, fPhi, fTheta, R);

}

void MyInt::SetTheta( double theta){
    fTheta= theta;
}
void MyInt::SetPhi(double phi){
    fPhi= phi;
}

void MyInt:: MyClear(){
    fX=0;
    fY=0;
    fZ=0;
    fPhi=0;
    fTheta=0;
}

//------------------------------------------PRIVATE METHODS-----------------------
double MyInt::Calc1(double f, double t){
    double cosf= TMath::Cos(f);
    double sint= TMath::Sin(t);
    return cosf*sint;
}

double MyInt::Calc2(double f, double t){
    double sinf= TMath::Sin(f);
    double sint= TMath::Sin(t);
    return sint*sinf;
}

double MyInt::Calc3(double t){
    return TMath::Cos(t);
}

double MyInt::tInt(double x, double y, double z, double f, double t, double R){
    double result;
    double d= ( (x*Calc1(f,t)+y*Calc2(f,t) ) *( x*Calc1(f,t)+y*Calc2(f,t) ) ) - ( Calc1(f,t)*Calc1(f,t) + Calc2(f,t)*Calc2(f,t) )*( x*x + y*y - R*R );
    double dr=TMath::Sqrt(d);
    if(dr < 1e-10 ) dr=0;
    double t1= -( (x*Calc1(f,t)+y*Calc2(f,t)) + dr ) / ( Calc1(f,t)*Calc1(f,t) + Calc2(f,t)*Calc2(f,t) );
    double t2= -( (x*Calc1(f,t)+y*Calc2(f,t)) - dr ) / ( Calc1(f,t)*Calc1(f,t) + Calc2(f,t)*Calc2(f,t) );

    if(t1>0) result=t1;
    if(t2>0) result=t2;

    return result;
}
    
//-----------------------------------------------------------------------SCATTERING--------------
double MyInt::tInt(double x, double y, double z, double* u, double R ){
    double result;
    double d= ( (x*u[0] + y*u[1] ) *( x*u[0]+y*u[1] ) ) - ( u[0]*u[0] + u[1]*u[1] )*( x*x + y*y - R*R );
    double dr=TMath::Sqrt(d);
    if(dr < 1e-10 ) dr=0;
    double t1= -( (x*u[0] + y*u[1] ) + dr ) / ( u[0]*u[0] + u[1]*u[1] );
    double t2= -( (x*u[0] + y*u[1] ) - dr ) / ( u[0]*u[0] + u[1]*u[1] );

    if(t1>0) result=t1;
    if(t2>0) result=t2;
    return result; 
}
//-------------------------------------------------------------------------------------------------

bool MyInt::CheckZ(){
    if( TMath::Abs(fZ) <= 13.5 ) return 1;
    else return 0;
}