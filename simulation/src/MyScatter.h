#ifndef MYSCATTER_H
#define MYSCATTER_H

class MyScatter: public TObject{
    public:
    MyScatter();
    MyScatter( double theta, double phi, double thetap, double phip);
    
    void StartScatter(double theta, double phi, double thetap, double phip);
    void ClearScatter();
    void SetPhi(double* u);
    void SetTheta(double* u);
    void SetPhip(double phip);
    void SetThetap(double thetap);

    double GetTheta(){ return ftheta;}
    double GetPhi(){return fphi;}
    void rotate(double *u);
   
    private:
    double ftheta;
    double fphi;
    double fthetap;
    double fphip;

    ClassDef(MyScatter,1)
};

#endif

