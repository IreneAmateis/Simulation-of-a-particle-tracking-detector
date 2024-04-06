#ifndef MYINT_H
#define MYINT_H

//classe che calcola le intersezioni, ossia se la particella intercetta o no i vari layer

class MyInt : public TObject{
    public:
        MyInt();
        
        double GetX() const {return fX;}
        double GetY() const {return fY;}
        double GetZ() const {return fZ;}
        double GetTheta() const {return fTheta;}
        double GetPhi() const {return fPhi;}

        void StartPoint(double x, double y, double z);
        bool CheckZ();
    
        void SetX(double R, double* u, bool ms);
        void SetY(double R, double* u, bool ms);
        void SetZ(double R, double* u, bool ms);
        void SetTheta(double theta);
        void SetPhi(double phi);
        void MyClear();

    private:
        double fX;
        double fY;
        double fZ;

        double fTheta;
        double fPhi;

        double tInt(double x, double y, double z, double f, double t, double R);
        double tInt(double x, double y, double z, double* u, double R );
        
        double Calc1(double f, double t);
        double Calc2(double f, double t);
        double Calc3(double t);

    ClassDef(MyInt,1)
};

#endif