#ifndef MYINDEX_H
#define MYINDEX_H

class MyIndex: public TObject{
    public:
    MyIndex();
    MyIndex(double z, double phi, int index);

    void SetData(double z, double phi, int index);
    void SetData(double z, double phi);
    
    double Getz(){return fz;} //usate per debug 
    double Getphi(){return fphi;}
    int Getindex(){return findex;}

    private:
    double fz;
    double fphi;
    int findex;


    ClassDef(MyIndex,1)
};

#endif