#include<TObject.h>
#include"MyIndex.h"

ClassImp(MyIndex)

MyIndex::MyIndex():TObject(),
fz(0),
fphi(0),
findex(0){
    //DEFAULT CONSTRUCTOR
}

MyIndex::MyIndex(double z, double phi, int index):TObject(),
fz(z),
fphi(phi),
findex(index){
    //STANDARD CONSTRUCTOR
}

void MyIndex::SetData(double z, double phi, int index){
    fz=z;
    fphi=phi;
    findex=index;
}

void MyIndex::SetData(double z, double phi){
    fz= fz+z;
    fphi= fphi+phi;
}