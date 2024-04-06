#include<Riostream.h>
#include<TMath.h>
#include<TTree.h>
#include<TFile.h>
#include<TClonesArray.h>
#include<TBranch.h>
#include<TRandom3.h>
#include"MyGen.h"
#include"MyInt.h"
#include"MyScatter.h"
#include"MyIndex.h"
#include<TStopwatch.h>

double const Rb = 3.; //cm
double const R1 = 4.; //cm
double const R2 = 7.; //cm
double const H = 27.; //cm

void MySimulation( int Nevents=100000, bool multScatt= true, unsigned int seed=78540){

//-------------------------------------------------------------------------------------------modificare qui per scegliere il tipo di simulazione

    //TFile hfile("simulazioneUni.root","RECREATE");                                        //molteplicità fissa e Z uniforme
    //TTree *tree = new TTree("TU","TTree di simulazione 2");

    TFile hfile("simulazioneDistro.root","RECREATE");                                    //molteplicità da distribuzione assegnata e Z gaussiano
    TTree *tree = new TTree("TD","TTree di simulazione Distr. Assegnata");

    //TFile hfile("simulazione.root", "RECREATE");                                           //molteplicità uniforme e Z gaussiano
    //TTree* tree= new TTree("T", "TTree di simulazione");
//---------------------------------------------------------------------------------------------

    TStopwatch timer;

    typedef struct{
        double X,Y,Z;
        int mult;
    } VTX;

    static VTX point;

    delete gRandom;
    MyGen* gen= new MyGen();
    gRandom=gen;
    gen->SetSeed(seed);

    MyInt* tInterz= new MyInt();
    MyScatter* scatter= new MyScatter();

    double* u= new double [3]; //array per il multiple scattering

    double noise1[]={0,0}; //array per i punti di noise
    double noise2[]={0,0};

    TClonesArray* ptrhits1= new TClonesArray("MyIndex", 200);
    TClonesArray &hits1= *ptrhits1;
    TClonesArray* ptrhits2= new TClonesArray("MyIndex", 200);
    TClonesArray &hits2= *ptrhits2;
    tree->Branch("Vertice",&point.X, "X/D:Y:Z:mult/I");
    tree->Branch("Layer1", &ptrhits1);
    tree->Branch("Layer2", &ptrhits2);

    timer.Start();
    int tothits1;
    int tothits2;
    double phi;
    double theta;
    double phip;
    double thetap;
    double phip1;
    double thetap1;

    for(int Nvertex=0; Nvertex< Nevents; Nvertex++){

        for(int i=0;i<3;i++) u[i]=0; //settiamo a zero il puntatore u

//------------------------------------------------------------------------------VERTICE - modificare qui per scegliere il tipo di simulazione
        point.X= gen->GenX();
        point.Y = gen->GenY();
        point.Z= gen->GenZ();
        //point.Z = gen->GenZUni();
//------------------------------------------------------------------------------MOLTEPLICITA' - modificare qui per scegliere il tipo di simulazione
        //point.mult = gen->GenMultUni();                                                                                 
        //point.mult = gen->GenMultFixed();                                               
        point.mult=gen->GenMultDistro();                                              
//----------------------------------------------------------------------------------------
        
        tothits1=0; //numero di intersezioni per layer per vertice generato
        tothits2=0;//numero di intersezioni per layer per vertice generato

        cout<<"interazione numero " << Nvertex+1<<" molteplicità "<< point.mult <<endl ;
        //<< "VERTICE    X   " << point.X << "  Y "<< point.Y << "   Z  "<< point.Z << endl;

        for( int i=0; i<point.mult ; i++){

            phi = gen->GenPhi(); 
            theta = gen->GenTheta();
            tInterz->StartPoint(point.X, point.Y, point.Z); //inizializza all'oggetto di tInterz la posizione del vertice
            tInterz->SetPhi(phi);
            tInterz->SetTheta(theta);
           
        //-----------BEAM PIPE-------------;
            
            tInterz->SetX(Rb, u, false);
            tInterz->SetY(Rb, u, false);
            tInterz->SetZ(Rb, u, false);

        //    cout<< "INTERSEZIONE  X " << tInterz->GetX()<< " Y "<<  tInterz->GetY ()<< "  Z " <<tInterz->GetZ() << endl;

            //--------------multiple scattering con beam pipe-----------------
            if(multScatt){  
               phip=   gen->GenPhi();    //estratti nel sdR'
               thetap= gen->GenThetap(); //estratti nel sdR'  
               scatter->StartScatter(theta,phi,thetap,phip); 
               scatter->rotate(u); //popola u con la nuova direzione
               scatter->SetTheta(u); //calcoliamo theta nuovo nel sdr LAB
               scatter->SetPhi(u);   //calcoliamo phi nuovo nel srd LAB
               tInterz->SetTheta( scatter->GetTheta() );   //impostiamo i nuovi angoli all'oggetto MyInt
               tInterz->SetPhi( scatter->GetPhi() );      //impostiamo i nuovi angoli all'oggetto MyInt
            }

        //-----------LAYER 1-----------;

            tInterz->SetX(R1, u , multScatt);
            tInterz->SetY(R1, u , multScatt);
            tInterz->SetZ(R1, u , multScatt);

            if( tInterz->CheckZ()==1 ) {

                new( hits1[tothits1] )MyIndex( tInterz->GetZ() , tInterz->GetPhi() , i);
                tothits1++;
                
                //-----------------multiple scatter con layer 1------------------------------
                if(multScatt){
                    phip1=   gen->GenPhi(); 
                    thetap1= gen->GenThetap();
                    scatter->SetThetap(thetap1);
                    scatter->SetPhip(phip1);
                    scatter->rotate(u);
                    scatter->SetTheta(u);
                    scatter->SetPhi(u);
                    tInterz->SetTheta( scatter->GetTheta() );
                    tInterz->SetPhi( scatter->GetPhi() );
                }
            }    
            
        //-----------LAYER 2---------------;

            tInterz->SetX(R2, u, multScatt);
            tInterz->SetY(R2, u, multScatt);
            tInterz->SetZ(R2, u, multScatt);

            if( tInterz->CheckZ() ==1 ) {
                new( hits2[tothits2] ) MyIndex( tInterz->GetZ() ,tInterz->GetPhi() , i) ; 
                tothits2++;
            }

            tInterz->MyClear();
            scatter->ClearScatter();

        } //molteplicità 

        //generazione dei punti di noise

        if((Nvertex%5000)==0.){
            gen->GenNoise(noise1);
            gen->GenNoise(noise2);
            new( hits1[tothits1] ) MyIndex( noise1[0] , noise1[1] , Nvertex);
            new( hits2[tothits2] ) MyIndex( noise2[0] , noise2[1] , Nvertex);
            tothits1++;
            tothits2++;
        }
    
        tree->Fill();
        ptrhits1->Clear();
        ptrhits2->Clear();

    } //vertice

    hfile.Write();
    hfile.Close();
    timer.Stop();
    timer.Print();

    delete gen;
    delete tInterz;
    delete ptrhits1;
    delete ptrhits2;
    delete scatter;
   

} //macro
