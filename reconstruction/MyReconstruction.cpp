//macro di ricostruzione

#include <Riostream.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TFile.h>
#include <TStopwatch.h>
#include <TCanvas.h>
#include <vector>
#include<TEfficiency.h>
#include "MyIndex.h"
#include "MyRec.h"

const double R1=4.; //cm
const double R2=7.; //cm
const double dPhiMax=0.006; //rad

void MyReconstruction(unsigned int seed=1234){
    
    //lettura file di input
    typedef struct{ 
        double X,Y,Z;
        int mult;
    } VTX;

    static VTX point;

    TClonesArray* hits1= new TClonesArray("MyIndex", 300);
    TClonesArray* hits2= new TClonesArray("MyIndex", 300);
//--------------LETTURA FILE---------------------------------------------------------------------modificare qui per scegliere il tipo di simulazione
    //TFile hfile("../simulazione/simulazioneUni.root");       //molteplicità fissata ad un valore (50)
    //TTree* tree= (TTree*) hfile.Get("TU");

    //TFile hfile("../simulazione/simulazione.root");         //distribuz di molteplicità uniforme tra 1 e 80 
    //TTree* tree= (TTree*) hfile.Get("T");
    
    TFile hfile("../simulazione/simulazioneDistro.root");    //distribuz assegnata di molteplicità 
    TTree* tree= (TTree*) hfile.Get("TD");
//----------------------------------------------------------------------------------------
    TBranch* b1= tree->GetBranch("Vertice");
    TBranch* b2= tree->GetBranch("Layer1");
    TBranch* b3= tree->GetBranch("Layer2");
    b1->SetAddress( &point.X );
    b2->SetAddress(&hits1);
    b3->SetAddress(&hits2);

    vector<MyIndex> L1,L2;    
    vector<double> vectorVtx;
    vector<int>vectorBin; //vector con elementi gli indici dei bin con maggior entries
    vector<double> vectorCentroids;

    //salvataggio su file di output
    double Zrec;
    double Zsim;
    int mul;
    int totv;
    //double epsilon;
    bool w; //per i grafici di efficienza

//---------------------------------------------------------------------------------------modificare qui per scegliere il tipo di simulazione
    //TFile rfile("ricostruzioneUni.root","RECREATE");                     //z uniforme molteplicità fissa
    //TTree *treeO = new TTree("RU2","TTree di ricostruzione Uni");

    //TFile rfile("ricostruzione.root","RECREATE");                        //molteplicità con distribuz uniforme , z gaussiano
    //TTree *treeO = new TTree("R2","TTree di ricostruzione");

    TFile rfile("ricostruzioneDistro.root","RECREATE");                     //molteplicità con distribuz assegnata , z gaussiano
    TTree *treeO = new TTree("RD2","TTree di ricostruzione Distr. Assegnata");
//-----------------------------------------------------------------------------------------    

    treeO->Branch("mul",&mul,"mul/I");
    treeO->Branch("Zsim",&Zsim,"Zsim/D");
    treeO->Branch("Zrec",&Zrec,"Zrec/D");
    double xmin;
    double xmax;
    int nbins;
    double countBin; //per calcolare il numero di bin
    int count;      //numero di entries in un bin
    int countMax;   //numero di entries nel bin più popolato
    
    delete gRandom;
    MyRec* rec= new MyRec();
    gRandom=rec;
    rec->SetSeed(seed); 

    TH1D* histo= new TH1D();

//---debug---contare il numero di eventi per i diversi casi possibili
    double zero=0;
    double uno=0;
    double due=0;
    int duegruppi=0; //quante volte ci sono più di due gruppi equipollenti
    int problems=0;

    TStopwatch timer;
    timer.Start();

//-------------------------------------------------------------grafico effienza--------------------modificare qui per scegliere il tipo di simulazione
    
    TEfficiency* eff_mult= new TEfficiency("eff","molteplicita' assegnata;mul;#epsilon",15,0,60);    //multDistro    
    //TEfficiency* eff_mult= new TEfficiency("eff","molteplicita' uniforme;mul;#epsilon",20,0,80); //multUni
//--------------------------------------------------------------------------------------------------------------------------------

    for(int ev=0; ev<tree->GetEntries(); ev++){

        totv=0; //conta quanti vertici ho calcolato per evento tramite le tracklet
        nbins=0;      

        tree->GetEvent(ev);

        for(int i=0; i<hits1->GetEntries(); i++){ 
            MyIndex* layer1= (MyIndex*) hits1->At(i); //creaz di Cloni di hits 1
            layer1->SetData(rec->GenZr(), rec->GenPhir(R1)); //smearing su L1
            L1.push_back(*layer1);
        } 

        for(int j=0; j<hits2->GetEntries(); j++) {
            MyIndex* layer2= (MyIndex*) hits2->At(j);//creaz di Cloni di hits 2
            layer2->SetData(rec->GenZr(), rec->GenPhir(R2)); //smearing su L2
            L2.push_back(*layer2);
        }
        
        mul=point.mult;
        Zsim=point.Z;
        double Zvertex;
        TAxis* Xaxis=histo->GetXaxis();

        for(int l1=0; l1< static_cast<int> (L1.size() ); l1++){ //loop layer 1
            for(int l2=0; l2< static_cast<int> ( L2.size() ) ; l2++){ //loop layer 2
                double dPhi=TMath::Abs(  ( L1[l1].GetPhi() )- ( L2[l2].GetPhi() ) ) ;

                if(dPhi< dPhiMax){
                    Zvertex=rec->TrackZ( L1[l1].GetZ() , R1 , L2[l2].GetZ() , R2 ) ; //tracklets
                    vectorVtx.push_back(Zvertex);
                    if(l1==0){
                        xmax= Zvertex;
                        xmin= Zvertex;
                    }
                    if( Zvertex >= xmax) xmax=Zvertex;
                    if( Zvertex <= xmin) xmin=Zvertex;
                    totv++; 
                } //if dPhi           
            } // l2
        }//l1

        countBin=TMath::Abs( ( xmax-xmin) /0.2 ); 
        if(countBin<=1) nbins=1;
        if(countBin>1) nbins=static_cast<int>(countBin);

        histo->SetBins( nbins, xmin , xmax );
        for(int el=0;el<totv;el++) histo->Fill(vectorVtx[el]); //istogramma delle tracklets   

        //-----------------------------------------------------------------------------------------------
        
        count=0;
        countMax=0;

        if(totv==0) { //non ci sono trackelts
            Zrec=50.; //non ricostruisce se non trova eventi
            w=false;
            zero++; //numero di volte in cui questo accade
        }
        if( totv==1) { //una sola tracklets
            Zrec= vectorVtx[0]; //se ha un solo evento ricostruisce quello 
            w=true;
            uno++; //numero di volte in cui questo accade
        }

        if(totv>1){ //ho almeno due tracklets
            due++;
            for(int ibin=1; ibin <= nbins; ibin++){ //cicliamo sui bin
                count= histo->GetBinContent(ibin); 

                if( count == countMax) vectorBin.push_back( ibin );

                if( count > countMax){
                    vectorBin.clear(); 
                    countMax = count;   
                    vectorBin.push_back( ibin );                     
                } 

            }//ibin

            //calcolo dei centroidi per ogni binMax selezionato

            int index; //indice di un bin
            double center;
            double center2;
            double sum;
            double sum2;
            double e1; //counts di index
            double c1; //centro di index
            double height;
            double height2;
            double weight;
            double weight2;
            
            double centroide;
            double totsum;
            double totsumMax=0;

            //------------algoritmo di somma-calcolo della media pesata sui bin contigui--------------
            
            for( int q=0; q< static_cast<int>(vectorBin.size()) ; q++){ //cicli sui binMax
                
                index= vectorBin[q];
                weight=0;
                weight2=0;
                sum=0;
                sum2=0;
                e1=histo->GetBinContent(index);
                c1=histo->GetBinCenter(index);

                for(int s= index+1; s <=nbins; s++){ //parti dall'index-esimo(escluso) bin e vai a destra
                    height= histo->GetBinContent(s) ; 
                    if(height > 0 ){ //se non è vuoto
                        center=histo->GetBinCenter(s); 
                        sum= sum + center*height;
                        weight= weight + height;
                    } 
                    else s= nbins+1; //esci dal for
                }

                for( int r= index-1; r> 0; r--){ 
                    height2= histo->GetBinContent(r);
                    if( height2>0){
                        center2= histo->GetBinCenter(r);
                        sum2= sum2+ center2*height2;
                        weight2= weight2 + height2;
                    }
                    else r=0; //esci dal for
                }

                centroide= (sum+sum2+e1*c1)/(weight+weight2+e1); //sommiamo i contributi a destra e a sinistra del picco max
                totsum= weight+weight2+e1; // devo vedere quale gruppo ha il peso maggiore

                if( totsum==totsumMax) vectorCentroids.push_back(centroide); //centroidi con numero di elementi uguali all'interno della distribuzione dell'histo

                if( totsum > totsumMax){                  
                    vectorCentroids.clear();
                    totsumMax= totsum;
                    vectorCentroids.push_back(centroide);
                }
            }//for binMax  

//---------------debug-------------------------
            /*
            if( static_cast<int>(vectorCentroids.size()) == 0 ) {
                cout << "VECTOR CENTROIDS HA DIM ZERO " << endl;
                problems++;
                }
            */


            //-------calcolo di Zvertice             

            int sizeCent =static_cast<int>(vectorCentroids.size());
            if(sizeCent>2) duegruppi++; //debug

            sort( vectorCentroids.begin(), vectorCentroids.end() );  //così posso confrontarne due per volta in ordine

            double centro=100.; //valore di controllo 

            if(sizeCent==1) centro=vectorCentroids[0];
            if(sizeCent>1){
                for(int c=0; c<sizeCent;c++){                
                    if(vectorCentroids[c]==vectorCentroids[c+1]) centro=vectorCentroids[c];
                }//altrimenti centro rimane uguale a 100
            }

            if(centro>99.){ 
                Zrec=50.; //può succedere che i gruppi siano tutti distanti tra loro, non si ricostruisce      
                w=false;
            }

            int nmean=0;
            double somma=0.;
            
            if(centro<99.){
                for( int k=0; k< totv; k++){               
                    if( TMath::Abs(vectorVtx[k] - centro) <=0.15 ){
                        nmean=nmean+1;
                        somma=somma+vectorVtx[k];
                    }
                }
                
                if(nmean==0) {
                    Zrec=50.; 
                    w=false;
                }
                if(nmean>0) {
                    Zrec=somma/nmean;
                    w=true;
                }
            }
        } //totv>1
//------------------------------------------------------------------------

        cout<<"Il vertice ricostruito per l'evento "<< ev+1 << " è: "<< Zrec <<endl;

        eff_mult->Fill(w, mul);

        L1.clear();
        L2.clear();
        vectorVtx.clear();
        vectorCentroids.clear();
        vectorBin.clear();
        treeO->Fill();
        histo->Reset();
    }//eventi
    

//-------------------------------------------------file di output per i grafici----------modificare qui per scegliere il tipo di simulazione
    TFile Ofile("effDistro.root", "RECREATE"); 
    //TFile Ofile("effMultUni.root", "RECREATE");

//-----Graphics----------------
    TCanvas* c2= new TCanvas();
    c2->cd();
    eff_mult->Draw("ALP");
    eff_mult->Write("Effvsmul");
    Ofile.Close();


//----debug---------------------
    //cout<<"totv=0 "<<zero<<endl<<"totv=1 "<<uno<<endl;
    //cout<<"più di due gruppi equipollenti: "<<duegruppi<<endl;
    
    timer.Stop();
    timer.Print();
    rfile.Write();
    rfile.Close();
    hfile.Close();

}//macro