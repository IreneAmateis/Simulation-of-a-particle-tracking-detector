//ANALISI DATI CON MOLTEPLICITA' UNIFORME

#include <Riostream.h>
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include <TAxis.h>
#include <TH1D.h>
#include<TFile.h>
#include<TStyle.h>
#include<TGraphErrors.h>
#include<TCanvas.h>
#include<vector>
#include<TF1.h>

void analisiMultUni(){
    gStyle->SetOptFit(1);
    double Zsim;
    double Zrec;
    int mul;
    int totv;

    TFile hfile("../ricostruzione/ricostruzione.root");      //molteplicità uniforme e z gaussiano
    TTree* tree= (TTree*) hfile.Get("R2");

    TBranch* b1=tree->GetBranch("Zsim");
    TBranch* b2=tree->GetBranch("Zrec");
    TBranch* b3=tree->GetBranch("mul");
    b1->SetAddress(&Zsim );
    b2->SetAddress(&Zrec);
    b3->SetAddress(&mul);

    double Nint;
    double Nriv;
    double xmin=0.;
    double xmax=0.;
    int nbins;
    double countBin;
    double res;
    vector<double> Zres;
    double var;

    int ms=14;
    double multiplicity[]={3,4,5,7,10,13,17,22,27,33,40,47,62,80}; //è generata uniformemente tra 1 e 80
    double step[]={.5,.5,.5,1.5,1.5,2.,2.5,2.5,2.5,3.,3.,3.,4.,7.};

    double efficiency[ms];
    double resolution[ms];
    double sx[ms];
    double sy[ms];
    double sy2[ms];

    TH1D* histo= new TH1D();

//------------------RISOLUZIONE VS MOLTEPLICITA'

    for(int j=0; j<ms; j++ ){

        for( int i=0; i< tree->GetEntries(); i++){
            tree->GetEvent(i);
            if ( mul<=multiplicity[j]+step[j] && mul>=multiplicity[j]-step[j] ){
                if( Zrec < 49. && TMath::Abs(Zsim)<=5.3){
                    res=Zsim-Zrec;
                    if(i==0){
                        xmax=res;
                        xmin=res;
                    }
                    if( res>xmax) xmax=res;
                    if( res<xmin) xmin=res;
                    Zres.push_back(res);
                }
            }
        }

        countBin= TMath::Abs((xmax-xmin)/(0.05)) ; 
        if(countBin<=1) nbins=1;
        if(countBin>1) nbins=static_cast<int>(countBin);
        histo->SetBins(nbins, xmin, xmax);
        for( int k=0; k< static_cast<int> (Zres.size() ) ; k++) histo->Fill( Zres[k] );
        resolution[j]=histo->GetStdDev()*10000;
        sy2[j]=histo->GetStdDevError()*10000;
        sx[j]=step[j];
        histo->Reset();
        Zres.clear();
        
    }
    TGraphErrors* MyPlotr= new TGraphErrors( ms, multiplicity, resolution, sx, sy2);

    MyPlotr->GetXaxis()->SetTitle("Molteplicity ");
    MyPlotr->GetYaxis()->SetTitle(" Resolution [#mum]");
    MyPlotr->SetTitle("Risoluzione vs Moltiplicita' uniforme");
    MyPlotr->SetMarkerStyle(20);

    TCanvas* c2=new TCanvas();
    MyPlotr->Draw("ALP");

    TFile Ofile("analisiMultUni.root", "RECREATE");
    MyPlotr->Write("resvsmul");

}
