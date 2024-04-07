//ANALISI DATI Z UNIFORME E MOLTEPLICITA' A 50

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

void analisiZUni(){
    gStyle->SetOptFit(1);
    double Zsim;
    double Zrec;
    int mul;
    double res;
    vector<double> Zres;
    double var;

    TFile hfile("../ricostruzione/ricostruzioneUni.root");      //molteplicità fissa e z uniforme
    TTree* tree= (TTree*) hfile.Get("RU2");

    TBranch* b1= tree->GetBranch("Zsim");
    TBranch* b2=tree->GetBranch("Zrec");
    TBranch* b3=tree->GetBranch("mul");
    b1->SetAddress(&Zsim );
    b2->SetAddress(&Zrec);
    b3->SetAddress(&mul);

    int ms=16;
    double Ztrue[]={-22.,-18.,-16.,-15.,-13.,-7.,-4.,-1.,1.,4.,7.,13.,15.,16.,18.,22.};
    double step[]={0.5,0.5,1.,1.5,1.5,1.5,0.5,0.5,0.5,0.5,1.5,1.5,1.5,1.,0.5,0.5};

    int ms2=18;
    double Ztrue2[]={-13.5,-13.,-12.,-11,-9,-8,-7.,-4.,-1.,1.,4.,7.,8.,9.,11.,12.,13.,13.5};
    double step2[]={1.,1.,1.,1.,1.,1.,1.5,1.5,0.5,0.5,0.5,0.5,1.5,1.5,1.,1.,1.,1.,1.,1.};
    
    double resolution[ms];
    double sx[ms];
    double sx2[ms2];
    for(int k=0;k<ms;k++) sx[k]=step[k];
    for(int k=0;k<ms2;k++) sx2[k]=step2[k];
    double efficiency[ms];
    double sy[ms];
    double sy2[ms2];
    double xmin=0.;
    double xmax=0.;
    int nbins;
    double countBin;

    double Nint;
    double Nriv;
    
    TH1D* histo1[ms];
    TCanvas* d[ms];
    TH1D* histo= new TH1D();

//----------------------------EFFICIENZA VS ZTRUE--------------------------------------------------------------
    for(int j=0; j< ms; j++ ){
        Nint=0;
        Nriv=0;
        for( int i=0; i< tree->GetEntries(); i++){
            tree->GetEvent(i);
            if (Zsim >= (Ztrue[j] - step[j]) &&  Zsim<= (Ztrue[j] + step[j])){
                Nint++;
                if( Zrec < 49.) Nriv++;
            }
        }
        cout<<"Nint "<< Nint<<"     Nriv "<<Nriv<<endl;
        if(Nint>0) efficiency[j]= Nriv / Nint;
        if(Nint==0) efficiency[j]=0;
        if(efficiency[j]>0. && efficiency[j]<1.) sy[j]= TMath::Sqrt( efficiency[j]*( 1- efficiency[j] ) / Nint );
        if(efficiency[j]== 1.) sy[j]= 1/Nint;

    }
   

//------------------------------RISOLUZIONE VS ZTRUE----------------------------------------------------------
    TF1* fgaus= new TF1("fgaus", "gaus", xmin,xmax);
    
    for(int j=0; j<ms2; j++){
       
        for(int i=0; i< tree->GetEntries(); i++){
            tree->GetEvent(i);
            if( Zsim >= (Ztrue2[j] - step2[j]) && Zsim <= (Ztrue2[j] + step2[j]) ){
                if(Zrec < 49.){
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
        countBin= TMath::Abs((xmax-xmin)/(0.01)) ; //ogni classe 100 micron
        if(countBin<=1) nbins=1;
        if(countBin>1) nbins=static_cast<int>(countBin);
        histo->SetBins(nbins, xmin, xmax);
        for( int k=0; k< static_cast<int> (Zres.size() ) ; k++) histo->Fill( Zres[k] ); 
        histo->Fit("fgaus"); 
        resolution[j]=fgaus->GetParameter(2)*10000;
        sy2[j]= fgaus->GetParError(2)*10000;
        histo->Reset();
        Zres.clear();
    }

    TGraphErrors* MyPlot= new TGraphErrors( ms, Ztrue, efficiency, sx, sy); 
    MyPlot->GetXaxis()->SetTitle(" Ztrue [cm]");
    MyPlot->GetYaxis()->SetTitle(" Efficiency ");
    MyPlot->SetTitle("Efficienza vs Ztrue uniforme");
    MyPlot->GetYaxis()->SetRangeUser(0., 1.05);
    TCanvas* c1= new TCanvas();
    MyPlot->Draw("ALP");

    TGraphErrors* MyPlot2= new TGraphErrors( ms2, Ztrue2, resolution, sx2, sy2);
    MyPlot2->SetMarkerStyle(50);
    MyPlot2->GetXaxis()->SetTitle(" Ztrue [cm] ");
    MyPlot2->GetYaxis()->SetTitle(" Resolution [#mum]");
    MyPlot2->SetTitle("Risoluzione vs Ztrue uniforme");
    TCanvas* c2=new TCanvas();
    MyPlot2->Draw("ALP");

}


