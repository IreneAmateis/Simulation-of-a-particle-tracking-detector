//ANALISI DATI CON DISTRIBUZIONE ASSEGNATA DI MOLTEPLICITA'

#include <Riostream.h>
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include <TAxis.h>
#include<TF1.h>
#include <TH1D.h>
#include<TFile.h>
#include<TStyle.h>
#include<TGraphErrors.h>
#include<TCanvas.h>
#include<vector>
//ANALISI FILE CON MOLTEPLICITA' ASSSEGNATA

void analisiDistro(){
    double Zsim;
    double Zrec;
    int mul;

    TFile hfile("../ricostruzione/ricostruzioneDistro.root");      //molteplicità da distibuzione assegnata
    TTree* tree= (TTree*) hfile.Get("RD2");

    TBranch* b1= tree->GetBranch("Zsim");
    TBranch* b2=tree->GetBranch("Zrec");
    TBranch* b3=tree->GetBranch("mul");
    b1->SetAddress(&Zsim );
    b2->SetAddress(&Zrec);
    b3->SetAddress(&mul);

    double xmin=0.;
    double xmax=0.;
    int nbins;
    double countBin; 
    double res;
    double sigma;
    double mean;
    vector<double> Zres1;
    vector<double> Zres2;

    int msa=12;
    double multiplicity[]={2,3,4,6,10,13,17,22,27,33,40,47};
    double step[]={1.,1.,1.,1.5,2.5,2.5,2.5,2,2,2,2,2};

    int msb=10;
    double Ztrue[]={-15.5,-13.,-8.,-4.,-1.5,1.5,4.,8.,13.,15.5};
    double step2[]={3., 3., 3., 1.5, 1.5, 1.5, 1.5, 3, 3, 3};

    double sxa[msa]; //errore molteplicità
    double sxb[msb]; //errore Ztrue
    for(int k=0; k<msa;k++) sxa[k]=step[k]; 
    for(int k=0; k<msb;k++) sxb[k]=step2[k];
    double resolution1[msa];
    double resolution2[msb];
    double efficiency2[msb];
    double sy1r[msa]; //errori risoluzione
    double sy2r[msb];
    double sy2e[msb]; //errore efficienza
    double var;

    double Nint;
    double Nriv;

    TFile fileOut("analisiDistro.root","RECREATE");
    TH1D* histo[msa];
    TH1D* histo2r[msb];
    
    char funzione[50];
    char titolo[50];

    //----------RISOLUZIONE VS MOLTEPLICITA' ASSEGNATA-------------------------------------------------
    for(int j=0; j<msa; j++){
        sprintf(funzione, "fgaus%d", j);
        sprintf(titolo, "istogramma ris vs mult  %d", j);
        TF1* fgaus= new TF1("fgaus", "gaus", -0.6, 0.6);
        histo[j]= new TH1D();

        for( int i=0; i< tree->GetEntries(); i++){
            tree->GetEvent(i);
            if ( mul> (multiplicity[j]-step[j]) && mul <(multiplicity[j]+step[j]) ){
                if( Zrec < 49.&& TMath::Abs(Zsim)<=15.9){ 
                    res=Zsim-Zrec;
                    if(i==0){
                        xmax=res;
                        xmin=res;
                    }
                    if( res>xmax) xmax=res;
                    if( res<xmin) xmin=res;
                    Zres1.push_back(res);
                }
            }
        } //distro residui 

        countBin= TMath::Abs ((xmax-xmin)/(0.01)) ; 
        if(countBin<=1) nbins=1;
        if(countBin>1) nbins=static_cast<int>(countBin);
        histo[j]->SetTitle(titolo);
        histo[j]->SetBins(nbins, xmin, xmax);
        for( int k=0; k< static_cast<int> (Zres1.size() ) ; k++) histo[j]->Fill( Zres1[k] );
        
        sigma=histo[j]->GetStdDev(1);
        mean= histo[j]->GetMean(1);
        double Norm= histo[j]->GetMaximum(xmax);
        
        fgaus->SetLineColor(kGreen);
        fgaus->SetParameter(0, Norm);
        fgaus->SetParameter(1, mean);
        fgaus->SetParameter(2, sigma);
        histo[j]->Fit("fgaus","N"); 

        resolution1[j]=fgaus->GetParameter(2)*10000;
        sy1r[j]= fgaus->GetParError(2)*10000;
        Zres1.clear();
    }


//----------------RISOLUZIONE VS ZTRUE, MOLTEPLICITA' DISTRIBUZIONE ASSEGNATA----------------------------
    for(int j=0; j<msb; j++){
        TF1* fgaus2= new TF1("fgaus2", "gaus", -0.6,0.6);
        histo2r[j]=new TH1D(); 
        for(int i=0; i< tree->GetEntries(); i++){
            tree->GetEvent(i);
            if( Zsim >= Ztrue[j] - step2[j] && Zsim <= Ztrue[j] + step2[j]){ //prendi il valore di residui che corrisponde a quel determinato Zvert
                if(Zrec < 49. && TMath::Abs(Zsim)<=15.9){ //Zsim entro 3sigma
                    res=Zsim-Zrec;
                    if(i==0){
                        xmax=res;
                        xmin=res;
                    }
                    if( res>xmax) xmax=res;
                    if( res<xmin) xmin=res;
                    Zres2.push_back(res);
                }
            }
        }
        
        countBin= TMath::Abs ((xmax-xmin)/(0.01)) ;
        if(countBin<=1) nbins=1;
        if(countBin>1) nbins=static_cast<int>(countBin);
        histo2r[j]->SetBins(nbins, xmin, xmax);
        for( int k=0; k< static_cast<int> (Zres2.size() ) ; k++) histo2r[j]->Fill( Zres2[k] );
        
        histo2r[j]->SetFillColor(kCyan-10);
        double max= histo2r[j]->GetMaximum(xmax);
        fgaus2->SetParameter(0, max);
        fgaus2->SetParameter(1, mean);
        fgaus2->SetParameter(2, sigma);
        histo2r[j]->Fit("fgaus2","N");
        resolution2[j]=fgaus2->GetParameter(2)*10000;
        sy2r[j]= fgaus2->GetParError(2)*10000;
        
        Zres2.clear();
    }
//-----------------------------------------------------------------------------------------------------------------

//-------EFFICIENZA VS ZTRUE ---------------------

    for(int i=0; i<msb; i++){
        Nint=0.;
        Nriv=0.;
        for(int j=0; j< tree->GetEntries(); j++){
            tree->GetEvent(j);
            if( Zsim >= Ztrue[i] - step2[i] && Zsim <= Ztrue[i] + step2[i] && TMath::Abs(Zsim)<=15.9){ 
                Nint++;
                if(Zrec < 49.) Nriv++;
            } //events
        } //tree

        if(Nint>0) efficiency2[i]= Nriv / Nint;
        if( Nint==0 ) efficiency2[i]=0;
        if(efficiency2[i]==0. || efficiency2[i]==1.) sy2e[i]=1/Nint;  
        if(efficiency2[i]>0. && efficiency2[i]<1.) sy2e[i]= TMath::Sqrt( efficiency2[i]*( 1- efficiency2[i] ) / Nint );
    }  
    
//---------------------------GRAFICA-----------------------------------

    TGraphErrors* MyPlot1r= new TGraphErrors( msa, multiplicity, resolution1, sxa, sy1r); //resolution vs multiplicity
    MyPlot1r->GetYaxis()->SetTitle(" Resolution [ #mum ] ");
    MyPlot1r->GetXaxis()->SetTitle("Moltiplicity ");
    MyPlot1r->SetTitle("Risoluzione vs Molteplicita' assegnata");
    TCanvas* c1r=new TCanvas("c1r");
    MyPlot1r->Draw("ALP");

    TGraphErrors* MyPlot2r= new TGraphErrors( msb, Ztrue, resolution2, sxb, sy2r); //resolution vs ztrue
    MyPlot2r->GetXaxis()->SetTitle(" Ztrue [cm]");
    MyPlot2r->GetYaxis()->SetTitle(" Resolution [ #mum ] ");
    MyPlot2r->SetTitle("Risoluzione vs Ztrue, con molteplicita' assegnata");
    TCanvas* c2r=new TCanvas("c2r");
    MyPlot2r->Draw("ALP");

    TGraphErrors* MyPlot2e= new TGraphErrors( msb, Ztrue, efficiency2, sxb, sy2e); //efficienza vs Ztrue
    MyPlot2e->GetXaxis()->SetTitle(" Ztrue [cm] ");
    MyPlot2e->GetYaxis()->SetTitle(" efficiency ");
    MyPlot2e->SetTitle("Efficienza vs Ztrue, con molteplicita' asseganta");
    TCanvas* c2e= new TCanvas("c2e");
    MyPlot2e->Draw("ALP");

    MyPlot2e->Write("ZtruevsEff");
    MyPlot1r->Write("resvsmul");
    MyPlot2r->Write("Ztruevsres");
    

    Zres1.clear();
    Zres2.clear();

    fileOut.Close();
}