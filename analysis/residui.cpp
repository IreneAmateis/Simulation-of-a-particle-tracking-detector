// file clone di analisi.cpp
// in questo file provo delle modifiche al file utilizzato come analisi dati 
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


bool CheckMult(int mult);
    
void residui(){
    //gStyle->SetOptFit(1);
    double Zsim;
    double Zrec;
    int mul;

    TFile hfile("../ricostruzione/ricostruzione.root");        //molteplicità uniforme
    TTree* tree = (TTree*) hfile.Get("R2");

//--------------------------------------------------------------------------------------------------
 
    //char titolo[] = "residui molteplicita' uniforme - con tutte le molteplicita'";
    char titolo[]= " residui con molteplicita' uniforme nell'intervallo 63.5-64.5" ;
    TBranch* b1= tree->GetBranch("Zsim");
    TBranch* b2=tree->GetBranch("Zrec");
    TBranch* b3=tree->GetBranch("mul");
    b1->SetAddress(&Zsim );
    b2->SetAddress(&Zrec);
    b3->SetAddress(&mul);

    double xmin=0.;
    double xmax=0.;
    int nbins;
    double res;
    vector<double> Zres;
    vector<double> Zres2;

    TH1D* histo= new TH1D();

    for( int i=0; i< tree->GetEntries() ; i++){
        tree->GetEvent(i);
        if ( CheckMult(mul) ){
            if( Zrec < 49.){
                res=Zsim-Zrec;
                if( res>xmax) xmax=res;
                if( res<xmin) xmin=res;
                Zres.push_back(res);
            }
        }
    }

    nbins= static_cast<int> ((xmax-xmin)/(0.01)) ;
    histo->SetBins(nbins, xmin, xmax);
    histo->SetTitle("istogramma dei residui ");
    for( int j=0; j< static_cast<int> (Zres.size() ) ; j++)  histo->Fill( Zres[j] );
    histo->SetFillColor(kRed-10);

    double mean=histo->GetMean(1);
    double sigma= histo->GetStdDev(1);
    cout << " mean " << mean << " sigma " << sigma <<endl;
    TAxis* Xaxis1= histo->GetXaxis();
    Xaxis1->SetRangeUser(mean-0.5, mean +1); //riscalo il grafico e mi ricalcolo mean e sigma
    mean=histo->GetMean(1);
    sigma= histo->GetStdDev(1);
    cout << " mean " << mean << " sigma " << sigma <<endl;
    //cout<<"ASIMMETRIA "<< histo->GetSkewness(1)<<endl;

//---------------------------------------SELEZIONE 3 SIGMA ----------------------------> toglie tutto il fondo delle code
    TH1D* histo2= new TH1D();
    double xmin2= mean-3*sigma;
    double xmax2 = mean+ 3*sigma;
    int nbins2= (xmax2-xmin2)/0.005;
    for( int i=0; i< tree->GetEntries() ; i++){
        tree->GetEvent(i);
        if ( CheckMult(mul) && Zrec<49. ){
        res=Zsim-Zrec;
        if( res > mean-3*sigma && res < mean + 3*sigma ) Zres2.push_back(res);
        }
    }

    histo2->SetBins( nbins2, xmin2, xmax2);
    histo2->SetTitle(titolo);
    for( int j=0; j< static_cast<int> (Zres2.size() ); j++) histo2->Fill( Zres2[j]);
    double mean2=histo2->GetMean(1);
    double sigma2= histo2->GetStdDev(1);
    histo2->SetOption("E");

    histo2->SetMarkerStyle(kFullCircle);
    cout<<"ASIMMETRIA ISTOGRAMMA 3SIGMA "<< histo2->GetSkewness(1)<<endl;
    

    TF1* fgaus= new TF1("fgaus", "gaus", xmin, xmax);
    fgaus->SetParameter(1, mean);
    fgaus->SetParameter(2, sigma);
    histo->Fit("fgaus", "R");
    //cout <<"fgaus : "<< "Chi^2: " << fgaus->GetChisquare() << ", number DoF: " << fgaus->GetNDF() << ", (Probability: " << fgaus->GetProb() << ")." << endl;


    TF1* fgaus2= new TF1("fgaus2", "gaus", xmin2, xmax2 );
    fgaus2->SetParameter(1, mean2);
    fgaus2->SetParameter(2, sigma2);
    histo2->Fit("fgaus2", "R");
    cout <<"fgaus2 : "<< "Chi^2: " << fgaus2->GetChisquare() << ", number DoF: " << fgaus2->GetNDF() << ", (Probability: " << fgaus2->GetProb() << ")." << endl;

    histo2->GetXaxis()->SetTitle(" Zrec - Zsim [cm] ");

    
//-------------------------------------------------------------------------------------------

    TFile MyFile("residui.root","RECREATE");
    histo->Write();
    histo2->Write();

}
//-----------------------------------------------------------------------------------------------------

bool CheckMult(int mult){
    double min=63.5; //impostare limiti molteplicità 
    double max=64.5; //impostare limiti molteplicità
    if( mult> min && mult< max) return true;
    else return false;
}


