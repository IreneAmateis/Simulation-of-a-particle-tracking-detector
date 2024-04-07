//macro per sovrapporre i grafici di risoluzione vs molteplicità e confrontare tra molteplicità uniforme e assegnata
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
#include<TMultiGraph.h>
#include<TCanvas.h>
#include<vector>

void MyGraph(){
    TFile hfile("analisiDistro.root");
    TGraphErrors* Myplot1r= (TGraphErrors*) hfile.Get("resvsmul");

    Myplot1r->SetFillStyle(1);
    Myplot1r->SetMarkerColor(1);
    Myplot1r->SetMarkerStyle(20);
    Myplot1r->SetLineColor(1);
    
    TFile infile("analisiMultUni.root");
    TGraphErrors* Myplot2r= (TGraphErrors*) infile.Get("resvsmul");
    Myplot2r->SetFillStyle(2);
    Myplot2r->SetMarkerColor(2);
    Myplot2r->SetLineColor(2);

    Myplot2r->SetMarkerStyle(20);
    
    TCanvas* d= new TCanvas();
    TMultiGraph* mg2= new TMultiGraph();
    mg2->Add(Myplot1r);
    mg2->Add(Myplot2r);
    mg2->SetTitle("Risoluzione vs Molteplicita'");
    mg2->GetXaxis()->SetTitle("Molteplicita' ");
    mg2->GetYaxis()->SetTitle(" Risoluzione [#mum]");
    mg2->Draw("ap");

    TLegend* leg1 = new TLegend(0.3,0.1,0.7,0.5);
    leg1->AddEntry(Myplot1r,"distr. assegnata");
    leg1->AddEntry(Myplot2r,"distr. uniforme");
    leg1->Draw("same");

    TFile ofile("MyGraph.root","recreate");
    mg2->Write("resvsmul");

}