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
#include<TEfficiency.h>
#include<TLegend.h>
#include<TMultiGraph.h>
#include<TCanvas.h>
#include<vector>

void readGraph(){

    TFile hfile("effDistro.root");
    TEfficiency* eff1= (TEfficiency*) hfile.Get("Effvsmul");
    eff1->SetMarkerStyle(20);
    eff1->SetName("Mult Distro");
    
    TFile infile("effMultUni.root");
    TEfficiency* eff2= (TEfficiency*) infile.Get("Effvsmul");
    eff2->SetName("Mult uniforme");
    eff2->SetTitle("Efficienza vs Molteplicita'");
    eff2->SetMarkerStyle(20);
    eff2->SetMarkerColor(2);
    eff2->SetFillStyle(2);
    eff2->SetLineColor(2);

    TCanvas* c= new TCanvas("c","Efficienza vs Molteplicità ");
    c->SetName(" Efficienza vs Molteplicità");
    eff2->Draw();
    eff1->Draw("same");
    TLegend* leg1 = new TLegend(0.3,0.1,0.7,0.5);
    leg1->AddEntry(eff1,"distr. assegnata");
    leg1->AddEntry(eff2,"distr. uniforme");
    leg1->Draw("same");

}