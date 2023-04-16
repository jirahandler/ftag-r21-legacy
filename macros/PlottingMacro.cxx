#include <fstream>
#include <stdint.h>
#include <TLatex.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>
#include <cmath>
#include <string>
#include <vector>
#include <stdint.h>
#include "assert.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"

#include "atlasstyle-00-03-05/AtlasLabels.C"

using namespace std;

void PlottingMacro(){

  
   //define canvas and legend
   TCanvas *C = new TCanvas("C","",0.,0.,800,600);
   C->Divide(1,1);

   TLegend *leg = new TLegend(0.45,0.72,0.60,0.87); 
   leg->SetTextFont(42);
   leg->SetTextSize(0.04);
   leg->SetLineColor(0);
   leg->SetLineStyle(1);
   leg->SetLineWidth(12);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   

  
   //Light Flavor Studies
   std::cout << "VHF Light Flavor Studies" << std::endl;

   TFile ttbarExample("/afs/cern.ch/user/v/vcairo/Work/bTagging/xAODAthena/run/flav_Akt4EMTo_withMV2C20Flip.root");
   
   TTree *tree;
   
   ttbarExample.GetObject("bTag_AntiKt4EMTopoJets",tree);
   
   TH1F *mv2c20_LF	= new TH1F("mv2c20_LF", "mv2c20_LF", 100, -2, 2.);
   TH1F *mv2c20_B	= new TH1F("mv2c20_B",  "mv2c20_B",  100, -2, 2.);
   
   tree->Draw("jet_mv2c20>>mv2c20_LF","jet_truthflav==0");
   mv2c20_LF->GetXaxis()->SetTitle("mv2c20");
   mv2c20_LF->GetYaxis()->SetTitle("Entries");
   mv2c20_LF->SetTitle("");
   mv2c20_LF->GetYaxis()->SetTitleOffset(1.5);
   mv2c20_LF->SetStats(0);
   mv2c20_LF->SetLineColor(1);
   mv2c20_LF->SetLineWidth(2);
   //xProfile->SetMaximum(10000);
   leg->AddEntry(mv2c20_LF,"Light Flavor jets", "l");

   tree->Draw("jet_mv2c20>>mv2c20_B","jet_truthflav==5");
   mv2c20_B->GetXaxis()->SetTitle("mv2c20");
   mv2c20_B->GetYaxis()->SetTitle("Entries");
   mv2c20_B->SetLineColor(2);
   mv2c20_B->SetLineWidth(2);
   leg->AddEntry(mv2c20_B,"B jets", "l");

   
   C->cd(1);
   mv2c20_LF ->Draw("hist");
   mv2c20_B  ->Draw("hist,same");
   ATLASLabelSqrt(0.62,0.49, "Internal");
   leg->Draw();
   C->Print("Profiles.pdf");
   C->Clear();  
}
