
void extract_prj(string inFile, string saveStr, bool saveRoot = false){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* inF = new TFile(inFile.c_str());
  TFile* outF = nullptr;
  if(saveRoot) outF = new TFile((saveStr+".root").c_str(),"recreate");

  bool DoSyst = false;

  TH1D* sys_up_h;
  TH1D* sys_dw_h;

  TH1D* nominal_h = (TH1D*)inF->Get("nominal_h");
  if(DoSyst){
   sys_up_h = (TH1D*)inF->Get("sys_envelope_up_h");
   sys_dw_h = (TH1D*)inF->Get("sys_envelope_dw_h");
  }

  TLatex* lA = new TLatex(0.58,0.86,"ATLAS Internal");
  lA->SetNDC(kTRUE); lA->SetTextSize(0.04); lA->SetTextFont(72);

  TLatex* lSq2;
  if (saveStr.find("WP70") != std::string::npos) {
   lSq2 = new TLatex(0.58,0.78, "MV2c10 FixedCut 70%WP");
  }else if (saveStr.find("WP77") != std::string::npos) {
   lSq2 = new TLatex(0.58,0.78, "MV2c10 FixedCut 77%WP");
  }else if (saveStr.find("WP85") != std::string::npos) {
   lSq2 = new TLatex(0.58,0.78, "MV2c10 FixedCut 85%WP");
  }

  TLatex* lSq3;
  if (saveStr.find("PF") != std::string::npos) {
   lSq3 = new TLatex(0.58,0.74, "AntiKt4EMPFlowJet");
  }else if (saveStr.find("Topo") != std::string::npos) {
   lSq3 = new TLatex(0.58,0.74, "AntiKt4EMTopoJet");
  }
  //Would need to add condition to them if we use mc16d & e in the future
  TLatex* lSq = new TLatex(0.58,0.82,"#sqrt{s} = 13 TeV, 36.1 fb^{-1}");
  TLatex* lSq4 = new TLatex(0.58,0.70,"FTAG1 t#bar{t} sample, mc16a");
  lSq->SetNDC(kTRUE); lSq->SetTextSize(0.03); lSq->SetTextFont(42);
  lSq2->SetNDC(kTRUE); lSq2->SetTextSize(0.03); lSq2->SetTextFont(42);
  lSq3->SetNDC(kTRUE); lSq3->SetTextSize(0.03); lSq3->SetTextFont(42);
  lSq4->SetNDC(kTRUE); lSq4->SetTextSize(0.03); lSq4->SetTextFont(42);

  TCanvas* cen_cv = new TCanvas("cen_cv", "cen_cv",800,700);
  cen_cv->SetLogx();
  cen_cv->SetGrid();

  nominal_h->SetMinimum(0);
  nominal_h->SetMaximum(3);
  nominal_h->GetYaxis()->SetTitle("LF-jet SF");
  nominal_h->SetMarkerStyle(2);
  nominal_h->SetLineColor(kBlack);
  nominal_h->GetYaxis()->SetRangeUser(0.8,1.5);
  nominal_h->Draw("P");
  if (DoSyst){
    sys_up_h->SetFillColor(kAzure-4);
    sys_up_h->SetLineColor(kWhite);
    sys_up_h->SetFillStyle(3144);
    sys_dw_h->SetLineColor(kWhite);
    sys_dw_h->SetFillColor(kWhite-1);
    sys_dw_h->SetFillStyle(3004);
    sys_up_h->Draw("same hist");
    sys_dw_h->Draw("same hist");
    nominal_h->Draw("same P");
  }
  gPad->RedrawAxis();

  TLegend* leg =new TLegend(0.12, 0.86, 0.4, 0.75); leg->SetFillStyle(0);leg->SetFillColor(0);leg->SetBorderSize(0);
  leg->AddEntry(nominal_h, "nominal SF","lp");
  leg->Draw();
  if (DoSyst){
    leg->AddEntry(sys_up_h, "Syst. uncertainty","f");
  }

  lA->Draw(); lSq->Draw();lSq2->Draw();lSq3->Draw();lSq4->Draw();

  cen_cv->Print((saveStr+".pdf").c_str());
  cen_cv->Print((saveStr+".eps").c_str());
  cen_cv->Print((saveStr+".png").c_str());


  if(saveRoot){
    outF->cd();
    nominal_h->Write();
    cen_cv->Write();
    if (DoSyst){
    sys_up_h->Write();
    sys_dw_h->Write();
    }
  }

}
