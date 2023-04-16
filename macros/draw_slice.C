
void draw_slice(){

  TFile f("trkSys_SF_85wp.root");
  gStyle->SetOptStat(0);

  TCanvas cent_cv("cent_cv","cent_cv",600,600);
  cent_cv.SetLogx();
  TH1D* cent_h = ((TH2F*) f.Get("nominal_h")) ->ProjectionX("cent",1,1);
  cent_h->GetYaxis()->SetTitle("LF  SF");
  cent_h->SetMinimum(0.5); cent_h->SetMaximum(2.5);
  cent_h->Draw("h");
  TLatex p_cent; 
  p_cent.SetNDC(); p_cent.SetTextFont(42); p_cent.SetTextSize(0.04); p_cent.SetTextColor(kBlack);
  p_cent.DrawLatex(0.20,0.23,Form("MV2c10 85%% working point, |#eta|<1.2"));
  cent_cv.Print("wp85_LF_SF_central_eta.pdf");


  TCanvas forw_cv("forw_cv","forw_cv",600,600);
  forw_cv.SetLogx();
  TH1D* forw_h =  ((TH2F*)f.Get("nominal_h"))->ProjectionX("forw",2,2);
  forw_h->GetYaxis()->SetTitle("LF  SF");
  forw_h->SetMinimum(0.5); forw_h->SetMaximum(2.5);
  forw_h->Draw("h");

  TLatex p_forw; 
  p_forw.SetNDC(); p_forw.SetTextFont(42); p_forw.SetTextSize(0.04); p_forw.SetTextColor(kBlack);
  p_forw.DrawLatex(0.20,0.23,Form("MV2c10 85%% working point, 1.2<|#eta|<2.5"));
  forw_cv.Print("wp85_LF_SF_forward_eta.pdf");
}
