#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>

using namespace std;

int doSyst = 0;//draw sys breakdown if doSys>1
bool removeFakes = false; 
bool SmearingOnly = true; //Only use D0 and Z0, testing in Nominal only
bool useJETType = false; //Used if jets are flagged as hadInt/strange etc, in place that dropping tracks
int correlationSys = 1; //0 == no corrSys, 1 == D0Z0 corrSys, 2 == D0Z0Corl corrSys


void setAtOne(TH1* h1){
  for( int bx = 1; bx < h1->GetNbinsX()+1; ++bx ) {
    for(int by = 1; by < h1->GetNbinsY()+1; ++by) {
      h1->SetBinContent(bx, by, 1);
      h1->SetBinError(bx, by, 0);
    }
  }
}

void deltaOne_QuadAdd(TH1* h1, TH1* add_h){
  for( int bx = 1; bx < h1->GetNbinsX()+1; ++bx ) {
    for(int by = 1; by < h1->GetNbinsY()+1; ++by) {
      double newBinValue = 1 + sqrt( pow( (h1->GetBinContent(bx, by) -1) ,2) + 
				     pow( (add_h->GetBinContent(bx, by) -1), 2) );
      h1->SetBinContent(bx, by, newBinValue);
      h1->SetBinError(bx, by, h1->GetBinError(bx, by) ); //++++++++TO BE IMPROVED!!!
    }
  }
}

void deltaOne_Add(TH1* h1, TH1* add_h){
  for( int bx = 1; bx < h1->GetNbinsX()+1; ++bx ) {
    for(int by = 1; by < h1->GetNbinsY()+1; ++by) {
      double newBinValue = 1 + (h1->GetBinContent(bx, by) -1) + (add_h->GetBinContent(bx, by) -1) ;
      h1->SetBinContent(bx, by, newBinValue);
      h1->SetBinError(bx, by, sqrt(pow(h1->GetBinError(bx, by),2) + pow(add_h->GetBinError(bx, by),2) ) );
    }
  }
}

void multiply(TH1* h1, TH1* multiply_h){
  for( int bx = 1; bx < h1->GetNbinsX()+1; ++bx ) {
    for(int by = 1; by < h1->GetNbinsY()+1; ++by) {
      double newBinValue = h1->GetBinContent(bx, by) * multiply_h->GetBinContent(bx, by) ;
      h1->SetBinContent(bx, by, newBinValue);
      h1->SetBinError(bx, by, sqrt(pow(multiply_h->GetBinContent(bx, by)*h1->GetBinError(bx, by),2) + pow(h1->GetBinContent(bx, by)*multiply_h->GetBinError(bx, by),2) ) );
    }
  }
}

TH1* symmetryze_wrt_one(TH1* h1){
  TH1* symm_h = (TH1*) h1->Clone( Form("symm_%s",h1->GetName()));
  symm_h->Reset();
  for(int bx = 1; bx < h1->GetNbinsX()+1; ++bx ){
      for(int by = 1; by < h1->GetNbinsY()+1; ++by) {
        double newBinValue = 2 - h1->GetBinContent(bx, by);
        symm_h->SetBinContent(bx, by, newBinValue);
      }
  }
  return symm_h;
}

TH1* symmetryze_wrt_h(TH1* h1, TH1* h_baseline){
  TH1* symm_h = (TH1*) h1->Clone( Form("symm_%s",h1->GetName()));
  symm_h->Reset();

  for(int bx = 1; bx < h1->GetNbinsX()+1; ++bx ){
      for(int by = 1; by < h1->GetNbinsY()+1; ++by) {
        double baselineBinValue = h_baseline->GetBinContent(bx, by);
        double delta = h1->GetBinContent(bx, by)  - baselineBinValue ;
        //symmetrize w.r.t. baseline
        if(delta>0) symm_h->SetBinContent(bx, by, baselineBinValue - fabs(delta));
        else  symm_h->SetBinContent(bx, by, baselineBinValue + fabs(delta));
      }
  }

  return symm_h;
}



TH1* combine_hists(TFile* inF, string histBaseName, vector<string> multiply_v, 
                   TFile* inF2 = nullptr, vector<string> multiplyF2_v = vector<string>() ){
  string startName = multiply_v[0];

  TH1* ret_h = (TH1*)(inF->Get( (histBaseName+"_"+startName).c_str()))->Clone("ret_h");
  ret_h->Reset(); 
  setAtOne(ret_h);

  for(uint mult_i = 0; mult_i< multiply_v.size(); ++mult_i){
    string effectStr = multiply_v[mult_i];
    cout<<"Multiply values of hist: "<<effectStr<<endl;
    TH1* tmp_h = (TH1*) inF->Get( (histBaseName+"_"+effectStr).c_str());
    if(!tmp_h) cerr<<"ERROR: can not retrieve hist:" <<histBaseName+"_"+effectStr<<endl;

    if(effectStr.find("FAKE_RATE") != string::npos || 
       effectStr.find("TAIL") != string::npos || 
       effectStr.find("STRANGEHAD") != string::npos || 
       effectStr.find("CONV_HADINT") != string::npos ) {
      cout<<"Symmetrization around 1 of hist: "<<effectStr<<endl;
      tmp_h = symmetryze_wrt_one(tmp_h);
    }

    multiply(ret_h, tmp_h);
  }

  if(inF2){
    for(uint mult_i = 0; mult_i< multiplyF2_v.size(); ++mult_i){
      string effectStr = multiplyF2_v[mult_i];
      cout<<"File2, multiply values of hist: "<<effectStr<<endl;
      TH1* tmp_h = (TH1*) inF2->Get( (histBaseName+"_"+effectStr).c_str());
      if(!tmp_h) cerr<<"ERROR: can not retrieve hist:" <<histBaseName+"_"+effectStr<<endl;

      if(effectStr.find("FAKE_RATE") != string::npos || 
       effectStr.find("TAIL") != string::npos || 
       effectStr.find("STRANGEHAD") != string::npos || 
       effectStr.find("CONV_HADINT") != string::npos ) {
        cout<<"Symmetrization around 1 of hist: "<<effectStr<<endl;
        tmp_h = symmetryze_wrt_one(tmp_h);
      }
      
      multiply(ret_h, tmp_h);
    }
  }

  return ret_h;
}


void combine_effects(string inFileStr, string outName, string histName, string CDItxt = "", string plotStr = "", string generatorSysFile = "", string nominalFile = "", bool doIPSymm = false){
  TH1::SetDefaultSumw2(kTRUE);
  TH1::AddDirectory(kFALSE);

  //This file contains all syst and nominal
  TFile f_in(inFileStr.c_str());
  TFile* f_nominalOnly = nullptr;
  if(nominalFile != "") f_nominalOnly = new TFile(nominalFile.c_str());

  string STRANGE_EFF = useJETType ? "TYPE_SHadrons_var" : "STRANGEHAD_RATE_var";
  string CONV_HADINT_EFF = useJETType ? "TYPE_ConvOrHadInt_var" : "CONV_HADINT_RATE_var";

  string resStr = "RES"; //it was "RES_TIGHT"
  vector<string>  multiply_v = vector<string>(); //multiplication is the right way to combine SFs!
  multiply_v = {resStr+"_D0_MEAS_var",resStr+"_Z0_MEAS_var"};
  if (!SmearingOnly){
    multiply_v.push_back(STRANGE_EFF);
    multiply_v.push_back("TAIL_RATE_var");
  };
  if(!removeFakes)  multiply_v.push_back("FAKE_RATE_LOOSE_var");

  TH1* nominal_h = f_nominalOnly ? //Load from f_nominalOnly if I gave it before
    combine_hists(f_nominalOnly, histName, multiply_v ):
    combine_hists(&f_in, histName, multiply_v );
  nominal_h->SetNameTitle("nominal_h","nominal_h");


  /////////////  Systematic effects ////////////////
  vector<TH1*> sys_h_v;
  if(doSyst){    
    if(generatorSysFile!= ""){
      //Generator comparison and sytematics, done with Herwig
      cout<<"Generator Syst."<<endl;
      TFile generator_sys_f(generatorSysFile.c_str()); 
      multiply_v = {resStr+"_D0_MEAS_var",resStr+"_Z0_MEAS_var","TAIL_RATE_var"};
      if(!removeFakes)  multiply_v.push_back("FAKE_RATE_LOOSE_var");
      vector<string> multiplyF2_v = {STRANGE_EFF};//Take this from nominal file
      TH1* generator_sys_h = combine_hists(&generator_sys_f, histName, multiply_v,
                                           &f_in, multiplyF2_v);
      generator_sys_h->SetNameTitle("FT_EFF_generator_sys_UP","FT_EFF_generator_sys_UP");
      generator_sys_h->SetLineColor(kAzure+1);
      TH1* generator_sys_sym_h = symmetryze_wrt_h(generator_sys_h, nominal_h);
      generator_sys_sym_h->SetNameTitle("FT_EFF_generator_sys_DOWN","FT_EFF_generator_sys_DOWN");
      sys_h_v.push_back(generator_sys_h); 
      sys_h_v.push_back(generator_sys_sym_h); 
    }

    cout<<"Smearing Syst."<<endl;
    multiply_v = {resStr+"_D0_MEAS_UP_var", resStr+"_Z0_MEAS_var", STRANGE_EFF, "TAIL_RATE_var"};
    if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
    TH1* d0_smearing_UP_h = combine_hists(&f_in, histName, multiply_v );
    d0_smearing_UP_h->SetNameTitle("FT_EFF_d0_smearing_UP","FT_EFF_d0_smearing_UP");
    d0_smearing_UP_h->SetLineColor(8);
    sys_h_v.push_back(d0_smearing_UP_h);

    multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_UP_var", STRANGE_EFF, "TAIL_RATE_var"};
    if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
    TH1* z0_smearing_UP_h = combine_hists(&f_in, histName, multiply_v );
    z0_smearing_UP_h->SetNameTitle("FT_EFF_z0_smearing_UP","FT_EFF_z0_smearing_UP");
    z0_smearing_UP_h->SetLineColor(kGreen+3);
    sys_h_v.push_back(z0_smearing_UP_h);

    if(doIPSymm){
      TH1* d0_smearing_sym_h = symmetryze_wrt_h(d0_smearing_UP_h, nominal_h);
      d0_smearing_sym_h->SetNameTitle("FT_EFF_d0_smearing_DOWN"," FT_EFF_d0_smearing_DOWN");
      sys_h_v.push_back(d0_smearing_sym_h); 

      TH1* z0_smearing_sym_h = symmetryze_wrt_h(z0_smearing_UP_h, nominal_h);
      z0_smearing_sym_h->SetNameTitle("FT_EFF_z0_smearing_DOWN"," FT_EFF_z0_smearing_DOWN");
      sys_h_v.push_back(z0_smearing_sym_h); 
    } else {
      multiply_v = {resStr+"_D0_MEAS_DOWN_var", resStr+"_Z0_MEAS_var", STRANGE_EFF, "TAIL_RATE_var"};
      if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
      TH1* d0_smearing_DOWN_h = combine_hists(&f_in, histName, multiply_v );
      d0_smearing_DOWN_h->SetNameTitle("FT_EFF_d0_smearing_DOWN","FT_EFF_d0_smearing_DOWN");
      d0_smearing_DOWN_h->SetLineColor(8);
      sys_h_v.push_back(d0_smearing_DOWN_h);
      
      multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_DOWN_var", STRANGE_EFF, "TAIL_RATE_var"};
      if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
      TH1* z0_smearing_DOWN_h = combine_hists(&f_in, histName, multiply_v );
      z0_smearing_DOWN_h->SetNameTitle("FT_EFF_z0_smearing_DOWN","FT_EFF_z0_smearing_DOWN");
      z0_smearing_DOWN_h->SetLineColor(kGreen+3);
      sys_h_v.push_back(z0_smearing_DOWN_h);
      
      multiply_v = {resStr+"_D0_MEAS_Zmm_var", resStr+"_Z0_MEAS_var", STRANGE_EFF, "TAIL_RATE_var"};
      if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
      TH1* d0_smearing_Zmm_h = combine_hists(&f_in, histName, multiply_v );
      d0_smearing_Zmm_h->SetNameTitle("FT_EFF_d0_smearing_Zmm","FT_EFF_d0_smearing_Zmm");
      d0_smearing_Zmm_h->SetLineColor(kCyan+1);
      sys_h_v.push_back(d0_smearing_Zmm_h);

      if(CDItxt != ""){//CDI needs an UP variation to be symmetrized, here I produce it
        TH1* d0_smearing_Zmm_sym_h = symmetryze_wrt_h(d0_smearing_Zmm_h, nominal_h);
        d0_smearing_Zmm_sym_h->SetNameTitle("FT_EFF_d0_smearing_Zmm_UP"," FT_EFF_d0_smearing_Zmm_UP");
        sys_h_v.push_back(d0_smearing_Zmm_sym_h); 
      }

      multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_Zmm_var", STRANGE_EFF, "TAIL_RATE_var"};
      if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
      TH1* z0_smearing_Zmm_h = combine_hists(&f_in, histName, multiply_v );
      z0_smearing_Zmm_h->SetNameTitle("FT_EFF_z0_smearing_Zmm","FT_EFF_z0_smearing_Zmm");
      z0_smearing_Zmm_h->SetLineColor(kBlue);
      sys_h_v.push_back(z0_smearing_Zmm_h);

      if(CDItxt != ""){//CDI needs an UP variation to be symmetrized, here I produce it
        TH1* z0_smearing_Zmm_sym_h = symmetryze_wrt_h(z0_smearing_Zmm_h, nominal_h);
        z0_smearing_Zmm_sym_h->SetNameTitle("FT_EFF_z0_smearing_Zmm_UP"," FT_EFF_z0_smearing_Zmm_UP");
        sys_h_v.push_back(z0_smearing_Zmm_sym_h); 
      }
    }

    if(!removeFakes){//Use tight fakes for syst.
      cout<<"Fake Syst."<<endl;
      multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_var", STRANGE_EFF, "TAIL_RATE_var"};
      multiply_v.push_back("FAKE_RATE_TIGHT_var");
      TH1* fake_h = combine_hists(&f_in, histName, multiply_v );
      fake_h->SetNameTitle("FT_EFF_fake_UP","FT_EFF_fake_UP");
      fake_h->SetLineColor(94);
      TH1* fake_sym_h = symmetryze_wrt_h(fake_h, nominal_h);
      fake_sym_h->SetNameTitle("FT_EFF_fake_DOWN","FT_EFF_fake_DOWN");
      sys_h_v.push_back(fake_h); 
      sys_h_v.push_back(fake_sym_h); 
    }


    cout<<"Strange-hadrons Syst."<<endl;
    multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_var", "TAIL_RATE_var"};
    if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");    
    TH1* sHad_h = combine_hists(&f_in, histName, multiply_v );
    sHad_h->SetNameTitle("FT_EFF_sHad_UP","FT_EFF_sHad_UP");
    sHad_h->SetLineColor(kYellow+1);
    TH1* sHad_sym_h = symmetryze_wrt_h(sHad_h, nominal_h);
    sHad_sym_h->SetNameTitle("FT_EFF_sHad_DOWN","FT_EFF_sHad_DOWN");    
    sys_h_v.push_back(sHad_h); 
    sys_h_v.push_back(sHad_sym_h); 

    cout<<"Material interaction Syst."<<endl;
    multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_var", STRANGE_EFF, "TAIL_RATE_var", CONV_HADINT_EFF};
    if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");    
    TH1* matInt_h = combine_hists(&f_in, histName, multiply_v );
    matInt_h->SetNameTitle("FT_EFF_matInt_UP","FT_EFF_matInt_UP");
    matInt_h->SetLineColor(kPink+1);
    TH1* matInt_sym_h = symmetryze_wrt_h(matInt_h, nominal_h);
    matInt_sym_h->SetNameTitle("FT_EFF_matInt_DOWN","FT_EFF_matInt_DOWN");   
    sys_h_v.push_back(matInt_h); 
    sys_h_v.push_back(matInt_sym_h); 

    cout<<"TAIL RATE Syst."<<endl;
    multiply_v = {resStr+"_D0_MEAS_var", resStr+"_Z0_MEAS_var", STRANGE_EFF,  "TAIL_RATE_SYS_var"};
    if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
    TH1* tailRateSys_h = combine_hists(&f_in, histName, multiply_v );
    tailRateSys_h->SetNameTitle("FT_EFF_tailRateSys_UP","FT_EFF_tailRateSys_UP");
    tailRateSys_h->SetLineColor(kRed);
    TH1* tailRateSys_sym_h = symmetryze_wrt_h(tailRateSys_h, nominal_h);
    tailRateSys_sym_h->SetNameTitle("FT_EFF_tailRateSys_DOWN","FT_EFF_tailRateSys_DOWN");   
    sys_h_v.push_back(tailRateSys_h); 
    sys_h_v.push_back(tailRateSys_sym_h); 

    if(correlationSys == 1){
      cout<<"Correlation Syst."<<endl;
      multiply_v = {resStr+"_D0Z0_MEAS_var", STRANGE_EFF, "TAIL_RATE_var"};
      if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
      TH1* correl1_h = combine_hists(&f_in, histName, multiply_v );
      correl1_h->SetNameTitle("FT_EFF_correl_UP","FT_EFF_correl_UP");
      correl1_h->SetLineColor(100);
      TH1* correl1_sym_h = symmetryze_wrt_h(correl1_h, nominal_h);
      correl1_sym_h->SetNameTitle("FT_EFF_correl_DOWN","FT_EFF_correl_DOWN");
      sys_h_v.push_back(correl1_h);
      sys_h_v.push_back(correl1_sym_h);
    }
    
    if(correlationSys == 2){
      cout<<"Correlation Syst. with simultaneous smearing"<<endl;
      multiply_v = {resStr+"_D0Z0Corl_MEAS_var", STRANGE_EFF,  "TAIL_RATE_var"};
      if(!removeFakes) multiply_v.push_back("FAKE_RATE_LOOSE_var");
      TH1* correl_h = combine_hists(&f_in, histName, multiply_v );
      correl_h->SetNameTitle("FT_EFF_corlSmr_UP","FT_EFF_corlSmr_UP");
      correl_h->SetLineColor(kMagenta);
      TH1* correl_sym_h = symmetryze_wrt_h(correl_h, nominal_h);
      correl_sym_h->SetNameTitle("FT_EFF_corlSmr_DOWN","FT_EFF_corlSmr_DOWN");
      sys_h_v.push_back(correl_h);
      sys_h_v.push_back(correl_sym_h);
    }
  }

  //Save histos
  TFile f_out((outName+".root").c_str(), "recreate");
  nominal_h->Write();
  TH1* sys_envelope_up_h =nullptr, *sys_envelope_dw_h = nullptr;
  if(doSyst) {
    //Prepare and save the up and down envelope variations
    sys_envelope_up_h = (TH1*)nominal_h->Clone();
    sys_envelope_up_h->SetNameTitle("sys_envelope_up_h","sys_envelope_up_h");
    sys_envelope_dw_h = (TH1*)nominal_h->Clone();
    sys_envelope_dw_h->SetNameTitle("sys_envelope_dw_h","sys_envelope_dw_h");

    for(uint i=0; i<sys_h_v.size(); ++i) sys_h_v[i]->Write();
   
    for(int bx = 1; bx < nominal_h->GetNbinsX()+1; ++bx ){
      for(int by = 1; by < nominal_h->GetNbinsY()+1; ++by) {
        double bin_val = nominal_h->GetBinContent(bx, by);
        double binSqrtAdd_up = 0;
        double binSqrtAdd_dw = 0;

        for(uint i=0; i<sys_h_v.size(); ++i) {
          string sysName = sys_h_v[i]->GetName();
          
          double sys_bin_val = sys_h_v[i]->GetBinContent(bx, by);
          if(sys_bin_val > bin_val) binSqrtAdd_up += pow((sys_bin_val - bin_val), 2);
          else binSqrtAdd_dw += pow((sys_bin_val - bin_val), 2);

        }
        //Add MC error from nominal to the envelope
        double bin_err_sq = pow(nominal_h->GetBinError(bx, by),2); 
        sys_envelope_up_h->SetBinContent(bx, by, bin_val+sqrt(binSqrtAdd_up+bin_err_sq));
        sys_envelope_dw_h->SetBinContent(bx, by, bin_val-sqrt(binSqrtAdd_dw+bin_err_sq));
        //cout<<"Nominal = "<<bin_val<<" , MCStatErr = "<<nominal_h->GetBinError(bx, by)<<" , SysErr = "<<sqrt(binSqrtAdd_up)<<endl;
      }
    }

    sys_envelope_up_h->Write();
    sys_envelope_dw_h->Write();  
  }

  if(CDItxt != ""){

    FILE * CDI_f;
    char CDIName[100];
    CDI_f = fopen (CDItxt.c_str(),"w");

    string WP ="", calibPoint = "";

    if(outName.find("50") != string::npos){WP = "50"; calibPoint = "0.976933";}
    if(outName.find("60") != string::npos){WP = "60"; calibPoint = "0.934906";}
    if(outName.find("70") != string::npos){WP = "70"; calibPoint = "0.824427";}
    if(outName.find("77") != string::npos){WP = "77"; calibPoint = "0.645925";}
    if(outName.find("85") != string::npos){WP = "85"; calibPoint = "0.175847";}

    //header
    fprintf(CDI_f,"Analysis(AdjustedMC,light,MV2c10,FixedCutBEff_%s,AntiKt4EMTopoJets){\n",WP.c_str());

    fprintf(CDI_f,"meta_data_s (Hadronization, Pythia8EvtGen)\n");
    fprintf(CDI_f,"meta_data_s (OperatingPoint, %s)\n", calibPoint.c_str());

    TAxis* x_ax =nominal_h->GetXaxis();
    TAxis* y_ax =nominal_h->GetYaxis();

    for(int bx = 1; bx < nominal_h->GetNbinsX()+1; ++bx ) {
      for(int by = 1; by < nominal_h->GetNbinsY()+1; ++by) {

        fprintf(CDI_f,"bin(%.1f<pt<%.1f,%.2f<abseta<%.2f)\n", 
                x_ax->GetBinLowEdge(bx), x_ax->GetBinUpEdge(bx),
                y_ax->GetBinLowEdge(by), y_ax->GetBinUpEdge(by));
        fprintf(CDI_f,"{\n");
        fprintf(CDI_f,"    central_value(%.3f,%.3f)\n", 
                nominal_h->GetBinContent(bx, by), nominal_h->GetBinError(bx, by));
        //++++++        fprintf("    meta_data(N_jets tagged,5780.8,131.8)");
        //++++++        fprintf("    meta_data(N_jets total,11379.6,292.4)");
        //+++++    fprintf("    usys(FT_EFF_MC statistics,1.5%)");
        
        //Doing this only for the UP variation because the CDI file uses symmetric sys
        for(uint i=0; i<sys_h_v.size(); ++i){
          string sysName = sys_h_v[i]->GetName();
          size_t position = sysName.find("_UP");
          if( position != string::npos) {
            sysName.erase(position, position+3);
            double sysVar = 100.*(sys_h_v[i]->GetBinContent(bx, by) - nominal_h->GetBinContent(bx,by))/nominal_h->GetBinContent(bx,by);
            fprintf(CDI_f,"    sys(%s,%.2f%%)\n", sysName.c_str(), sysVar);
          }
        }
        fprintf(CDI_f,"}\n");
      }
    }
    fprintf(CDI_f,"}\n");

    fclose (CDI_f);
  }

  if(plotStr != ""){
    TCanvas* plot_cv = new TCanvas("plot_cv", "plot_cv",800,700);
    if(histName.find("pt")!=string::npos) plot_cv->SetLogx();
    plot_cv->SetRightMargin(0.07);
    plot_cv->SetTopMargin(0.07);
    plot_cv->SetBottomMargin(0.12);
    plot_cv->SetLeftMargin(0.11);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadTickY(1);
    gStyle->SetPadTickX(1); 
    nominal_h->GetXaxis()->SetTitleOffset(1.3);
    nominal_h->GetYaxis()->SetTitleOffset(1.3);
    nominal_h->SetLineWidth(2);
    nominal_h->SetMinimum(0.55);
    nominal_h->SetMaximum(3.9);
    if( plotStr.find("_C") != string::npos ) nominal_h->SetMaximum(1.6);
    nominal_h->GetYaxis()->SetTitle("LF-jet SF");

    if(nominal_h->GetNbinsY()>1) {     
      nominal_h->GetYaxis()->SetTitle("jet |#eta|");
      gStyle->SetPalette(57);
      gStyle->SetNumberContours(60);
      gStyle->SetPaintTextFormat("2.2f");
      
      nominal_h->SetMaximum(1.5);
      nominal_h->SetMarkerSize(1.5);
      
      nominal_h->Draw("COLZTEXT");
      TLatex* lab = new TLatex(0.88,0.92,"LF-SF");
      lab->SetNDC(kTRUE); lab->SetTextSize(0.04); lab->SetTextFont(72);
      lab->Draw();  
    } else {
      nominal_h->Draw();
    }

    TLatex* tl = new TLatex(0.63,0.87,plotStr.c_str());
    tl->SetNDC(kTRUE);  tl->SetTextSize(0.05);  tl->SetTextFont(42);
    TLatex* lA = new TLatex(0.63,0.79,"ATLAS Internal");
    lA->SetNDC(kTRUE); lA->SetTextSize(0.04); lA->SetTextFont(72);
    TLatex* lSq = new TLatex(0.63,0.735,"#sqrt{s} = 13 TeV");
    lSq->SetNDC(kTRUE); lSq->SetTextSize(0.04); lSq->SetTextFont(42);
    tl->Draw(); lA->Draw(); lSq->Draw();

    if(doSyst){
      TLegend* legSys =new TLegend(0.18, 0.14, 0.58, 0.38); 
      legSys->SetFillStyle(0);legSys->SetFillColor(0);legSys->SetBorderSize(0);
      sys_envelope_up_h->SetLineWidth(2);
      sys_envelope_dw_h->SetLineWidth(2);
      sys_envelope_up_h->Draw("samehist");
      sys_envelope_dw_h->Draw("samehist"); 

      legSys->AddEntry(sys_envelope_up_h, "TOTAL Syst. (envelope)","lp");

      if(doSyst>1){//draw the sys breakdown
        for(uint i=0; i<sys_h_v.size(); ++i){
          sys_h_v[i]->SetLineWidth(2);
          sys_h_v[i]->SetLineStyle(2);
          sys_h_v[i]->Draw("samehist");
          string sysLeg = sys_h_v[i]->GetTitle();
          size_t position_up = sysLeg.find("_UP");
          size_t position_dw = sysLeg.find("_DOWN");
          if( position_up != string::npos) {
            sysLeg.erase(position_up, position_up+3);
            legSys->AddEntry(sys_h_v[i], sysLeg.c_str(),"l");
          } else if( position_dw != string::npos) {
            continue; //Do not add up and down separately in the legend
          } else {
            legSys->AddEntry(sys_h_v[i], sysLeg.c_str(),"l");
          }
        }
      }
      nominal_h->Draw("axis same");
      legSys->Draw();
    }
    plot_cv->Print((plotStr+".eps").c_str());
    //    plot_cv->Print((plotStr+".pdf").c_str());
  }
  f_out.Close();
}
