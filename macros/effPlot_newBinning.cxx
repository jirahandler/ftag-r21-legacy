#include <fstream>
#include <stdint.h>
#include <TLatex.h>
#include <iostream>
#include <iomanip>
#include <TMath.h>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <stdint.h>
#include "assert.h"
#include "TString.h"
#include <unistd.h>
#include <dirent.h>
#include <cstdlib>
#include <time.h>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TChain.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"

#include "atlasstyle-00-03-05/AtlasLabels.C"

//set to 0 to speed up execution
#define _MORE_H 0
#define _DOTRK 0
#define _DOD0TailsJetSF 0 //This is the ad-hoc jet scaling study, now superseeded by the use of proper track variations

using namespace std;

double gD0TailSF = 1., gD0TailSF_UP = 1., gD0TailSF_DOWN = 1.;
//Scan files in a Directory
void ScanDir (string file1, TChain* myT_1, string sysName="") {
  TList* hList = new TList();
  cout << "Input is a directory: going fancy: " << endl;
  DIR * dir, * dir2;
  dirent * pdir, * pdir2;
  dir = opendir( file1.c_str() );     // open current directory
  while ((pdir = readdir(dir)))  {     // loop on directory contents
    string foldName=pdir->d_name;
    if (foldName == "..") continue;
    if (!sysName.empty() && sysName != foldName) continue;
    cout << pdir->d_name << endl;
    dir2 = opendir( (file1+"/"+foldName).c_str() );     // open current directory
    int nFiles = 0;
    while ((pdir2 = readdir(dir2)))  {
      string fName=pdir2->d_name;
      if (fName.find("root")==string::npos) continue;
      ++nFiles;
      //      cout << fName << endl;
      myT_1->Add( (file1+"/"+foldName+"/"+fName).c_str() );
    }
    cout <<"NFiles = "<< nFiles << endl;
  }
}

enum flavId{ LF_id = 0, C_id = 4, B_id =5 };


//container of the histograms used for calibration:
// Nominal, vector of variations, and vector of variation ratios, one for each flavour
class calib_h{

  string baseName;
  bool is1Dim;
  //maps of the flavors for the hists I am intested in
  map< int, TH1*> nominal_h_m;
  //map of flavor and sys names
  map< int, map< string, TH1*> > sys_flav_h_m;
  map< int, map< string, TH1*> > sys_flav_corl_h_m;//only correlated events
  map< int, map< string, TH1*> > sys_flav_uncorl_h_m;//only UN-correlated events
  map< int, map< string, TH1*> > sys_flav_variation_h_m;

public:
  calib_h(map<int, string> flav_m, vector<string> sys_v, TH1* exampleHist,
          string inputFile = ""){
    is1Dim = (exampleHist->GetNbinsY() == 1); //this is 2Dim?
    baseName = exampleHist->GetName();

    TFile* inputF = nullptr;
    inputFile=inputFile.substr(0,inputFile.find(".root"));
    if(inputFile != "") inputF = new TFile((inputFile+".root").c_str(), "read");

    for(auto flav_i : flav_m){
      int flavId = flav_i.first;
      string flavStr = flav_i.second;//should be LF, B, C, etc.
      string nameStr = baseName+"_"+flavStr;

      map<string, TH1*> sys_h_m;
      map<string, TH1*> sys_corl_h_m;
      map<string, TH1*> sys_uncorl_h_m;

      if(inputFile == ""){
        string labelX = exampleHist->GetXaxis()->GetTitle();
        string labelY = exampleHist->GetYaxis()->GetTitle();
        if(is1Dim) nominal_h_m[flavId] = (TH1D*)exampleHist->Clone(nameStr.c_str());
        else nominal_h_m[flavId] = (TH2D*)exampleHist->Clone(nameStr.c_str());
        nominal_h_m[flavId]->SetTitle(nameStr.c_str());
        nominal_h_m[flavId]->GetXaxis()->SetTitle( (labelX+" ("+flavStr+" jets)").c_str());
        if(is1Dim) {
          nominal_h_m[flavId]->GetYaxis()->SetTitle( "Entries");
        } else {
          nominal_h_m[flavId]->GetYaxis()->SetTitle( (labelY+" ("+flavStr+" jets)").c_str());
          nominal_h_m[flavId]->GetZaxis()->SetTitle( "Entries");
        }
      } else { 
        //Initialization from saved file. Inputs MUST be consistent
        if(is1Dim) nominal_h_m[flavId] = (TH1D*)(inputF->Get(nameStr.c_str()))->Clone();
        else nominal_h_m[flavId] = (TH2D*)(inputF->Get(nameStr.c_str()))->Clone();
      }

      for(uint s_i=0; s_i < sys_v.size(); ++s_i){
        string sysStr = sys_v[s_i];
        nameStr =  baseName+"_"+flavStr+"_"+sysStr;

        if(inputFile == ""){
          string labelX = exampleHist->GetXaxis()->GetTitle();
          string labelY = exampleHist->GetYaxis()->GetTitle();
          if(is1Dim) {
            sys_h_m[sys_v[s_i]] = (TH1D*)exampleHist->Clone(nameStr.c_str());
            sys_corl_h_m[sys_v[s_i]] = (TH1D*)exampleHist->Clone((nameStr+"_corl").c_str());
            sys_uncorl_h_m[sys_v[s_i]] = (TH1D*)exampleHist->Clone((nameStr+"_uncorl").c_str());
          } else {
            sys_h_m[sys_v[s_i]] = (TH2D*)exampleHist->Clone(nameStr.c_str());
            sys_corl_h_m[sys_v[s_i]] = (TH2D*)exampleHist->Clone((nameStr+"_corl").c_str());
            sys_uncorl_h_m[sys_v[s_i]] = (TH2D*)exampleHist->Clone((nameStr+"_uncorl").c_str());
          }
          sys_h_m[sys_v[s_i]]->SetTitle(nameStr.c_str());
          sys_h_m[sys_v[s_i]]->GetXaxis()->SetTitle( (labelX+" ("+flavStr+" jets)").c_str());
          sys_corl_h_m[sys_v[s_i]]->SetTitle((nameStr+"_corl").c_str());
          sys_corl_h_m[sys_v[s_i]]->GetXaxis()->SetTitle( (labelX+" ("+flavStr+" jets)").c_str());
          sys_uncorl_h_m[sys_v[s_i]]->SetTitle((nameStr+"_uncorl").c_str());
          sys_uncorl_h_m[sys_v[s_i]]->GetXaxis()->SetTitle( (labelX+" ("+flavStr+" jets)").c_str());
          if(is1Dim) {
            sys_h_m[sys_v[s_i]]->GetYaxis()->SetTitle( "Entries");
            sys_corl_h_m[sys_v[s_i]]->GetYaxis()->SetTitle( "Entries");
            sys_uncorl_h_m[sys_v[s_i]]->GetYaxis()->SetTitle( "Entries");
          } else {
            sys_h_m[sys_v[s_i]]->GetYaxis()->SetTitle( (labelY+" ("+flavStr+" jets)").c_str());
            sys_h_m[sys_v[s_i]]->GetZaxis()->SetTitle( "Entries");
            sys_corl_h_m[sys_v[s_i]]->GetYaxis()->SetTitle( (labelY+" ("+flavStr+" jets)").c_str());
            sys_corl_h_m[sys_v[s_i]]->GetZaxis()->SetTitle( "Entries");
            sys_uncorl_h_m[sys_v[s_i]]->GetYaxis()->SetTitle( (labelY+" ("+flavStr+" jets)").c_str());
            sys_uncorl_h_m[sys_v[s_i]]->GetZaxis()->SetTitle( "Entries");
          }
        } else {
          //Initialization from saved file. Inputs MUST be consistent
          if(is1Dim) {
            sys_h_m[sys_v[s_i]] = (TH1D*)(inputF->Get(nameStr.c_str()))->Clone();
            sys_corl_h_m[sys_v[s_i]] = (TH1D*)(inputF->Get((nameStr+"_corl").c_str()))->Clone();
            sys_uncorl_h_m[sys_v[s_i]] = (TH1D*)(inputF->Get((nameStr+"_uncorl").c_str()))->Clone();
          } else {
            sys_h_m[sys_v[s_i]] = (TH2D*)(inputF->Get(nameStr.c_str()))->Clone();
            sys_corl_h_m[sys_v[s_i]] = (TH2D*)(inputF->Get((nameStr+"_corl").c_str()))->Clone();
            sys_uncorl_h_m[sys_v[s_i]] = (TH2D*)(inputF->Get((nameStr+"_uncorl").c_str()))->Clone();
          }
        }     
      }
      sys_flav_h_m[flavId] = sys_h_m;
      sys_flav_corl_h_m[flavId] = sys_corl_h_m;
      sys_flav_uncorl_h_m[flavId] = sys_uncorl_h_m;
    }
    if(inputF) inputF->Close();
  };


  void fillNominal(int flav, double w, double fill_x, double fill_y=-99999){
    TH1D* h1_tmp = nullptr;
    TH2D* h2_tmp = nullptr;
    if(is1Dim) h1_tmp = (TH1D*) nominal_h_m[flav];
    else h2_tmp = (TH2D*)  nominal_h_m[flav];

    if(is1Dim) h1_tmp->Fill(fill_x, w);
    else h2_tmp->Fill(fill_x, fill_y, w);
  }

  void fillNominal_wFlags(int flav, double w, int LFCalibType_flag,
                          double fill_x, double fill_y=-99999){

    fillNominal(flav, w, fill_x, fill_y);

#if  _DOD0TailsJetSF == 1
    fillSys(flav, gD0TailSF*w, "D0Tails", fill_x, fill_y);//++++++++TESTING WITH A HORRIBLE GLOBAL...
    fillSys(flav, gD0TailSF_UP*w, "D0Tails_UP", fill_x, fill_y);//++++++++TESTING WITH A HORRIBLE GLOBAL...
    fillSys(flav, gD0TailSF_DOWN*w, "D0Tails_DOWN", fill_x, fill_y);//++++++++TESTING WITH A HORRIBLE GLOBAL...
#endif

    if(LFCalibType_flag >-1) {//Do nothing if the flag is <0
      const int KShort_OR_Lambda = 0x1 | 0x2;
      const int Conversion_OR_HadMatInterac = 0x4 | 0x8;
      
      if( (LFCalibType_flag & KShort_OR_Lambda) ){
        fillSys(flav, 1.3*w, "TYPE_SHadrons", fill_x, fill_y);
      } else {
        fillSys(flav, 1.0*w, "TYPE_SHadrons", fill_x, fill_y);
      }
      
      if( (LFCalibType_flag & Conversion_OR_HadMatInterac) ){
        fillSys(flav, 1.1*w, "TYPE_ConvOrHadInt", fill_x, fill_y);
      } else {
        fillSys(flav, 1.0*w, "TYPE_ConvOrHadInt", fill_x, fill_y);
      }
    }
  }

  void fillSys_wCorl(int flav, double w, string sysStr, bool isCorl_event, double fill_x, double fill_y=-99999){
    fillSys(flav, w, sysStr, fill_x, fill_y);
    if(isCorl_event){
      TH1D* h1_tmp = nullptr;
      TH2D* h2_tmp = nullptr;
      if(is1Dim) h1_tmp = (TH1D*) sys_flav_corl_h_m[flav][sysStr];
      else h2_tmp = (TH2D*) sys_flav_corl_h_m[flav][sysStr];

      if(is1Dim) h1_tmp->Fill(fill_x, w);
      else h2_tmp->Fill(fill_x, fill_y, w);
    } else {
      TH1D* h1_tmp = nullptr;
      TH2D* h2_tmp = nullptr;
      if(is1Dim) h1_tmp = (TH1D*) sys_flav_uncorl_h_m[flav][sysStr];
      else h2_tmp = (TH2D*) sys_flav_uncorl_h_m[flav][sysStr];

      if(is1Dim) h1_tmp->Fill(fill_x, w);
      else h2_tmp->Fill(fill_x, fill_y, w);
    }
  }

  void fillSys(int flav, double w, string sysStr, double fill_x, double fill_y=-99999){
    TH1D* h1_tmp = nullptr;
    TH2D* h2_tmp = nullptr;

    if(is1Dim) h1_tmp = (TH1D*) sys_flav_h_m[flav][sysStr];
    else h2_tmp = (TH2D*) sys_flav_h_m[flav][sysStr];

    if(is1Dim) h1_tmp->Fill(fill_x, w);
    else h2_tmp->Fill(fill_x, fill_y, w);
  }

  void makeVar_h(){
    for(auto nominal_i : nominal_h_m){//Loop on flavors
      int flav_id = nominal_i.first;
      TH1* nominal_h = nominal_i.second;

      map<string, TH1*> sys_variation_h_m;
      map<string, TH1*> sys_h_m = sys_flav_h_m[flav_id];
      map<string, TH1*> sys_corl_h_m = sys_flav_corl_h_m[flav_id];
      map<string, TH1*> sys_uncorl_h_m = sys_flav_uncorl_h_m[flav_id];

      for(auto sys_i : sys_h_m){
        string sysStr = sys_i.first;
        string nameStr = sys_i.second->GetName(); nameStr += "_var";
        TH1* sys_h = sys_i.second;
        //  cout<<sysStr<<" for "<<nameStr<<endl;

        //Prepare syst variation hist
        if(is1Dim) sys_variation_h_m[sysStr] = (TH1D*) sys_h->Clone(nameStr.c_str());
        else sys_variation_h_m[sysStr] = (TH2D*) sys_h->Clone(nameStr.c_str());
        sys_variation_h_m[sysStr]->SetTitle(nameStr.c_str());
        if(is1Dim) sys_variation_h_m[sysStr]->GetYaxis()->SetTitle("SF");
        else sys_variation_h_m[sysStr]->GetZaxis()->SetTitle("SF");

        sys_variation_h_m[sysStr]->Reset();

        //Now do the SF_var= sys_h/nominal_h
        // NB: using the basic division  sys_variation_h_m[sysStr]->Divide(nominal_h);
        // is not correct for stat errs because all events in sys_h are also in nominal_h
        for(int xi=0; xi <= sys_variation_h_m[sysStr]->GetNbinsX(); ++xi){
          for(int yi=0; yi <= sys_variation_h_m[sysStr]->GetNbinsY(); ++yi){

            // NB: I might have weighted events so that the bin content is n-entries
            // and the corresponding stat error is sqrt(sum weight^2)
            double sys_Value = sys_h->GetBinContent(xi,yi);
            double nom_Value = nominal_h->GetBinContent(xi,yi);
            double sys_StatEr = sys_h->GetBinError(xi,yi);
            double nom_StatEr = nominal_h->GetBinError(xi,yi);

            if(nom_Value != 0){
              double sys_var = sys_Value/nom_Value;
              sys_variation_h_m[sysStr]->SetBinContent(xi, yi, sys_var);

              double corl_Value = sys_corl_h_m[sysStr]->GetBinContent(xi,yi);
              double uncorl_Value = sys_uncorl_h_m[sysStr]->GetBinContent(xi,yi);
              double corl_StatEr = sys_corl_h_m[sysStr]->GetBinError(xi,yi);
              double uncorl_StatEr = sys_uncorl_h_m[sysStr]->GetBinError(xi,yi);

              double sys_var_err_uncorl = (uncorl_Value) == 0 ? 0 :
                (uncorl_Value/nom_Value) * sqrt( uncorl_StatEr*uncorl_StatEr/(uncorl_Value*uncorl_Value) 
                                                 + nom_StatEr*nom_StatEr/(nom_Value*nom_Value) );

              double sys_var_err_corl  = (corl_Value) == 0 ? 0 :
                (corl_Value/nom_Value) * sqrt( corl_StatEr*corl_StatEr/(corl_Value*corl_Value) 
                                               + nom_StatEr*nom_StatEr/(nom_Value*nom_Value)
                                               - 2* corl_StatEr*nom_StatEr/(corl_Value*nom_Value) );

              double sys_var_err = sqrt(sys_var_err_corl*sys_var_err_corl + sys_var_err_uncorl*sys_var_err_uncorl);
              sys_variation_h_m[sysStr]->SetBinError(xi, yi, sys_var_err);
              //cout<<"sys_var = "<<sys_var<<endl<<"sys_var_err = "<<sys_var_err<<" , sys_var_err_uncorl = "<<sys_var_err_uncorl<<" , sys_var_err_corl = "<<sys_var_err_corl<<endl;
            }
          }
        }
      }
      sys_flav_variation_h_m[flav_id] = sys_variation_h_m;
    }
  };


  void write_h(TFile* outFile, const bool wVariation= true){
    outFile->cd();
    for(auto nominal_i : nominal_h_m){//loop on flavors
      int flav_id = nominal_i.first;
      nominal_i.second->Write();

      map<string, TH1*> sys_h_m = sys_flav_h_m[flav_id];
      for(auto sys_i : sys_h_m){//loop on sys
        string sysStr = sys_i.first;
        sys_flav_h_m[flav_id][sysStr]->Write();
        sys_flav_corl_h_m[flav_id][sysStr]->Write();
        sys_flav_uncorl_h_m[flav_id][sysStr]->Write();
        if(wVariation) sys_flav_variation_h_m[flav_id][sysStr]->Write();
      }
    }
  };

};

//+++++++ TESTED in
// eosmount eosLink
// asetup 20.7.8.7, here
//++++++++++++

//mv2c10 working points listed here:
// https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks#MV2c10_tagger_added_on_11th_May
// FixedCutBEff_30  0.997715  29.99   99.95   1406.40   4891.95   51290.94
// FixedCutBEff_50  0.976933  50.05   99.62   106.41  614.79  5701.28
// FixedCutBEff_60  0.934906  60.03   99.00   34.54   183.98  1538.78
// FixedCutBEff_70  0.824427  69.97   97.46   12.17   54.72   381.32
// FixedCutBEff_77  0.645925  76.97   95.17   6.21  22.04   134.34
// FixedCutBEff_85  0.175847  84.95   89.66   3.10  8.17  33.53

void effPlot_newBinning(string outFileID = "effPlot_FixedCutBEff_60", const double bTagCut = 0.934906,
                        string filesDir =  "/hepstore/awychan/LF_test_1117/ttbar/",
                        const int nEventMax = -1,
                        const float MV2c100cut= -999., const float MV2cl100cut= -999.,
                        bool isJZxW_Py8 = true, //if false then is Herwig
                        int nSplit = -1, int nSubJ = -1, string fileToMergeInput = ""
                        ){

  TH1::SetDefaultSumw2(kTRUE);
  TH1::AddDirectory(kFALSE);
  const int nEventPrint = 1000;
  //Light Flavor Studies
  std::cout << "VHF Light Flavor Studies" << std::endl;

  const bool doC = false, doB = false;
  const bool doLFCalibType_sys = false;//Based on the TYPE flagging of entire jets
  const bool studyCorrelation = true;
  const bool doJVT = true;
  const bool useMCWeights = false;
  const int skipJZ0and1 = 0;//0 = no skip, 1=JZ0 skip, 2= JZ0,JZ1 skip
  const int doOnlyDSID = -1;//if != -1 then only this DSID is analyzed
  const int  checkPU = 1;//0 = no check,  1 = check for DR<0.6, 2 = check for DR<0.3
  const int verb = 0;

  if(doOnlyDSID > -1) cout<<endl<<"WARNING: doing only DSID = "<<doOnlyDSID<<endl<<endl;

  //Trying script parallel processing 
  string subJStr = ""; 
  bool isMultiJob = false;
  if(nSplit != -1 || nSubJ != -1){
    cout<<"Multi-job processing: nSplit = "<<nSplit<<" nSubJ = "<<nSubJ<<endl;
    if( nSplit < 2 || nSubJ < 1 || nSubJ>nSplit ) {cerr<<"Config ERROR!"<<endl; exit(-1);}
    subJStr = Form("_JobN%dof%d",nSubJ, nSplit);
    isMultiJob = true;
  }

  clock_t startTime = clock();
  cout<<endl<<"--> Start filling hists at: "<<double(startTime)/CLOCKS_PER_SEC<<" <--"<<endl;
  TFile outF((outFileID+subJStr+".root").c_str(), "recreate");
  TChain stdFile("bTag_AntiKt4EMTopoJets");
  //TChain stdFile("flav_Akt4EMTo");
  ScanDir ( filesDir.c_str(), &stdFile, "Nominal");
  cout <<"Dir: " << filesDir.c_str() << endl;

  vector<string> trk_sys_v = {
    "RES_D0_MEAS", "RES_Z0_MEAS",
    "FAKE_RATE_LOOSE",
    //"FAKE_RATE_TIGHT",
    //"RES_D0Z0Corl_MEAS",
    //"TAIL_RATE_SYS",
    //"TAIL_PAR_SYS",
    //"RES_D0_MEAS_UP", "RES_Z0_MEAS_UP",
    //"RES_D0_MEAS_DOWN", "RES_Z0_MEAS_DOWN",
    //"RES_D0_MEAS_Zmm", "RES_Z0_MEAS_Zmm",
    //"RES_D0Z0_MEAS",
    //"CONV_HADINT_RATE",
    "STRANGEHAD_RATE",
    "TAIL_RATE",

    //"CONV_HADINT_RATE_SYS",
    //"STRANGEHAD_RATE_SYS",
    // // "EFF_LOOSE_TIDE",
    // // "BIAS_D0_WM", "BIAS_Z0_WM", "BIAS_QOVERP_SAGITTA_WM",
    // // "EFF_LOOSE_GLOBAL", "EFF_LOOSE_IBL", "EFF_LOOSE_PP0",
    // // "EFF_TIGHT_GLOBAL", "EFF_TIGHT_IBL", "EFF_TIGHT_PP0",
    // // "RES_D0_DEAD", "RES_Z0_DEAD",
  };

  if(!isJZxW_Py8){ //For Herwig studies
    trk_sys_v.clear();
    trk_sys_v.push_back("RES_D0_MEAS");
    trk_sys_v.push_back("RES_Z0_MEAS");
    trk_sys_v.push_back("TAIL_RATE");
    trk_sys_v.push_back("FAKE_RATE_LOOSE");

  }

  if(doLFCalibType_sys) {
    //Consider also strange-hadron +30% and hadronic material interaction or Conversions + 10%
    trk_sys_v.push_back("TYPE_SHadrons");
    trk_sys_v.push_back("TYPE_ConvOrHadInt");
#if  _DOD0TailsJetSF == 1
    trk_sys_v.push_back("D0Tails");
    trk_sys_v.push_back("D0Tails_UP");
    trk_sys_v.push_back("D0Tails_DOWN");
#endif
  }
  
  //For each DSID, fill all events and jets you see in nominal sample (use event_num*10^3 + jet_i)
  // and compare if they appear also in the syst.
  // Correlated and uncorr stat uncertainties can be done looking at it
  map<int, set<long int> >same_event_checker;

  double ptbins[] = {20, 30, 60, 100, 180, 300, 500, 1000, 3000}; //to be used for both 1D and 2D
  int   n_jetpt_bins = sizeof(ptbins)/sizeof(double) - 1;

  //Corresponding to Matthias' negative-LF calibration
  double ptbinsRebin[] = {20, 60, 100, 200, 300, 500, 750, 1000, 3000}; //to be used for both 1D and 2D
  int   n_jetptRebin_bins = sizeof(ptbinsRebin)/sizeof(double) - 1;

  double etabins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.5}; //to be used for both 1D and 2D
  int   n_jeteta_bins = sizeof(etabins)/sizeof(double) - 1;

  double etabinsRebin[] = {0.0,1.2,2.5}; //to be used for both 1D and 2D
  int   n_jetetaRebin_bins = sizeof(etabinsRebin)/sizeof(double) - 1;

  int n_jetpt_bins_2D  =  n_jetpt_bins;
  int n_jeteta_bins_2D =  n_jeteta_bins;
  int n_jetptRebin_bins_2D  =  n_jetptRebin_bins;
  int n_jetetaRebin_bins_2D =  n_jetetaRebin_bins;

  //  double mv2c10_bins[] = {-1.0, 0.175847, 0.645925, 0.824427, 0.934906, 0.976933, 0.997715, 1};
  double mv2c10_bins[] = {-1.0, 0.175847, 0.645925, 0.824427, 0.934906, 0.976933, 1};//No 30% WP
  int n_mv2c10_bins = sizeof(mv2c10_bins)/sizeof(double) - 1;

  map<int, string> flav_map;
  flav_map[LF_id] = "LF";
  if(doB) flav_map[B_id] = "B";
  if(doC) flav_map[C_id] = "C";

  TH1D *mv2c10_h   = new TH1D("mv2c10",  "mv2c10:mv2c10", 202, -1.01, 1.01);
  calib_h mv2c10_ch(flav_map, trk_sys_v, mv2c10_h, fileToMergeInput);

  TH1D *mv2c10wp_h   = new TH1D("mv2c10wp",  "mv2c10wp:mv2c10",  n_mv2c10_bins, mv2c10_bins);
  calib_h mv2c10wp_ch(flav_map, trk_sys_v, mv2c10wp_h, fileToMergeInput);

  TH1D *pt_h   = new TH1D("pt",  "pt;jet p_{T} [GeV]", n_jetpt_bins, ptbins);
  calib_h pt_ch(flav_map, trk_sys_v, pt_h, fileToMergeInput);

  TH1D *eta_h  = new TH1D("eta", "eta;jet #eta", n_jeteta_bins, etabins);
  calib_h eta_ch(flav_map, trk_sys_v, eta_h, fileToMergeInput);

  TH2D *etapt_h   = new TH2D("etapt",  "etapt;jet p_{T} [GeV];jet #eta", n_jetpt_bins_2D, ptbins, n_jetetaRebin_bins_2D, etabinsRebin );
  calib_h etapt_ch(flav_map, trk_sys_v, etapt_h, fileToMergeInput);

  //Use best pT binning for pseudo-continous calib
  TH2D *cenEta_pt_mv2c10wp_h   = new TH2D("cenEta_pt_mv2c10wp", "cenEta_pt_mv2c10wp;jet p_{T} [GeV];mv2c10", n_jetpt_bins_2D, ptbins,  n_mv2c10_bins, mv2c10_bins);
  calib_h cenEta_pt_mv2c10wp_ch(flav_map, trk_sys_v, cenEta_pt_mv2c10wp_h, fileToMergeInput);
  TH2D *frwEta_pt_mv2c10wp_h   = new TH2D("frwEta_pt_mv2c10wp", "frwEta_pt_mv2c10wp;jet p_{T} [GeV];mv2c10", n_jetpt_bins_2D, ptbins,  n_mv2c10_bins, mv2c10_bins);
  calib_h frwEta_pt_mv2c10wp_ch(flav_map, trk_sys_v, frwEta_pt_mv2c10wp_h, fileToMergeInput);

  TH1D *ptFine_h   = new TH1D("ptFine",  "ptFine;jet p_{T} [GeV]", 300, 0, 1200);
  calib_h ptFine_ch(flav_map, trk_sys_v, ptFine_h, fileToMergeInput);

#if _MORE_H > 0
  TH2D *etaptRebin_h = new TH2D("etaptRebin", "etaptRebin;jet p_{T} [GeV];jet #eta", n_jetptRebin_bins_2D, ptbinsRebin, n_jetetaRebin_bins_2D, etabinsRebin );
  calib_h etaptRebin_ch(flav_map, trk_sys_v, etaptRebin_h, fileToMergeInput);

  TH1D *mv2c10Flip_h   = new TH1D("mv2c10Flip",  "mv2c10Flip;mv2c10 flip", 102, -1.02, 1.02);
  calib_h mv2c10Flip_ch(flav_map, trk_sys_v, mv2c10Flip_h, fileToMergeInput);

  TH1D *mv2c10wpFlip_h   = new TH1D("mv2c10wpFlip",  "mv2c10wpFlip;mv2c10 flip",   n_mv2c10_bins, mv2c10_bins);
  calib_h mv2c10wpFlip_ch(flav_map, trk_sys_v, mv2c10wpFlip_h, fileToMergeInput);

  //These are for Matthias, I will use his binning
  TH2D *cenEta_pt_mv2c10wpFlip_h   = new TH2D("cenEta_pt_mv2c10wpFlip", "cenEta_pt_mv2c10wpFlip;jet p_{T} [GeV];mv2c10Flip", n_jetptRebin_bins_2D, ptbinsRebin,  n_mv2c10_bins, mv2c10_bins);
  calib_h cenEta_pt_mv2c10wpFlip_ch(flav_map, trk_sys_v, cenEta_pt_mv2c10wpFlip_h, fileToMergeInput);
  TH2D *frwEta_pt_mv2c10wpFlip_h   = new TH2D("frwEta_pt_mv2c10wpFlip", "frwEta_pt_mv2c10wpFlip;jet p_{T} [GeV];mv2c10Flip", n_jetptRebin_bins_2D, ptbinsRebin,  n_mv2c10_bins, mv2c10_bins);
  calib_h frwEta_pt_mv2c10wpFlip_ch(flav_map, trk_sys_v, frwEta_pt_mv2c10wpFlip_h, fileToMergeInput);
#endif
#if _MORE_H > 1
  TH1D *ip3d_llr_h   = new TH1D("ip3d_llr", "ip3d_llr;ip3d llr", 100, -10, 30);
  calib_h ip3d_llr_ch(flav_map, trk_sys_v, ip3d_llr_h, fileToMergeInput );

  TH1D *sv1_llr_h   = new TH1D("sv1_llr", "sv1_llr;sv1 llr", 100, -10, 30);
  calib_h sv1_llr_ch(flav_map, trk_sys_v, sv1_llr_h, fileToMergeInput);

  TH1D *jf_llr_h   = new TH1D("jf_llr", "jf_llr;jf llr", 100, -10, 30);
  calib_h jf_llr_ch(flav_map, trk_sys_v, sv1_llr_h, fileToMergeInput);
#endif
#if _DOTRK == 1
  TH1D *ntrk_h     = new TH1D("ntrk", "ntrk; n. tracks", 40, -0.5, 39.5);
  calib_h ntrk_ch(flav_map, trk_sys_v, ntrk_h, fileToMergeInput);

  TH1D *trk_pt_h   = new TH1D("trk_pt_LF_std", "trk_pt;leading track p_{T} [GeV]", 100, 0., 100);
  calib_h trk_pt_ch(flav_map, trk_sys_v, trk_pt_h, fileToMergeInput);

  TH1D *trk_eta_h  = new TH1D("trk_eta_LF_std", "trk_eta;leading track #eta", 50, -2.5, 2.5);
  calib_h trk_eta_ch(flav_map, trk_sys_v, trk_eta_h, fileToMergeInput);

  TH1D *trk_d0_h   = new TH1D("trk_d0", "trk_d0;leading track d_{0}^{BL} [mm]", 100, -2.5, 2.5);
  calib_h trk_d0_ch(flav_map, trk_sys_v, trk_d0_h, fileToMergeInput);

  TH1D *trk_z0_h   = new TH1D("trk_z0", "trk_z0;leading track z_{0}^{BL} [mm]", 60, -3., 3.);
  calib_h trk_z0_ch(flav_map, trk_sys_v, trk_z0_h, fileToMergeInput);
#endif
  
  if(fileToMergeInput == "") { //Skip all and just go to the merging phase
    TTreeReader reader(&stdFile);
    TTreeReaderValue< int                   >  eventnb         (reader, "eventnb");
    TTreeReaderValue< int                   >  mcchan          (reader, "mcchan");
    //TTreeReaderValue< int                   >  njets           (reader, "njets");
    TTreeReaderValue< float                 >  mcwg            (reader, "mcwg");
    TTreeReaderValue< vector<float>         >  jet_pt          (reader, "jet_pt");
    TTreeReaderValue< vector<float>         >  jet_eta         (reader, "jet_eta");
    TTreeReaderValue< vector<float>         >  jet_JVT         (reader, "jet_JVT");
    TTreeReaderValue< vector<int>           >  jet_isPU        (reader, "jet_isPU");
    TTreeReaderValue< vector<int>           >  jet_truthMatch  (reader, "jet_truthMatch");
    //TTreeReaderValue< vector<int>           >  jet_LFCalibType (reader, "jet_LFCalibType");
    TTreeReaderValue< vector<int>           >  jet_LabDr_HadF  (reader, "jet_LabDr_HadF");
    TTreeReaderValue< vector<int>           >  jet_aliveAfterOR(reader, "jet_aliveAfterOR");
    TTreeReaderValue< vector<int>           >  jet_aliveAfterORmu(reader, "jet_aliveAfterORmu");
    TTreeReaderValue< vector<double>         >  jet_mv2c10      (reader, "jet_mv2c10");
    TTreeReaderValue< vector<double>         >  jet_mv2c100      (reader, "jet_mv2c100");
    TTreeReaderValue< vector<double>         >  jet_mv2cl100     (reader, "jet_mv2cl100");
#if _MORE_H > 0
    TTreeReaderValue< vector<float>         >  jet_mv2c10Flip  (reader, "jet_mv2c10Flip");
#endif
#if _MORE_H > 1
    TTreeReaderValue< vector<float>         >  jet_ip3d_llr    (reader, "jet_ip3d_llr");
    TTreeReaderValue< vector<float>         >  jet_sv1_llr     (reader, "jet_sv1_llr");
    TTreeReaderValue< vector<float>         >  jet_jf_llr      (reader, "jet_jf_llr");
#endif
#if _DOTRK == 1
    TTreeReaderValue< vector<int>           >  jet_btag_ntrk   (reader, "jet_btag_ntrk");
    TTreeReaderValue< vector<int>           >  jet_ip3d_ntrk   (reader, "jet_ip3d_ntrk");
#endif
#if _DOD0TailsJetSF == 1
    TTreeReaderValue< vector<vector<float>> >  jet_trk_pt      (reader, "jet_trk_pt");
    TTreeReaderValue< vector<vector<float>> >  jet_trk_eta     (reader, "jet_trk_eta");
    TTreeReaderValue< vector<vector<float>> >  jet_trk_d0      (reader, "jet_trk_d0");
    TTreeReaderValue< vector<vector<float>> >  jet_trk_z0      (reader, "jet_trk_z0");
#endif

    const float track_pt_cut = 0.7; //cut on pt leading track (higher than 700 MeV)
    double j_pt, j_eta, mv2c10, mv2c10Flip, mc_w;
    int j_i, had_label, mcDSID;
    int LFCalibType;

    //Stuff useful for tracks
    int nTracks, nTrkTails, ptEtaBin;
    double  trk_pt, trk_eta, trk_d0, trk_z0;
    double totD0TailSF, maxD0TailSF, minD0TailSF;

    map<int, double> JZxW_scale_m;
    char* workDirChar = getenv("TestArea");
    string workDir= "/user/egraham/ljets_MCCalibrationWorks/";
    //string workDir= "/hepstore/awychan/PhD/QT/LFUpdate_20171119/ljets_MCCalibration/";
    if(useMCWeights){
      string fileStr = "norm_diJet_mc_py8";
      int fistDSID = 361020;
      if(isJZxW_Py8){  //defaults
        cout<<"Normalization factors for Py8 MC..."<<endl;
      } else {
        cout<<"Normalization factors for Hw MC..."<<endl;
        fileStr = "norm_diJet_mc_hw";
        fistDSID = 426040;
      }
      TFile xsecFile((workDir+"/macros/"+fileStr+".root").c_str());
      TH1D* normScale_h = (TH1D*) xsecFile.Get("hmc_normScale");
      TH1D* xSection_h = (TH1D*) xsecFile.Get("hmc_xSection");
      for(int bi=0; bi<13; ++bi){
        JZxW_scale_m[ fistDSID+bi ] = normScale_h->GetBinContent(bi+1);
        cout<<"DSID = "<< (fistDSID+bi)<<" xSection "<<xSection_h->GetBinContent(bi+1)
            <<", normScale = "<<normScale_h->GetBinContent(bi+1)<<endl;
      }
    }

#if _DOD0TailsJetSF == 1
    TH2F* d0Tail_h, *d0OneSigmaCore_h;
    TFile tailFile((workDir+"/macros/average_tails_d0reco_1p5to4_RMS_Rebin2_allJets_maxTrkPT100.root").c_str());
    d0Tail_h = (TH2F*)(tailFile.Get("d0reco_PtEta_tails_ratio"))->Clone();
    d0OneSigmaCore_h = (TH2F*)(tailFile.Get("d0reco_PtEta_oneSigmaDistance"))->Clone();
#endif

    long int totEntries = (nEventMax > -1) ? nEventMax : reader.GetEntries(true);

    long int jobStart = 0, jobEnd = totEntries, subJNEntries = 0, subJRemainder =0;
    if(isMultiJob) { 
      subJNEntries = totEntries/nSplit;
      subJRemainder = totEntries % nSplit;
      jobStart = subJNEntries*(nSubJ - 1);
      jobEnd = (nSplit != nSubJ) ? subJNEntries*nSubJ : subJNEntries*nSubJ +subJRemainder;
      cout<<"Multi-job processing: jobStart = "<<jobStart<<", jobEnd = "<<jobEnd
          <<" (subJNEntries="<<subJNEntries<<" subJRemainder = "<<subJRemainder<<")"<<endl;
    }

    long int nEvent = 0;
    ///// EVENT LOOP ////////////
    for( nEvent=jobStart; nEvent<jobEnd; ++nEvent) {
      if ((nEvent%nEventPrint) == 00) cout<< "nEvent=" <<nEvent <<" over "<<totEntries<<endl;
      if(reader.SetEntry(nEvent) != TTreeReader::kEntryValid) {
        cerr<<"WARNINIG: entry "<<nEvent<<" isn't good!"<<endl;
        continue;
      }

      mcDSID = (*mcchan);
      if(verb) cout<<"DSID: "<<mcDSID<<" ev_w = "<<(*mcwg)<<" sampleNorm = "<<JZxW_scale_m[ mcDSID ]<<endl;
      if(doOnlyDSID != -1 && doOnlyDSID != mcDSID) continue;
      if(skipJZ0and1) {
        if(skipJZ0and1 > 0 && (mcDSID == 361020 || mcDSID == 426040)) continue;
        if(skipJZ0and1 > 1 && (mcDSID == 361021 || mcDSID == 426041)) continue;
      }

      mc_w = useMCWeights ? JZxW_scale_m[ mcDSID ] * (*mcwg) : 1.0;

      int njets = jet_pt->size();

      for(j_i  =0; j_i < njets; ++j_i){
        if(verb) cout<<"j_i = "<<j_i<<" truth_match = "<<(*jet_truthMatch)[j_i]<<" pT = "<<((*jet_pt)[j_i]/1000.)
                     <<" mcv2c10 = "<<(*jet_mv2c10)[j_i]<<" Flav = "<<(*jet_LabDr_HadF)[j_i]<<endl;
        if(checkPU >0 && (*jet_isPU)[j_i]) continue;
        if(checkPU >1 && !(*jet_truthMatch)[j_i]) continue;

        //Needed when analyzing samples with truth-ele or mu, which are removed
        if(!(*jet_aliveAfterOR)[j_i]) continue;
        if(!(*jet_aliveAfterORmu)[j_i]) continue;

        j_pt = (*jet_pt)[j_i]/1000.;
        j_eta = fabs( (*jet_eta)[j_i]);

        if(j_pt<20. || j_eta>2.5) continue;//Standard calibration cuts
        //Jet cleaning probably should be applied only to not-bad jets but the fraction is minimal on MC
        if(doJVT && j_pt<60. && j_eta<2.4 && (*jet_JVT)[j_i]<0.59) continue;

        had_label = (*jet_LabDr_HadF)[j_i];
        if( !(had_label == LF_id || 
              (doB && had_label == B_id ) ||
              (doC && had_label == C_id ) )) continue;

        mv2c10 = (*jet_mv2c10)[j_i];

        bool passed_cut = mv2c10 > bTagCut;

        if(MV2c100cut>-999){
          float mv2c100 = (*jet_mv2c100)[j_i];
          float mv2cl100 = (*jet_mv2cl100)[j_i];
          passed_cut = (mv2c100 < MV2c100cut && mv2cl100 > MV2cl100cut);
        }

        if(passed_cut){ 
          LFCalibType =  -1;//With -1 nothing gets filled except nominal
          //if(doLFCalibType_sys) LFCalibType = (*jet_LFCalibType)[j_i];
          if(studyCorrelation) same_event_checker[mcDSID].insert( (*eventnb)*1000 + j_i);

#if _DOD0TailsJetSF == 1
          //Check each track (used for b-tagging) matched to a jet and see if, for its pT/eta,
          // the reco-d0 of the track is classified as "tail" i.e. it lies beyond nSigma from the core of the d0
          // resolution in that bin. At that point look at the tail-SF for such track.
          // SF for different tracks are multiplied. The final result is taken as the SF for this jet
          const double nSigmaTail = 2.0;//Number of tails one
        
          //Reset stuff
          totD0TailSF = 0.; gD0TailSF = 1., gD0TailSF_UP = 1., gD0TailSF_DOWN = 10.;
          nTrkTails = 0; ptEtaBin = -1; 
          nTracks = (*jet_trk_pt)[j_i].size();
        
          for(uint trk_i =0; trk_i<(uint)nTracks; ++trk_i){
            trk_pt  = (*jet_trk_pt)[j_i][trk_i]/1000.;
            trk_eta  = (*jet_trk_eta)[j_i][trk_i];
            trk_d0   = (*jet_trk_d0)[j_i][trk_i];
            trk_z0   = (*jet_trk_z0)[j_i][trk_i];
            if(trk_pt > 99.) trk_pt = 99.; //limits on tail hist file
            if(trk_eta > 2.499) trk_eta = 2.499; 
            if(trk_eta < -2.499 ) trk_eta = -2.499; 
          
            ptEtaBin = d0OneSigmaCore_h->FindBin(trk_pt, trk_eta);
            //cout<<j_i<<" "<<trk_i<<" "<<trk_d0<<" "<<ptEtaBin<<" "<<d0OneSigmaCore_h->GetBinContent(ptEtaBin)<<endl;
            if(ptEtaBin<0) continue;
            if(fabs(trk_d0) < nSigmaTail*d0OneSigmaCore_h->GetBinContent(ptEtaBin)) continue;
            double localTailSF = d0Tail_h->GetBinContent(ptEtaBin);
            if(localTailSF > gD0TailSF_UP) gD0TailSF_UP = localTailSF;
            if(localTailSF < gD0TailSF_DOWN) gD0TailSF_DOWN = localTailSF;
          
            totD0TailSF += d0Tail_h->GetBinContent(ptEtaBin);
            nTrkTails += 1;
            // cout<<nTrkTails<<" "<<totD0TailSF<<endl;
          }
          //I tried the simple multiplication of the SF for each track but the resulting LF calib is huge
          // the mistake in that I am not weighting tracks by their importance in the tagging algorithms, 
          // an better approach is probably to take the average tail SF
          if(nTrkTails>0) {
            gD0TailSF = totD0TailSF/nTrkTails;//++++++++TESTING WITH AN HORRIBLE GLOBAL...
            gD0TailSF_UP = gD0TailSF_UP;
            gD0TailSF_DOWN = gD0TailSF_DOWN;
          } else {
            gD0TailSF = 1.;
            gD0TailSF_UP = 1.;
            gD0TailSF_DOWN = 1.;
          }      
#endif
          pt_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt);
          eta_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_eta );
          etapt_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt, j_eta);
          mv2c10_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, mv2c10);
          mv2c10wp_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, mv2c10);
          if(j_eta<1.2) {
            cenEta_pt_mv2c10wp_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt, mv2c10);
          } else {
            frwEta_pt_mv2c10wp_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt, mv2c10);
          }
          ptFine_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt);
#if _MORE_H > 0
          mv2c10Flip = (*jet_mv2c10Flip)[j_i];
          etaptRebin_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt, j_eta);
          mv2c10Flip_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, mv2c10Flip);
          mv2c10wpFlip_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, mv2c10Flip);
          if(j_eta<1.2) {
            cenEta_pt_mv2c10wpFlip_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt, mv2c10Flip);
          } else {
            frwEta_pt_mv2c10wpFlip_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, j_pt, mv2c10Flip);
          }
#endif
#if _MORE_H > 1
          ip3d_llr_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, (*jet_ip3d_llr)[j_i]);
          sv1_llr_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, (*jet_sv1_llr)[j_i]);
          jf_llr_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, (*jet_jf_llr)[j_i]);
#endif

#if _DOTRK == 1

          float lead_trk_pT  = -999.;
          float lead_trk_eta = -999.;
          float lead_trk_d0  = -999.;
          float lead_trk_z0  = -999.;
          int   trk_count    =  0;

          //for(int j_t  =0; j_t < ((*jet_btag_ntrk)[j_i]); ++j_t){
          //if (!(*jet_ip3d_ntrk)[j_i]) continue;
          for(int j_t  =0; j_t < ((*jet_ip3d_ntrk)[j_i]); ++j_t){
            if(((((*jet_trk_pt)[j_i])).at(j_t)/1000.) > track_pt_cut ){ //require pt leading track to be higher than 700 MeV
              trk_count++;
              if( ((((*jet_trk_pt)[j_i])).at(j_t)/1000.) > lead_trk_pT ){
                lead_trk_pT   = ((((*jet_trk_pt)[j_i])).at(j_t)/1000.);
                lead_trk_eta  = (((*jet_trk_eta)[j_i]).at(j_t));
                lead_trk_d0   = (((*jet_trk_d0)[j_i]).at(j_t));
                lead_trk_z0   = (((*jet_trk_z0)[j_i]).at(j_t));
              }
            }
          }
          ntrk_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, trk_count);
          trk_pt_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, lead_trk_pT);
          trk_eta_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, lead_trk_eta);
          trk_d0_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, lead_trk_d0);
          trk_z0_ch.fillNominal_wFlags(had_label, mc_w, LFCalibType, lead_trk_z0);
#endif
        }
      }
    }
    cout<<endl<<"--> Nominal done at: "<<double(clock())/CLOCKS_PER_SEC<<" <--"<<endl;


    /////////// VARIATIONS //////////
    for(uint i=0; i< trk_sys_v.size(); ++i){
      string sysStr = trk_sys_v[i];
      if(sysStr.find("TYPE") != string::npos ||
         sysStr.find("Tails") != string::npos ) continue; //Deal with these on the nominal jets
      TChain tmpFile("bTag_AntiKt4EMTopoJets");
      ScanDir ( filesDir.c_str(), &tmpFile, sysStr.c_str());

      TTreeReader reader_var(&tmpFile);
      TTreeReaderValue< int                   >  eventnb_var         (reader_var, "eventnb");
      TTreeReaderValue< int                   >  mcchan_var          (reader_var, "mcchan");
      //TTreeReaderValue< int                   >  njets_var           (reader_var, "njets");
      TTreeReaderValue< float                >  mcwg_var            (reader_var, "mcwg");
      TTreeReaderValue< vector<float>         >  jet_pt_var          (reader_var, "jet_pt");
      TTreeReaderValue< vector<float>         >  jet_eta_var         (reader_var, "jet_eta");
      TTreeReaderValue< vector<float>         >  jet_JVT_var         (reader_var, "jet_JVT");
      TTreeReaderValue< vector<int>           >  jet_isPU_var        (reader_var, "jet_isPU");
      TTreeReaderValue< vector<int>           >  jet_truthMatch_var  (reader_var, "jet_truthMatch") ;
      TTreeReaderValue< vector<int>           >  jet_LabDr_HadF_var  (reader_var, "jet_LabDr_HadF");
      TTreeReaderValue< vector<int>           >  jet_aliveAfterOR_var(reader_var, "jet_aliveAfterOR");
      TTreeReaderValue< vector<int>           >  jet_aliveAfterORmu_var(reader_var, "jet_aliveAfterORmu");
      TTreeReaderValue< vector<double>         >  jet_mv2c10_var      (reader_var, "jet_mv2c10");
      TTreeReaderValue< vector<double>         >  jet_mv2c100_var      (reader_var, "jet_mv2c100");
      TTreeReaderValue< vector<double>         >  jet_mv2cl100_var     (reader_var, "jet_mv2cl100");
#if _MORE_H > 0
      TTreeReaderValue< vector<float>         >  jet_mv2c10Flip_var  (reader_var, "jet_mv2c10Flip");
#endif
#if _MORE_H > 1
      TTreeReaderValue< vector<float>         >  jet_ip3d_llr_var    (reader_var, "jet_ip3d_llr");
      TTreeReaderValue< vector<float>         >  jet_sv1_llr_var     (reader_var, "jet_sv1_llr");
      TTreeReaderValue< vector<float>         >  jet_jf_llr_var      (reader_var, "jet_jf_llr");
#endif
#if _DOTRK == 1
      TTreeReaderValue< vector<int>           >  jet_btag_ntrk_var   (reader_var, "jet_btag_ntrk");
      TTreeReaderValue< vector<int>           >  jet_ip3d_ntrk_var   (reader_var, "jet_ip3d_ntrk");
      TTreeReaderValue< vector<vector<float>> >  jet_trk_pt_var      (reader_var, "jet_trk_pt");
      TTreeReaderValue< vector<vector<float>> >  jet_trk_eta_var     (reader_var, "jet_trk_eta");
      TTreeReaderValue< vector<vector<float>> >  jet_trk_d0_var      (reader_var, "jet_trk_d0");
      TTreeReaderValue< vector<vector<float>> >  jet_trk_z0_var      (reader_var, "jet_trk_z0");
#endif
      ///// SYST. EVENT LOOP ////////////
      double j_pt_var, j_eta_var, mv2c10_var, mv2c10Flip_var;
      int j_i, had_label;
      bool isCorlEv = false;
      long int totEntries_var = (nEventMax > -1) ? nEventMax : reader_var.GetEntries(true);

      if(totEntries_var != totEntries)
        cerr<<"WARNINIG: nominal entries differ from variation entries... you will get a bias!"<<endl;

      if(isMultiJob) { 
        subJNEntries = totEntries_var/nSplit;
        subJRemainder = totEntries_var % nSplit;
        jobStart = subJNEntries*(nSubJ - 1);
        jobEnd = (nSplit != nSubJ) ? subJNEntries*nSubJ : subJNEntries*nSubJ +subJRemainder;
        cout<<"Multi-job processing: jobStart = "<<jobStart<<", jobEnd = "<<jobEnd
            <<" (subJNEntries="<<subJNEntries<<" subJRemainder = "<<subJRemainder<<")"<<endl;
      }

      long int nEvent_var = 0;
      for( nEvent_var=jobStart; nEvent_var<jobEnd; ++nEvent_var) {
        if ((nEvent_var % nEventPrint)==00) cout<< "nEvent_var=" <<nEvent_var<<" over "<<totEntries_var<<endl;
        if(reader_var.SetEntry(nEvent_var) != TTreeReader::kEntryValid) {
          cerr<<"WARNINIG: entry "<<nEvent_var<<" isn't good!"<<endl;
          continue;
        }

        mcDSID = (*mcchan_var);
        if(verb) cout<<"DSID: "<<mcDSID<<" ev_w = "<<(*mcwg_var)<<" sampleNorm = "<<JZxW_scale_m[ mcDSID ]<<endl;

        if(doOnlyDSID != -1 && doOnlyDSID != mcDSID) continue;
        if(skipJZ0and1) {
          if(skipJZ0and1 > 0 && (mcDSID == 361020 || mcDSID == 426040)) continue;
          if(skipJZ0and1 > 1 && (mcDSID == 361021 || mcDSID == 426041)) continue;
        }
        
        mc_w = useMCWeights ? JZxW_scale_m[ mcDSID ] * (*mcwg_var) : 1.0;
	int njets_var = jet_pt_var->size();        

        for(int j_i  =0; j_i < njets_var; ++j_i){
          if(verb) cout<<"j_i = "<<j_i<<" truth_match = "<<(*jet_truthMatch_var)[j_i]<<" pT = "<<((*jet_pt_var)[j_i]/1000.)
                       <<" mcv2c10 = "<<(*jet_mv2c10_var)[j_i]<<" Flav = "<<(*jet_LabDr_HadF_var)[j_i]<<endl;

          if(checkPU >0 && (*jet_isPU_var)[j_i]) continue;
          if(checkPU >1 && !(*jet_truthMatch_var)[j_i]) continue;

          //Needed when analyzing samples with truth-ele or mu, which are removed
          if(!(*jet_aliveAfterOR_var)[j_i]) continue;
          if(!(*jet_aliveAfterORmu_var)[j_i]) continue;

          j_eta_var = fabs( (*jet_eta_var)[j_i]);
          j_pt_var = (*jet_pt_var)[j_i]/1000.;

          if(j_pt_var<20. || j_eta_var>2.5) continue;//Standard calibration cuts
          //Jet cleaning probably should be applied only to not-bad jets but the fraction is minimal on MC
          if(doJVT && j_pt_var<60. && j_eta_var<2.4 && (*jet_JVT_var)[j_i] < 0.59) continue;

          had_label = (*jet_LabDr_HadF_var)[j_i];
          if( !(had_label == LF_id || 
                (doB && had_label == B_id) ||
                (doC && had_label == C_id ) )) continue;

          mv2c10_var = (*jet_mv2c10_var)[j_i];

          bool passed_cut = mv2c10_var > bTagCut;

          if(MV2c100cut > -999){          
            float mv2c100_var = (*jet_mv2c100_var)[j_i];
            float mv2cl100_var = (*jet_mv2cl100_var)[j_i];
            passed_cut = (mv2c100_var < MV2c100cut && mv2cl100_var > MV2cl100cut);
          }

          if(passed_cut){
            isCorlEv = false;
            if(studyCorrelation){
              if(same_event_checker[mcDSID].find( (*eventnb_var) * 1000 + j_i ) != same_event_checker[mcDSID].end() ) {
                isCorlEv = true;
              }
            }
            mv2c10_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, mv2c10_var);
            mv2c10wp_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, mv2c10_var);
            pt_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var);
            eta_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_eta_var);
            etapt_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var,  j_eta_var);
            if(j_eta_var<1.2) {
              cenEta_pt_mv2c10wp_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var, mv2c10_var);
            } else {
              frwEta_pt_mv2c10wp_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var, mv2c10_var);            
            }
            ptFine_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var);
#if _MORE_H > 0
            etaptRebin_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var,  j_eta_var);
            mv2c10Flip_var = (*jet_mv2c10Flip_var)[j_i];
            mv2c10Flip_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, mv2c10Flip_var);
            mv2c10wpFlip_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, mv2c10Flip_var);
            if(j_eta_var<1.2) {
              cenEta_pt_mv2c10wpFlip_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var, mv2c10Flip_var);
            } else {
              frwEta_pt_mv2c10wpFlip_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, j_pt_var, mv2c10Flip_var);
            }
#endif
#if _MORE_H > 1
            ip3d_llr_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, (*jet_ip3d_llr_var)[j_i]);
            sv1_llr_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, (*jet_sv1_llr_var)[j_i]);
            jf_llr_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, (*jet_jf_llr_var)[j_i]);
#endif
#if _DOTRK == 1
            float lead_trk_pT  = -999.;
            float lead_trk_eta = -999.;
            float lead_trk_d0  = -999.;
            float lead_trk_z0  = -999.;
            int   trk_count    =  0;
            //if (jet_btag_ntrk_var->size() <= (j_i)) continue;
            //if (!(*jet_ip3d_ntrk_var)[j_i]) continue;
            //for(int j_t  =0; j_t < ((*jet_btag_ntrk_var)[j_i]); ++j_t){
            for(int j_t  =0; j_t < ((*jet_ip3d_ntrk_var)[j_i]); ++j_t){
              if(((((*jet_trk_pt_var)[j_i])).at(j_t)/1000.) > track_pt_cut ){ //require pt leading track to be higher than 700 MeV
                trk_count++;
                if( ((((*jet_trk_pt_var)[j_i])).at(j_t)/1000.) > lead_trk_pT ){
                  lead_trk_pT   = ((((*jet_trk_pt_var)[j_i])).at(j_t)/1000.);
                  lead_trk_eta  = (((*jet_trk_eta_var)[j_i]).at(j_t));
                  lead_trk_d0   = (((*jet_trk_d0_var)[j_i]).at(j_t));
                  lead_trk_z0   = (((*jet_trk_z0_var)[j_i]).at(j_t));
                }
              }
            }
            ntrk_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, trk_count);
            trk_pt_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv,lead_trk_pT);
            trk_eta_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, lead_trk_eta);
            trk_d0_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, lead_trk_d0);
            trk_z0_ch.fillSys_wCorl(had_label, mc_w, sysStr, isCorlEv, lead_trk_z0);
#endif
          }
        }
      }
      cout<<endl<<"--> "<<sysStr<<" done at: "<<double(clock())/CLOCKS_PER_SEC<<" <--"<<endl;
    }
  }

  if(!isMultiJob){
    //produce ratios
    mv2c10_ch.makeVar_h();
    mv2c10wp_ch.makeVar_h();
    pt_ch.makeVar_h();
    eta_ch.makeVar_h();
    etapt_ch.makeVar_h();
    cenEta_pt_mv2c10wp_ch.makeVar_h();
    frwEta_pt_mv2c10wp_ch.makeVar_h();
    ptFine_ch.makeVar_h();
#if _MORE_H > 0
    etaptRebin_ch.makeVar_h();
    mv2c10Flip_ch.makeVar_h();
    mv2c10wpFlip_ch.makeVar_h();
    cenEta_pt_mv2c10wpFlip_ch.makeVar_h();
    frwEta_pt_mv2c10wpFlip_ch.makeVar_h();
#endif
#if _MORE_H > 1
    ip3d_llr_ch.makeVar_h();
    sv1_llr_ch.makeVar_h();
    jf_llr_ch.makeVar_h();
#endif
#if _DOTRK == 1
    ntrk_ch.makeVar_h();
    trk_pt_ch.makeVar_h();
    trk_eta_ch.makeVar_h();
    trk_d0_ch.makeVar_h();
    trk_z0_ch.makeVar_h();
#endif
  }
  
  //store everything
  mv2c10_ch.write_h(&outF, !isMultiJob);
  mv2c10wp_ch.write_h(&outF, !isMultiJob);
  cenEta_pt_mv2c10wp_ch.write_h(&outF, !isMultiJob);
  frwEta_pt_mv2c10wp_ch.write_h(&outF, !isMultiJob);
  pt_ch.write_h(&outF, !isMultiJob);
  eta_ch.write_h(&outF, !isMultiJob);
  etapt_ch.write_h(&outF, !isMultiJob);
  ptFine_ch.write_h(&outF, !isMultiJob);
#if _MORE_H > 0
  etaptRebin_ch.write_h(&outF, !isMultiJob);
  mv2c10Flip_ch.write_h(&outF, !isMultiJob);
  mv2c10wpFlip_ch.write_h(&outF, !isMultiJob);
  cenEta_pt_mv2c10wpFlip_ch.write_h(&outF, !isMultiJob);
  frwEta_pt_mv2c10wpFlip_ch.write_h(&outF, !isMultiJob);
#endif
#if _MORE_H > 1
  ip3d_llr_ch.write_h(&outF, !isMultiJob);
  sv1_llr_ch.write_h(&outF, !isMultiJob);
  jf_llr_ch.write_h(&outF, !isMultiJob);
#endif
#if _DOTRK == 1
  ntrk_ch.write_h(&outF, !isMultiJob);
  trk_pt_ch.write_h(&outF, !isMultiJob);
  trk_eta_ch.write_h(&outF, !isMultiJob);
  trk_d0_ch.write_h(&outF, !isMultiJob);
  trk_z0_ch.write_h(&outF, !isMultiJob);
#endif
  
  outF.Close();
}
