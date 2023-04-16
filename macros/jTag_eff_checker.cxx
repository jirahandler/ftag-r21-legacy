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
#include <time.h>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
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
#include "TGraphAsymmErrors.h"

#include "atlasstyle-00-03-05/AtlasLabels.C"

#define _MOREH 0
#define _DOTRK 0

using namespace std;
const bool gDoJVT = true;
const int gCheckPU = 2;//0 = no check,  1 = check for DR<0.6, 2 = check for DR<0.3
const bool gUseMCWeights = false;

//Scan files in a Directory
void ScanDir (string file1, TChain* myT_1, string sysName="") {
  TList* hList = new TList();
  cout << "Input is a directory: going fancy: " << endl;
  DIR*     dir;
  dirent*  pdir;
  dir = opendir( file1.c_str() );     // open current directory
  while ((pdir = readdir(dir)))  {     // loop on directory contents
    string foldName=pdir->d_name;
    if (foldName == "..") continue;
    if (!sysName.empty() && sysName != foldName) continue;
    cout << pdir->d_name << endl;
    DIR*     dir2;
    dirent*  pdir2;
    dir2 = opendir( (file1+"/"+foldName).c_str() );     // open current directory
    //dir2 = opendir( (file1+"/"+foldName+").c_str() );     // open current directory
    int nFiles = 0;
    while ((pdir2 = readdir(dir2)))  {
      string fName=pdir2->d_name;
      if (fName.find("root")==string::npos) continue;
      ++nFiles;
      //      cout << fName << endl;
      myT_1->Add( (file1+"/"+foldName+"/"+fName).c_str() );
      //myT_1->Add( (file1+"/"+foldName+"/*/"+fName).c_str() );
    }
    cout <<"NFiles = "<< nFiles << endl;
  }  
}

enum flavId{ LF_id = 0, C_id = 4, B_id =5 };

//+++++++ TESTED in 
// eosmount eosLink
// asetup 20.7.8.7, here
//++++++++++++

  //mv2c10 working points listed here:
  // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/BTaggingBenchmarks#MV2c10_tagger_added_on_11th_May
// FixedCutBEff_30 	0.997715 	29.99 	99.95 	1406.40 	4891.95 	51290.94
// FixedCutBEff_50 	0.976933 	50.05 	99.62 	106.41 	614.79 	5701.28
// FixedCutBEff_60 	0.934906 	60.03 	99.00 	34.54 	183.98 	1538.78
// FixedCutBEff_70 	0.824427 	69.97 	97.46 	12.17 	54.72 	381.32
// FixedCutBEff_77 	0.645925 	76.97 	95.17 	6.21 	22.04 	134.34
// FixedCutBEff_85 	0.175847 	84.95 	89.66 	3.10 	8.17 	33.53 

void jTag_eff_checker(string outFileID = "effPlot",  
                      string filesDir =  "/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/egraham_v01-03/", 
                      const int nEventMax =  300000000, const bool doFlipTag = false, 
                      const double maxJetPt = -1, const string addDir = "Nominal"){

  TH1::SetDefaultSumw2(kTRUE);
  TH1::AddDirectory(kFALSE);
  const int nEventPrint = 100000;

  clock_t startTime = clock();
  cout<<endl<<"--> Start filling hists at: "<<double(startTime)/CLOCKS_PER_SEC<<" <--"<<endl;

  TFile outF((outFileID+".root").c_str(), "recreate");

  //Light Flavor Studies
  std::cout << "VHF Light Flavor Studies" << std::endl;
  TChain stdFile("bTag_AntiKt4EMTopoJets");
  ScanDir ( filesDir.c_str(), &stdFile, addDir);

  double ptbinsFine[] = {20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 120, 150, 200, 300, 500, 750, 1000, 3000}; //to be used for both 1D and 2D
  int   n_jetpt_binsFine = sizeof(ptbinsFine)/sizeof(double) - 1;  

  double ptbins[] = {20, 30, 60, 100, 150, 200, 500, 1000, 3000}; //to be used for both 1D and 2D
  int   n_jetpt_bins = sizeof(ptbins)/sizeof(double) - 1;  

  //Corresponding to Matthias' negative-LF calibration
  double ptbinsRebin[] = {20, 60, 100, 200, 300, 500, 750, 1000, 3000}; //to be used for both 1D and 2D
  int   n_jetptRebin_bins = sizeof(ptbinsRebin)/sizeof(double) - 1;  
  
  double etabins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.5}; //to be used for both 1D and 2D
  int   n_jeteta_bins = sizeof(etabins)/sizeof(double) - 1;

  double etabinsRebin[] = {0.0,1.2,2.5}; //to be used for both 1D and 2D
  int   n_jetetaRebin_bins = sizeof(etabinsRebin)/sizeof(double) - 1;

  int n_jetpt_binsFine_2D  =  n_jetpt_binsFine;
  int n_jetpt_bins_2D  =  n_jetpt_bins;
  int n_jeteta_bins_2D =  n_jeteta_bins;
  int n_jetptRebin_bins_2D  =  n_jetptRebin_bins;
  int n_jetetaRebin_bins_2D =  n_jetetaRebin_bins;

  //  float mv2c10_bins[] = {-1.0, 0.175847, 0.645925, 0.824427, 0.934906, 0.976933, 0.997715, 1};
  float mv2c10_bins[] = {-1.0, 0.175847, 0.645925, 0.824427, 0.934906, 0.976933, 1};//No 30% WP
  vector<string> wpStr_v;
  if(!doFlipTag) wpStr_v = {"noTag", "wp85", "wp77", "wp70", "wp60", "wp50"};
  else wpStr_v = {"noTag", "wpNeg85", "wpNeg77", "wpNeg70", "wpNeg60", "wpNeg50"};
  int n_mv2c10_bins = sizeof(mv2c10_bins)/sizeof(float) - 1;

  vector<string> typeStr_v = {"other", "conv", "matInt", "sHad", "fake"};
  map<int, string> flav_m;
  flav_m[LF_id] = "LF";
  flav_m[B_id] = "B";
  flav_m[C_id] = "C";


  //map of totals: Flav(B,C,LF)->WP(All, 85, ...)->Hist(pt, mv2c10, ...)
  map<string,  map<string, map<string,  TH1F* > > > total_h_m;
  map<string,  map<string, map<string,  TH2F* > > > total_2h_m;

  //map of: Flav(B,C,LF)->WP(NoTag, 85, ...)->TYPE(Any,Conv,HadInt, s-had,...)->Hist(pt, mv2c10, ...)
  map<string, map< string, map<string, map<string,  TH1F* > > > > classes_h_m;
  map<string, map< string, map<string, map<string,  TH2F* > > > > classes_2h_m;

  for(auto& flav_i : flav_m ){
    string flav= flav_i.second;
    for(uint wp_i=0; wp_i<wpStr_v.size(); ++wp_i){
      string wp = wpStr_v[wp_i];
      
      total_h_m[flav][wp]["ptFine"] = new TH1F(Form("ptFine_%s_%s",flav.c_str(),wp.c_str()),
                                               Form("ptFine_%s_%s",flav.c_str(),wp.c_str()), 
                                               n_jetpt_binsFine, ptbinsFine);
      
      total_h_m[flav][wp]["pt"] = new TH1F(Form("pt_%s_%s",flav.c_str(),wp.c_str()),
                                           Form("pt_%s_%s",flav.c_str(),wp.c_str()), 
                                           n_jetpt_bins, ptbins);
      
      total_h_m[flav][wp]["eta"] = new TH1F(Form("eta_%s_%s",flav.c_str(),wp.c_str()),
                                            Form("eta_%s_%s",flav.c_str(),wp.c_str()), 
                                            n_jeteta_bins, etabins);
      
      total_h_m[flav][wp]["mv2c10"] = new TH1F(Form("mv2c10_%s_%s",flav.c_str(),wp.c_str()),
                                               Form("mv2c10_%s_%s",flav.c_str(),wp.c_str()), 202, -1.01, 1.01);

      total_2h_m[flav][wp]["etaptFine"] = new TH2F(Form("etaptFine_%s_%s",flav.c_str(),wp.c_str()),
                                                   Form("etaptFine_%s_%s",flav.c_str(),wp.c_str()),  
                                                   n_jetpt_binsFine_2D, ptbinsFine, n_jeteta_bins_2D, etabins);
      
      total_2h_m[flav][wp]["etapt"] = new TH2F(Form("etapt_%s_%s",flav.c_str(),wp.c_str()),
                                               Form("etapt_%s_%s",flav.c_str(),wp.c_str()),  
                                               n_jetpt_bins_2D, ptbins, n_jetetaRebin_bins_2D, etabinsRebin);
      
      total_2h_m[flav][wp]["etaptRebin"] = new TH2F(Form("etaptRebin_%s_%s",flav.c_str(),wp.c_str()),
                                                    Form("etaptRebin_%s_%s",flav.c_str(),wp.c_str()),  
                                                    n_jetptRebin_bins_2D, ptbinsRebin, n_jetetaRebin_bins_2D, etabinsRebin);
      
      for(uint t_i=0; t_i<typeStr_v.size(); ++t_i){
        string type = typeStr_v[t_i];        

        //      cout<<Form("pt_%s_%s_%s",flav.c_str(),type.c_str(),wp.c_str())<<endl;
        classes_h_m[flav][wp][type]["ptFine"] = new TH1F(Form("ptFine_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                         Form("ptFine_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()), 
                                                         n_jetpt_binsFine, ptbinsFine);

        classes_h_m[flav][wp][type]["pt"] = new TH1F(Form("pt_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                 Form("pt_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()), 
                                                 n_jetpt_bins, ptbins);

        classes_h_m[flav][wp][type]["eta"] = new TH1F(Form("eta_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                  Form("eta_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()), 
                                                  n_jeteta_bins, etabins);

        classes_h_m[flav][wp][type]["mv2c10"] = new TH1F(Form("mv2c10_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                         Form("mv2c10_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()), 202, -1.01, 1.01);

        classes_2h_m[flav][wp][type]["etaptFine"] = new TH2F(Form("etaptFine_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                             Form("etapt_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),  
                                                             n_jetpt_binsFine_2D, ptbinsFine, n_jeteta_bins_2D, etabins);

        classes_2h_m[flav][wp][type]["etapt"] = new TH2F(Form("etaptRebin_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                         Form("etaptRebin_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),  
                                                         n_jetpt_bins_2D, ptbins, n_jetetaRebin_bins_2D, etabinsRebin);

        classes_2h_m[flav][wp][type]["etaptRebin"] = new TH2F(Form("etaptRebin_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),
                                                            Form("etaptRebin_%s_%s_%s",flav.c_str(),wp.c_str(),type.c_str()),  
                                                            n_jetptRebin_bins_2D, ptbinsRebin, n_jetetaRebin_bins_2D, etabinsRebin);
      }
    }
  }  

  TTreeReader reader(&stdFile);
  TTreeReaderValue< int                   >  eventnb         (reader, "eventnb"); 
  TTreeReaderValue< int                   >  njets           (reader, "njets"); 
  TTreeReaderValue< float                 >  mcwg            (reader, "mcwg"); 
  TTreeReaderValue< vector<float>         >  jet_pt          (reader, "jet_pt"); 
  TTreeReaderValue< vector<float>         >  jet_eta         (reader, "jet_eta");
  TTreeReaderValue< vector<float>         >  jet_JVT         (reader, "jet_JVT");
  TTreeReaderValue< vector<int>           >  jet_isPU        (reader, "jet_isPU"); 
  TTreeReaderValue< vector<int>           >  jet_truthMatch  (reader, "jet_truthMatch"); 
  TTreeReaderValue< vector<int>           >  jet_LFCalibType (reader, "jet_LFCalibType");
  TTreeReaderValue< vector<int>           >  jet_aliveAfterOR(reader, "jet_aliveAfterOR");
  TTreeReaderValue< vector<int>           >  jet_aliveAfterORmu(reader, "jet_aliveAfterORmu");
  TTreeReaderValue< vector<int>           >  jet_LabDr_HadF  (reader, "jet_LabDr_HadF");  
  TTreeReaderValue< vector<float>         >  jet_mv2c10      (reader, "jet_mv2c10");
  TTreeReaderValue< vector<float>         >  jet_mv2c10Flip  (reader, "jet_mv2c20Flip");
#if _MORE_H == 1
  TTreeReaderValue< vector<int>           >  jet_truthflav   (reader, "jet_truthflav");
#endif
#if _DOTRK == 1
  TTreeReaderValue< vector<int>           >  jet_btag_ntrk   (reader, "jet_btag_ntrk");
  TTreeReaderValue< vector <vector<float> > > jet_trk_truthMatchProbability(reader, "jet_trk_truthMatchProbability");
#endif

  int nEvent=0;
  const float track_pt_cut = 0.7; //cut on pt leading track (higher than 700 MeV)

  float j_pt, j_eta, mv2c10, mc_weight;
  int j_i, had_label;
  int LFCalibType, isPU;
  string flavStr, wpStr;

  ///// EVENT LOOP ////////////  
  while (reader.Next() && nEvent < nEventMax) {
    ++nEvent;
    if ((nEvent%nEventPrint)==0) cout<< "nEvent=" <<nEvent <<endl;
    mc_weight = gUseMCWeights ? (*mcwg) : 1.; 
    for(j_i  =0; j_i < *njets; ++j_i){
      if(gCheckPU >0 && (*jet_isPU)[j_i]) continue;
      if(gCheckPU >1 && !(*jet_truthMatch)[j_i]) continue;

      //Needed when analyzing samples with truth-ele or mu, which are removed
      if(!(*jet_aliveAfterOR)[j_i]) continue;
      if(!(*jet_aliveAfterORmu)[j_i]) continue;

      j_pt = (*jet_pt)[j_i]/1000.;
      j_eta = fabs( (*jet_eta)[j_i]);

      if(maxJetPt > -1 && j_pt > maxJetPt) continue;

      if(j_pt<20. || j_eta>2.5) continue;//Standard calibration cuts
      //Jet cleaning probably should be applied only to not-bad jets but the fraction is minimal on MC
      if(gDoJVT && j_pt<60. && j_eta<2.4 && (*jet_JVT)[j_i]<0.59) continue;

      //Check jet flavour
      had_label = (*jet_LabDr_HadF)[j_i];
      if( !(had_label == LF_id || had_label == B_id || had_label == C_id )) continue;
      flavStr = flav_m[had_label];

      //Check jet mv2c10 value and corresponding b-tag working point
      mv2c10 = (*jet_mv2c10)[j_i];
      if(doFlipTag) mv2c10 = (*jet_mv2c10Flip)[j_i];//Use the flipper tagger for these tests!

      LFCalibType = (*jet_LFCalibType)[j_i];

      float lowest_truthMatchProbability  = 999.;
#if _DOTRK == 1
      for(int j_t  =0; j_t < ((*jet_btag_ntrk)[j_i]); ++j_t){
        if( (*jet_trk_truthMatchProbability)[j_i][j_t] < lowest_truthMatchProbability ){
          lowest_truthMatchProbability = (*jet_trk_truthMatchProbability)[j_i][j_t];
        }
      }
      if(lowest_truthMatchProbability <1.) LFCalibType |= 0x20;
#endif

      for(uint wp_i=0; wp_i<wpStr_v.size(); ++wp_i){
        if(mv2c10 > mv2c10_bins[wp_i] ){
          wpStr = wpStr_v[wp_i]; 
          
          // cout<<"j_pt = "<<j_pt<<" , had_label = "<<had_label<<" , mv2c10 = "<<mv2c10<<" , LFCalibType = "<<LFCalibType<<endl;
          total_h_m[flavStr][wpStr]["ptFine"]->Fill(j_pt, mc_weight);
          total_h_m[flavStr][wpStr]["pt"]->Fill(j_pt, mc_weight);
          total_h_m[flavStr][wpStr]["eta"]->Fill(j_eta, mc_weight);
          total_h_m[flavStr][wpStr]["mv2c10"]->Fill(mv2c10, mc_weight);
          total_2h_m[flavStr][wpStr]["etaptFine"]->Fill(j_pt,j_eta,mc_weight);
          total_2h_m[flavStr][wpStr]["etapt"]->Fill(j_pt,j_eta, mc_weight);
          total_2h_m[flavStr][wpStr]["etaptRebin"]->Fill(j_pt, j_eta, mc_weight);


          if( (LFCalibType & (0x1 | 0x2)) == (0x1 | 0x2) ) {
            classes_h_m[flavStr][wpStr]["sHad"]["ptFine"]->Fill(j_pt, mc_weight);          
            classes_h_m[flavStr][wpStr]["sHad"]["pt"]->Fill(j_pt, mc_weight);
            classes_h_m[flavStr][wpStr]["sHad"]["eta"]->Fill(j_eta, mc_weight);
            classes_h_m[flavStr][wpStr]["sHad"]["mv2c10"]->Fill(mv2c10, mc_weight);
            classes_2h_m[flavStr][wpStr]["sHad"]["etaptFine"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["sHad"]["etapt"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["sHad"]["etaptRebin"]->Fill(j_pt, j_eta, mc_weight);
          }

          if( (LFCalibType & 0x4) == 0x4 ) {
            classes_h_m[flavStr][wpStr]["conv"]["ptFine"]->Fill(j_pt, mc_weight);          
            classes_h_m[flavStr][wpStr]["conv"]["pt"]->Fill(j_pt, mc_weight);
            classes_h_m[flavStr][wpStr]["conv"]["eta"]->Fill(j_eta, mc_weight);
            classes_h_m[flavStr][wpStr]["conv"]["mv2c10"]->Fill(mv2c10, mc_weight);
            classes_2h_m[flavStr][wpStr]["conv"]["etaptFine"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["conv"]["etapt"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["conv"]["etaptRebin"]->Fill(j_pt, j_eta, mc_weight);
          }

          if( (LFCalibType & 0x8) == 0x8 ) {
            classes_h_m[flavStr][wpStr]["matInt"]["ptFine"]->Fill(j_pt, mc_weight);          
            classes_h_m[flavStr][wpStr]["matInt"]["pt"]->Fill(j_pt, mc_weight);
            classes_h_m[flavStr][wpStr]["matInt"]["eta"]->Fill(j_eta, mc_weight);
            classes_h_m[flavStr][wpStr]["matInt"]["mv2c10"]->Fill(mv2c10, mc_weight);
            classes_2h_m[flavStr][wpStr]["matInt"]["etaptFine"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["matInt"]["etapt"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["matInt"]["etaptRebin"]->Fill(j_pt, j_eta, mc_weight);
          }

          if( (LFCalibType & 0x10) == 0x10 ) {
            classes_h_m[flavStr][wpStr]["fake"]["ptFine"]->Fill(j_pt, mc_weight);          
            classes_h_m[flavStr][wpStr]["fake"]["pt"]->Fill(j_pt, mc_weight);
            classes_h_m[flavStr][wpStr]["fake"]["eta"]->Fill(j_eta, mc_weight);
            classes_h_m[flavStr][wpStr]["fake"]["mv2c10"]->Fill(mv2c10, mc_weight);
            classes_2h_m[flavStr][wpStr]["fake"]["etaptFine"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["fake"]["etapt"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["fake"]["etaptRebin"]->Fill(j_pt, j_eta, mc_weight);
          }

          if( LFCalibType == 0 ) {
            classes_h_m[flavStr][wpStr]["other"]["ptFine"]->Fill(j_pt, mc_weight);          
            classes_h_m[flavStr][wpStr]["other"]["pt"]->Fill(j_pt, mc_weight);
            classes_h_m[flavStr][wpStr]["other"]["eta"]->Fill(j_eta, mc_weight);
            classes_h_m[flavStr][wpStr]["other"]["mv2c10"]->Fill(mv2c10, mc_weight);
            classes_2h_m[flavStr][wpStr]["other"]["etaptFine"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["other"]["etapt"]->Fill(j_pt,j_eta, mc_weight);
            classes_2h_m[flavStr][wpStr]["other"]["etaptRebin"]->Fill(j_pt, j_eta, mc_weight);
          }

        }
      }
    }
  }

  cout<<endl<<"--> Nominal done at: "<<double(clock())/CLOCKS_PER_SEC<<" <--"<<endl;  
  
  //store everything
  outF.cd();
  for(auto& flav_i : flav_m ){
    string flav= flav_i.second;
    for(uint wp_i=0; wp_i<wpStr_v.size(); ++wp_i){
      string wp = wpStr_v[wp_i];

      //Initial inputs      
      total_h_m[flav][wp]["ptFine"]->Write();
      total_h_m[flav][wp]["pt"]->Write();
      total_h_m[flav][wp]["eta"]->Write();
      total_h_m[flav][wp]["mv2c10"]->Write();
      total_2h_m[flav][wp]["etaptFine"]->Write();
      total_2h_m[flav][wp]["etapt"]->Write();
      total_2h_m[flav][wp]["etaptRebin"]->Write();

      for(uint t_i=0; t_i<typeStr_v.size(); ++t_i){
        string type = typeStr_v[t_i];        

        //Initial inputs of different type
        classes_h_m[flav][wp][type]["ptFine"]->Write();
        classes_h_m[flav][wp][type]["pt"]->Write();
        classes_h_m[flav][wp][type]["eta"]->Write();
        classes_h_m[flav][wp][type]["mv2c10"]->Write();
        classes_2h_m[flav][wp][type]["etaptFine"]->Write();
        classes_2h_m[flav][wp][type]["etapt"]->Write();
        classes_2h_m[flav][wp][type]["etaptRebin"]->Write();
	if(wp=="noTag") continue; 

        //Calculate various kind of tagging efficiencies

        //eff w.r.t. total LF, C, B jets
        TH1F* effTmp_ptFine    = (TH1F*) classes_h_m[flav][wp][type]["ptFine"]->Clone(Form("%s_TOTeff",classes_h_m[flav][wp][type]["ptFine"]->GetName()) );
        TH1F* effTmp_pt        = (TH1F*) classes_h_m[flav][wp][type]["pt"]->Clone(Form("%s_TOTeff",classes_h_m[flav][wp][type]["pt"]->GetName()) );
        TH1F* effTmp_eta       = (TH1F*) classes_h_m[flav][wp][type]["eta"]->Clone(Form("%s_TOTeff",classes_h_m[flav][wp][type]["eta"]->GetName()) );
        TH1F* effTmp_mv2c10    = (TH1F*) classes_h_m[flav][wp][type]["mv2c10"]->Clone(Form("%s_TOTeff",classes_h_m[flav][wp][type]["mv2c10"]->GetName()) );
        TH2F* effTmp_etaptFine = (TH2F*) classes_2h_m[flav][wp][type]["etaptFine"]->Clone(Form("%s_TOTeff",classes_2h_m[flav][wp][type]["etaptFine"]->GetName()) );
        TH2F* effTmp_etapt     = (TH2F*) classes_2h_m[flav][wp][type]["etapt"]->Clone(Form("%s_TOTeff",classes_2h_m[flav][wp][type]["etapt"]->GetName()) );
        TH2F* effTmp_etaptRebin= (TH2F*) classes_2h_m[flav][wp][type]["etaptRebin"]->Clone(Form("%s_TOTeff",classes_2h_m[flav][wp][type]["etaptRebin"]->GetName()) );

        effTmp_ptFine    ->Divide(total_h_m[flav]["noTag"]["ptFine"]); effTmp_ptFine->Write();
        effTmp_pt        ->Divide(total_h_m[flav]["noTag"]["pt"]); effTmp_pt->Write();
        effTmp_eta       ->Divide(total_h_m[flav]["noTag"]["eta"]); effTmp_eta->Write();
        effTmp_mv2c10    ->Divide(total_h_m[flav]["noTag"]["mv2c10"]); effTmp_mv2c10->Write();
        effTmp_etaptFine ->Divide(total_2h_m[flav]["noTag"]["etaptFine"]); effTmp_etaptFine->Write();
        effTmp_etapt     ->Divide(total_2h_m[flav]["noTag"]["etapt"]); effTmp_etapt->Write();
        effTmp_etaptRebin->Divide(total_2h_m[flav]["noTag"]["etaptRebin"]); effTmp_etaptRebin->Write();

        //eff w.r.t. LF, C, B jets of given "type"
        effTmp_ptFine    = (TH1F*) classes_h_m[flav][wp][type]["ptFine"]->Clone(Form("%s_eff",classes_h_m[flav][wp][type]["ptFine"]->GetName()) );
        effTmp_pt        = (TH1F*) classes_h_m[flav][wp][type]["pt"]->Clone(Form("%s_eff",classes_h_m[flav][wp][type]["pt"]->GetName()) );
        effTmp_eta       = (TH1F*) classes_h_m[flav][wp][type]["eta"]->Clone(Form("%s_eff",classes_h_m[flav][wp][type]["eta"]->GetName()) );
        effTmp_mv2c10    = (TH1F*) classes_h_m[flav][wp][type]["mv2c10"]->Clone(Form("%s_eff",classes_h_m[flav][wp][type]["mv2c10"]->GetName()) );
        effTmp_etaptFine = (TH2F*) classes_2h_m[flav][wp][type]["etaptFine"]->Clone(Form("%s_eff",classes_2h_m[flav][wp][type]["etaptFine"]->GetName()) );
        effTmp_etapt     = (TH2F*) classes_2h_m[flav][wp][type]["etapt"]->Clone(Form("%s_eff",classes_2h_m[flav][wp][type]["etapt"]->GetName()) );
        effTmp_etaptRebin= (TH2F*) classes_2h_m[flav][wp][type]["etaptRebin"]->Clone(Form("%s_eff",classes_2h_m[flav][wp][type]["etaptRebin"]->GetName()) );


        effTmp_ptFine    ->Divide(classes_h_m[flav]["noTag"][type]["ptFine"]); effTmp_ptFine->Write();
        effTmp_pt        ->Divide(classes_h_m[flav]["noTag"][type]["pt"]); effTmp_pt->Write();
        effTmp_eta       ->Divide(classes_h_m[flav]["noTag"][type]["eta"]); effTmp_eta->Write();
        effTmp_mv2c10    ->Divide(classes_h_m[flav]["noTag"][type]["mv2c10"]); effTmp_mv2c10->Write();
        effTmp_etaptFine ->Divide(classes_2h_m[flav]["noTag"][type]["etaptFine"]); effTmp_etaptFine->Write();
        effTmp_etapt     ->Divide(classes_2h_m[flav]["noTag"][type]["etapt"]); effTmp_etapt->Write();
        effTmp_etaptRebin->Divide(classes_2h_m[flav]["noTag"][type]["etaptRebin"]); effTmp_etaptRebin->Write();

        //Fraction of given type and flavour for each WP
        classes_h_m[flav][wp][type]["ptFine"]->Divide(total_h_m[flav][wp]["ptFine"]);
        classes_h_m[flav][wp][type]["ptFine"]->Write( Form("%s_frac",classes_h_m[flav][wp][type]["ptFine"]->GetName()));
        classes_h_m[flav][wp][type]["pt"]->Divide(total_h_m[flav][wp]["pt"]);
        classes_h_m[flav][wp][type]["pt"]->Write( Form("%s_frac",classes_h_m[flav][wp][type]["pt"]->GetName()));
        classes_h_m[flav][wp][type]["eta"]->Divide(total_h_m[flav][wp]["eta"]);
        classes_h_m[flav][wp][type]["eta"]->Write( Form("%s_frac",classes_h_m[flav][wp][type]["eta"]->GetName()));
        classes_h_m[flav][wp][type]["mv2c10"]->Divide(total_h_m[flav][wp]["mv2c10"]);
        classes_h_m[flav][wp][type]["mv2c10"]->Write( Form("%s_frac",classes_h_m[flav][wp][type]["mv2c10"]->GetName()));
        classes_2h_m[flav][wp][type]["etaptFine"]->Divide(total_2h_m[flav][wp]["etaptFine"]);
        classes_2h_m[flav][wp][type]["etaptFine"]->Write( Form("%s_frac",classes_2h_m[flav][wp][type]["etaptFine"]->GetName()));
        classes_2h_m[flav][wp][type]["etapt"]->Divide(total_2h_m[flav][wp]["etapt"]);
        classes_2h_m[flav][wp][type]["etapt"]->Write( Form("%s_frac",classes_2h_m[flav][wp][type]["etapt"]->GetName()));
        classes_2h_m[flav][wp][type]["etaptRebin"]->Divide(total_2h_m[flav][wp]["etaptRebin"]);
        classes_2h_m[flav][wp][type]["etaptRebin"]->Write( Form("%s_frac",classes_2h_m[flav][wp][type]["etaptRebin"]->GetName()));
      }

      if(wp=="noTag") continue; 
      //Tag fractions out of the box for different variables
      total_h_m[flav][wp]["ptFine"]->Divide(total_h_m[flav]["noTag"]["ptFine"]);
      total_h_m[flav][wp]["ptFine"]->Write( Form("%s_frac",total_h_m[flav][wp]["ptFine"]->GetName()));
      total_h_m[flav][wp]["pt"]->Divide(total_h_m[flav]["noTag"]["pt"]);
      total_h_m[flav][wp]["pt"]->Write( Form("%s_frac",total_h_m[flav][wp]["pt"]->GetName()));
      total_h_m[flav][wp]["eta"]->Divide(total_h_m[flav]["noTag"]["eta"]);
      total_h_m[flav][wp]["eta"]->Write( Form("%s_frac",total_h_m[flav][wp]["eta"]->GetName()));
      total_h_m[flav][wp]["mv2c10"]->Divide(total_h_m[flav]["noTag"]["mv2c10"]);
      total_h_m[flav][wp]["mv2c10"]->Write( Form("%s_frac",total_h_m[flav][wp]["mv2c10"]->GetName()));
      total_2h_m[flav][wp]["etaptFine"]->Divide(total_2h_m[flav]["noTag"]["etaptFine"]);
      total_2h_m[flav][wp]["etaptFine"]->Write( Form("%s_frac",total_2h_m[flav][wp]["etaptFine"]->GetName()));
      total_2h_m[flav][wp]["etapt"]->Divide(total_2h_m[flav]["noTag"]["etapt"]);
      total_2h_m[flav][wp]["etapt"]->Write( Form("%s_frac",total_2h_m[flav][wp]["etapt"]->GetName()));
      total_2h_m[flav][wp]["etaptRebin"]->Divide(total_2h_m[flav]["noTag"]["etaptRebin"]);
      total_2h_m[flav][wp]["etaptRebin"]->Write( Form("%s_frac",total_2h_m[flav][wp]["etaptRebin"]->GetName()));
    }
  }
  outF.Close();
}
