#include "TH1.h"
#include "TString.h"
#include "TMath.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TClass.h"
#include "TChain.h"

#include <unistd.h>
#include <dirent.h>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <string>
#include <iostream>
#include <map>

using namespace std;


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


void MakeDijetMCWeight( string filename, bool isPy8 = true, string filesDir = "/eos/atlas/atlascerngroupdisk/phys-sm/VHF/bTag_studies/v02-00/dijet_py8/"){

  const bool debug = false;
  TFile* infile = new TFile( (filename+".root").c_str(), "recreate" );
  
  TH1D* h_xSection = new TH1D("hmc_xSection", ";sample;cross-section", 13, -0.5, 12.5);
  TH1D* h_filterEff = new TH1D("hmc_filterEff", ";sample;filter efficiency", 13, -0.5, 12.5);
  TH1D* h_sliceScale = new TH1D("hmc_sliceScale", ";sample; sample normalization yield", 13, -0.5, 12.5);

  for( int i=1; i<=h_xSection->GetXaxis()->GetNbins(); i++){
    h_xSection->GetXaxis()->SetBinLabel( i, TString::Format("JZ%iW", i-1));
    h_filterEff->GetXaxis()->SetBinLabel( i, TString::Format("JZ%iW", i-1));
    h_sliceScale->GetXaxis()->SetBinLabel( i, TString::Format("JZ%iW", i-1));
  }

  if(isPy8){ // https://twiki.cern.ch/twiki/bin/view/AtlasProtected/JetEtmissMC15#Pythia8_dijet
    h_xSection->SetBinContent(1, 7.8420 * TMath::Power(10, 7) );
    h_xSection->SetBinContent(2, 7.8420 * TMath::Power(10, 7) );
    h_xSection->SetBinContent(3, 2.4334 * TMath::Power(10, 6) );
    h_xSection->SetBinContent(4, 2.6454 * TMath::Power(10, 4) );
    h_xSection->SetBinContent(5, 2.5464 * TMath::Power(10, 2) );
    h_xSection->SetBinContent(6, 4.5536 * TMath::Power(10, 0) );
    h_xSection->SetBinContent(7, 2.5752 * TMath::Power(10, -1) );
    h_xSection->SetBinContent(8, 1.6214 * TMath::Power(10, -2) );
    h_xSection->SetBinContent(9, 6.2505 * TMath::Power(10, -4) );
    h_xSection->SetBinContent(10,1.9640 * TMath::Power(10, -5) );
    h_xSection->SetBinContent(11,1.1961 * TMath::Power(10, -6) );
    h_xSection->SetBinContent(12,4.2260 * TMath::Power(10, -8) );
    h_xSection->SetBinContent(13,1.0370 * TMath::Power(10, -9) );
    
    h_filterEff->SetBinContent(1, 1.0240 * TMath::Power(10, 0) );
    h_filterEff->SetBinContent(2, 6.7198 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(3, 3.3264 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(4, 3.1953 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(5, 5.3009 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(6, 9.2325 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(7, 9.4016 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(8, 3.9282 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(9, 1.0162 * TMath::Power(10, -2) );
    h_filterEff->SetBinContent(10,1.2054 * TMath::Power(10, -2) );
    h_filterEff->SetBinContent(11,5.8935 * TMath::Power(10, -3) );
    h_filterEff->SetBinContent(12,2.7015 * TMath::Power(10, -3) );
    h_filterEff->SetBinContent(13,4.2502 * TMath::Power(10, -4) );
  } else { //Herwig
    h_xSection->SetBinContent(1,  43.569 * TMath::Power(10, 0) );
    h_xSection->SetBinContent(2,  46.655 * TMath::Power(10, 0) );
    h_xSection->SetBinContent(3,  1.6777 * TMath::Power(10, 6) );
    h_xSection->SetBinContent(4,  1.8118 * TMath::Power(10, 4) );
    h_xSection->SetBinContent(5,  1.6755 * TMath::Power(10, 2) );
    h_xSection->SetBinContent(6,  2.8717 * TMath::Power(10, 0) );
    h_xSection->SetBinContent(7,  1.5757 * TMath::Power(10, -1) );
    h_xSection->SetBinContent(8,  9.7880 * TMath::Power(10, -3) );
    h_xSection->SetBinContent(9,  3.8667 * TMath::Power(10, -4) );
    h_xSection->SetBinContent(10, 1.3056 * TMath::Power(10, -5) );
    h_xSection->SetBinContent(11, 8.4140 * TMath::Power(10, -7) );
    h_xSection->SetBinContent(12, 3.0273 * TMath::Power(10, -8) );
    h_xSection->SetBinContent(13, 7.0212 * TMath::Power(10, -10) );

    h_filterEff->SetBinContent(1, 9.8374 * TMath::Power(10, -1) );
    h_filterEff->SetBinContent(2, 3.6905 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(3, 2.6121 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(4, 2.3539 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(5, 3.9276 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(6, 7.0652 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(7, 7.5887 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(8, 3.4609 * TMath::Power(10, -4) );
    h_filterEff->SetBinContent(9, 8.9708 * TMath::Power(10, -3) );
    h_filterEff->SetBinContent(10,1.0164 * TMath::Power(10, -2) );
    h_filterEff->SetBinContent(11,4.5845 * TMath::Power(10, -3) );
    h_filterEff->SetBinContent(12,1.8749 * TMath::Power(10, -3) );
    h_filterEff->SetBinContent(13,2.8761 * TMath::Power(10, -4) );
  }

  //Light Flavor Studies
  std::cout << "VHF Light Flavor Studies" << std::endl;
  TChain stdFile("bTag_AntiKt4EMTopoJets");
  ScanDir ( filesDir.c_str(), &stdFile, "Nominal");

  TTreeReader reader(&stdFile);
  TTreeReaderValue< int                   >  eventnb (reader, "eventnb");
  TTreeReaderValue< int                   >  runnb   (reader, "runnb");
  TTreeReaderValue< int                   >  mcchan  (reader, "mcchan");
  TTreeReaderValue< float                 >  mcwg    (reader, "mcwg"); 

  map<int, double> JZxW_yield_m;
  map<int, double> JZxW_nEvts_m;
  double mc_w;  int mc_ch;

  long int nEvent =0, nEventPrint = 100000;
  long int totEntries =  reader.GetEntries(true);
  for( nEvent=0; nEvent<totEntries; ++nEvent) {
    if ((nEvent%nEventPrint)==00){ cout<< "nEvent=" <<nEvent <<endl;
      for (auto m : JZxW_yield_m) cout<<"DSID: "<<m.first<<": yield "<<m.second<<", nEvts "<<JZxW_nEvts_m[m.first]<<endl;
    }
    reader.SetEntry(nEvent);
    mc_w = *mcwg; 
    mc_ch = *mcchan; 
    JZxW_yield_m[mc_ch] += mc_w;
    JZxW_nEvts_m[mc_ch] += 1;
  }

  int startDSID = 361020;//Pythia8: from 361020 to 361032
  if(!isPy8) startDSID = 426040;//Herwig: from 4260040 to 426052
  for(int i=0; i<13; ++i){
    if(JZxW_yield_m[startDSID+i] != 0)
      h_sliceScale->SetBinContent(i+1, 1/JZxW_yield_m[startDSID+i] );
  }

  infile->cd();
  TH1D* h_normScale = (TH1D*)h_sliceScale->Clone("hmc_normScale");
  h_normScale->GetYaxis()->SetTitle("sample normalization scale-factor");
  h_normScale->Multiply(h_xSection);
  h_normScale->Multiply(h_filterEff);

  h_xSection->Write();
  h_filterEff->Write();
  h_sliceScale->Write();
  h_normScale->Write();
  infile->Close();

  return;
}
