#include "TFile.h"
#include "TH1.h"
#include "TKey.h"
#include "TObject.h"
#include <string>

void merge_f(string fInBase, string fInAdd1, string fInAdd2 = "", string fInAdd3 = ""){
  TH1::SetDefaultSumw2(kTRUE);
  TH1::AddDirectory(kFALSE);

  TFile fout(fInBase.c_str(),"create");//"effPlot_FixedCutBEff_60_dijetPy8_32M_new.root", "recreate");
  TFile *f1 = new TFile(fInAdd1.c_str()); //"effPlot_FixedCutBEff_60_dijetPy8_32M_1h.root");

  TFile* f2 = nullptr;
  if(fInAdd2 != "") f2 = new TFile(fInAdd2.c_str());
  TFile* f3 = nullptr;
  if(fInAdd3 != "") f3 = new TFile(fInAdd3.c_str());

  TObject *obj;
  TKey *key;
  
  TIter next1( f1->GetListOfKeys());
  while ((key = (TKey *) next1())) {
    obj = f1->Get(key->GetName()); // copy object to memory
    string name = key->GetName();
    if(name.find("_LF") != string::npos ){
      fout.cd();
      obj->Write();
    }
  }

  if(fInAdd2 != "") {  
    TIter next2( f2->GetListOfKeys());
    while ((key = (TKey *) next2())) {
      obj = f2->Get(key->GetName()); // copy object to memory
      string name = key->GetName();    
      if( name.find("_STRANGEHAD_RATE") != string::npos ) {
        fout.cd();
        obj->Write();
      }
    }
  }
  
  if(fInAdd3 != ""){
    TIter next3( f3->GetListOfKeys());
    while ((key = (TKey *) next3())) {
      obj = f3->Get(key->GetName()); // copy object to memory
      string name = key->GetName();
      if( name.find("_RATETAIL") != string::npos ||
          name.find("_RATE_SYS") != string::npos ||
          name.find("_BIAS") != string::npos ||
          name.find("_EFF") != string::npos ) continue;
    
      if( name.find("_CONV_HADINT_RATE") != string::npos ||
          name.find("_STRANGEHAD_RATE") != string::npos ||
          name.find("_MEAS_UP") != string::npos ||
          name.find("_MEAS_DOWN") != string::npos ||
          name.find("_MEAS_Zmm") != string::npos ||
          name.find("_FAKE_RATE") != string::npos ||
          name.find("_TAIL_RATE") != string::npos ||
          name.find("_D0Z0Corl") != string::npos ) {
        fout.cd();
        obj->Write();
      }
    }
  }

  fout.Close();
}
