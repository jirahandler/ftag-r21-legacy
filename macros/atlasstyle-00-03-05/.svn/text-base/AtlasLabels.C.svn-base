#include "AtlasLabels.h"

#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TPad.h"
#include "TMarker.h"


void ATLASLabel(Double_t x,Double_t y,const char* text,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"ATLAS");
  if (text) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}


void ATLASLabelLumi(Double_t x,Double_t y,const char* text,Color_t color)
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double dely = 0.027*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"ATLAS");

  TLatex p;
  p.SetNDC();
  p.SetTextFont(42);
  p.SetTextSize(0.03);
  p.SetTextColor(color);
  p.DrawLatex(x,y-dely-dely/2," #int#it{Ldt} = 20.34 fb^{-1}");

  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x,y-dely-dely-dely,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}

void ATLASLabelLumiSqrt(Double_t x,Double_t y,const char* text,Color_t color)
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double dely = 0.027*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"ATLAS");

  TLatex p;
  p.SetNDC();
  p.SetTextFont(42);
  p.SetTextSize(0.03);
  p.SetTextColor(color);
  p.DrawLatex(x,y-dely-dely/2," #int#it{Ldt} = 20.34 fb^{-1}, #sqrt{s} = 8 TeV");

  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x,y-dely-dely-dely,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}

void ATLASLabelSqrt(Double_t x,Double_t y,const char* text,Color_t color)
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  //l.SetTextSize(0.068); //to be used for canvas with ratio pad
  //l.SetTextSize(0.068); //to be used for canvas with ratio pad
  
  l.SetTextSize(0.045); //to be used for single canvas
  l.SetTextColor(color);

  //double dely = 0.035*696*gPad->GetWh()/(472*gPad->GetWw());
  double dely = 0.048*696*gPad->GetWh()/(472*gPad->GetWw());
  
  //double delx = 0.280*696*gPad->GetWh()/(472*gPad->GetWw());
  //double delx = 0.220*696*gPad->GetWh()/(472*gPad->GetWw());
  double delx = 0.230*696*gPad->GetWh()/(472*gPad->GetWw());
  //double delx = 0.180*696*gPad->GetWh()/(472*gPad->GetWw());

  //l.DrawLatex(x,y,"ATLAS Internal");
  l.DrawLatex(x,y,"ATLAS");

  TLatex p;
  p.SetNDC();
  //p.SetTextFont(72);
  p.SetTextFont(42);
  p.SetTextSize(0.045); //single canvas
  //p.SetTextSize(0.045);
  //p.SetTextSize(0.068); //ratio canvas
  //p.SetTextSize(0.03);
  p.SetTextColor(color);
  //p.DrawLatex(x,y-dely-dely/2," #sqrt{s} = 13 TeV");
  p.DrawLatex(x+delx/5,y-dely-dely/2," #sqrt{s} = 13 TeV"); //to draw the text in the middle

  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextSize(0.045); //single canvas
    //p.SetTextSize(0.068); //ratio canvas
    p.SetTextColor(color);
    //p.DrawLatex(x,y-dely-dely-dely,text);
    p.DrawLatex(x+delx/2,y,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}

void ATLASLabelSimulation(Double_t x,Double_t y,const char* text,Color_t color)
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double dely = 0.025*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"ATLAS Simulation");
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x,y-dely,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}


void ATLASLabelOld(Double_t x,Double_t y,bool Preliminary,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
  if (Preliminary) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+0.115,y,"Preliminary");
  }
}



void ATLASVersion(const char* version,Double_t x,Double_t y,Color_t color) 
{

  if (version) {
    char versionString[100];
    sprintf(versionString,"Version %s",version);
    TLatex l;
    l.SetTextAlign(22); 
    l.SetTextSize(0.04); 
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextColor(color);
    l.DrawLatex(x,y,versionString);
  }
}

