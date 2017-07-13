#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TUnfoldDensity.h"

void atributes(TH1D *histo, TString xtitle="", Int_t COLOR = 1, TString ytitle="Fraction", Int_t style = 1){
  histo->ResetAttLine();
  histo->ResetAttFill();
  histo->ResetAttMarker();
  histo->SetTitle("");
  histo->SetMarkerStyle(kFullCircle);
  histo->SetMarkerSize(0.8);
  histo->SetLineWidth(4);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetXaxis()->SetLabelFont  (   42);
  histo->GetXaxis()->SetLabelOffset(0.015);
  histo->GetXaxis()->SetLabelSize  (0.050);
  histo->GetXaxis()->SetNdivisions (  505);
  histo->GetXaxis()->SetTitleFont  (   42);
  histo->GetXaxis()->SetTitleOffset(  1.4);
  histo->GetXaxis()->SetTitleSize  (0.035);

  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetYaxis()->SetLabelFont  (   42);
  histo->GetYaxis()->SetLabelOffset(0.015);
  histo->GetYaxis()->SetLabelSize  (0.050);
  histo->GetYaxis()->SetNdivisions (  505);
  histo->GetYaxis()->SetTitleFont  (   42);
  histo->GetYaxis()->SetTitleOffset(  1.2);
  histo->GetYaxis()->SetTitleSize  (0.035);
  histo->SetLineColor(COLOR);
  histo->SetMarkerStyle(kFullDotLarge);
  histo->SetLineStyle(style);
}
void makeSystHist(int nsel = 0, int whichDY = 1, TString theHistName = "Pt"){
  const int nBinPt = 57; Float_t xbinsPt[nBinPt+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                       20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                                                       30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
                                                       40, 50, 60, 70, 80, 90,100,125,150,175, 
						      200,250,300,350,400,450,500,1000};

  const int nBinPt2 = 50; Float_t xbinsPt2[nBinPt2+1]; for(int i=0; i<=nBinPt2;i++) xbinsPt2[i] = i * 1.0;

  const int nBinRap = 24; Float_t xbinsRap[nBinRap+1]; for(int i=0; i<=nBinRap;i++) xbinsRap[i] = i * 0.1;

  const int nBinPtRap0 = 36; Float_t xbinsPtRap0[nBinPtRap0+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinPtRap1 = 36; Float_t xbinsPtRap1[nBinPtRap1+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinPtRap2 = 36; Float_t xbinsPtRap2[nBinPtRap2+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinPtRap3 = 36; Float_t xbinsPtRap3[nBinPtRap3+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                		   10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                		   20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                		  100,150,200,300,400,500,1000};
  const int nBinPtRap4 = 36; Float_t xbinsPtRap4[nBinPtRap4+1] = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                                  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                                                                  20, 25, 30, 35, 40, 50, 60, 70, 80, 90,
                                                                 100,150,200,300,400,500,1000};

  gStyle->SetOptStat(0);
  int version = 1; int alternative = 0;
  if     (whichDY == 0) { version = 0; alternative = 1;}
  else if(whichDY == 2) { version = 2; alternative = 1;}

  TH1D* histoSyst[7];
  if       (theHistName == "Pt"){
    histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPt, xbinsPt);
    histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPt, xbinsPt);
    histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPt, xbinsPt);
    histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPt, xbinsPt);
    histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPt, xbinsPt);
    histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPt, xbinsPt);
    histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPt, xbinsPt);
  } else if(theHistName == "Pt2"){
    histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPt2, xbinsPt2);
    histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPt2, xbinsPt2);
    histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPt2, xbinsPt2);
    histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPt2, xbinsPt2);
    histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPt2, xbinsPt2);
    histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPt2, xbinsPt2);
    histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPt2, xbinsPt2);
  } else if(theHistName == "Rap"){
    histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinRap, xbinsRap);
    histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinRap, xbinsRap);
    histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinRap, xbinsRap);
    histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinRap, xbinsRap);
    histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinRap, xbinsRap);
    histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinRap, xbinsRap);
    histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinRap, xbinsRap);
  } else if(theHistName == "PtRap0"){
   histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPtRap0, xbinsPtRap0);
   histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPtRap0, xbinsPtRap0);
   histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPtRap0, xbinsPtRap0);
   histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPtRap0, xbinsPtRap0);
   histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPtRap0, xbinsPtRap0);
   histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPtRap0, xbinsPtRap0);
   histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPtRap0, xbinsPtRap0);
  } else if(theHistName == "PtRap1"){
   histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPtRap1, xbinsPtRap1);
   histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPtRap1, xbinsPtRap1);
   histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPtRap1, xbinsPtRap1);
   histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPtRap1, xbinsPtRap1);
   histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPtRap1, xbinsPtRap1);
   histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPtRap1, xbinsPtRap1);
   histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPtRap1, xbinsPtRap1);
  } else if(theHistName == "PtRap2"){
   histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPtRap2, xbinsPtRap2);
   histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPtRap2, xbinsPtRap2);
   histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPtRap2, xbinsPtRap2);
   histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPtRap2, xbinsPtRap2);
   histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPtRap2, xbinsPtRap2);
   histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPtRap2, xbinsPtRap2);
   histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPtRap2, xbinsPtRap2);
  } else if(theHistName == "PtRap3"){
   histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPtRap3, xbinsPtRap3);
   histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPtRap3, xbinsPtRap3);
   histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPtRap3, xbinsPtRap3);
   histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPtRap3, xbinsPtRap3);
   histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPtRap3, xbinsPtRap3);
   histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPtRap3, xbinsPtRap3);
   histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPtRap3, xbinsPtRap3);
  } else if(theHistName == "PtRap4"){
   histoSyst[0] = new TH1D(Form("histoSyst_%d",0),  Form("histoSyst_%d",0),  nBinPtRap4, xbinsPtRap4);
   histoSyst[1] = new TH1D(Form("histoSyst_%d",1),  Form("histoSyst_%d",1),  nBinPtRap4, xbinsPtRap4);
   histoSyst[2] = new TH1D(Form("histoSyst_%d",2),  Form("histoSyst_%d",2),  nBinPtRap4, xbinsPtRap4);
   histoSyst[3] = new TH1D(Form("histoSyst_%d",3),  Form("histoSyst_%d",3),  nBinPtRap4, xbinsPtRap4);
   histoSyst[4] = new TH1D(Form("histoSyst_%d",4),  Form("histoSyst_%d",4),  nBinPtRap4, xbinsPtRap4);
   histoSyst[5] = new TH1D(Form("histoSyst_%d",5),  Form("histoSyst_%d",5),  nBinPtRap4, xbinsPtRap4);
   histoSyst[6] = new TH1D(Form("histoSyst_%d",6),  Form("histoSyst_%d",6),  nBinPtRap4, xbinsPtRap4);
  } else {
    printf("Wrong option!\n");
    return;
  }

  TFile *_file0 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_default.root",theHistName.Data(),nsel,version)); // default NLO/LO
  TFile *_file1 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_default.root",theHistName.Data(),nsel,alternative)); // alternative LO/NLO
  TFile *_file2 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_svd30.root",theHistName.Data(),nsel,version)); // SVD alt
  TFile *_file3 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_bayes10.root",theHistName.Data(),nsel,version)); // Bayes alt
  TFile *_file4 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_momres.root",theHistName.Data(),nsel,version)); // MonRes
  TFile *_file5 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_pdf.root",theHistName.Data(),nsel,version)); // PDF bkg
  TFile *_file6 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_qcd.root",theHistName.Data(),nsel,version)); // QCD bkg.
  TFile *_file7 = TFile::Open(Form("output/histoUnfolding%s_nsel%d_dy%d_rebin1_lepeff.root",theHistName.Data(),nsel,version)); // LepEff

  TH1D* histDef = (TH1D*)_file0->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  TH1D* histAlt[7];
  histAlt[0] = (TH1D*)_file1->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  histAlt[1] = (TH1D*)_file2->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  histAlt[2] = (TH1D*)_file3->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  histAlt[3] = (TH1D*)_file4->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  histAlt[4] = (TH1D*)_file5->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  histAlt[5] = (TH1D*)_file6->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));
  histAlt[6] = (TH1D*)_file7->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));

  double systVal[9][nBinPt],systTotalVal[nBinPt];

  for(int i=1; i<=histDef->GetNbinsX(); i++){
    for(int j=0; j<7; j++){
      systVal[j][i] = 100.*TMath::Abs(histDef->GetBinContent(i)-histAlt[j]->GetBinContent(i))/histDef->GetBinContent(i);
    }
    if      (histDef->GetBinCenter(i)+histDef->GetBinWidth(i)/2<=40  && systVal[0][i] > 1.0) systVal[0][i] = 1.0;
    else if (histDef->GetBinCenter(i)+histDef->GetBinWidth(i)/2<=100 && systVal[0][i] < 0.1) systVal[0][i] = 0.5;
    systVal[7][i] = 100.*histDef->GetBinError(i)/histDef->GetBinContent(i); // stat
    systVal[8][i] = 100.*0.023;
    printf("(%2d) %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",i,systVal[0][i],systVal[1][i],systVal[2][i],systVal[3][i],systVal[4][i],systVal[5][i],systVal[6][i],systVal[7][i]);
    systVal[1][i] = 0.0;
    systVal[2][i] = 0.0;
    histoSyst[0]->SetBinContent(i,sqrt(systVal[0][i]*systVal[0][i]+systVal[1][i]*systVal[1][i]+systVal[2][i]*systVal[2][i]+
                                       systVal[3][i]*systVal[3][i]+systVal[4][i]*systVal[4][i]+systVal[5][i]*systVal[5][i]+
                                       systVal[6][i]*systVal[6][i]+systVal[7][i]*systVal[7][i]+systVal[8][i]*systVal[8][i]
                                       ));
    histoSyst[1]->SetBinContent(i,sqrt(systVal[0][i]*systVal[0][i]+systVal[1][i]*systVal[1][i]+systVal[2][i]*systVal[2][i])); // unfolding
    histoSyst[2]->SetBinContent(i,systVal[3][i]); // momres
    histoSyst[3]->SetBinContent(i,sqrt(systVal[4][i]*systVal[4][i]+systVal[5][i]*systVal[5][i])); // Bkg.
    histoSyst[4]->SetBinContent(i,systVal[6][i]); // LepEff
    histoSyst[5]->SetBinContent(i,systVal[7][i]); // stat
    histoSyst[6]->SetBinContent(i,systVal[8][i]); // lumi
  }

  atributes(histoSyst[0],"Z p_{T} [GeV]",1,"Uncertainty (%)", 1);
  atributes(histoSyst[1],"Z p_{T} [GeV]",2,"Uncertainty (%)", 1);
  atributes(histoSyst[2],"Z p_{T} [GeV]",4,"Uncertainty (%)", 1);
  atributes(histoSyst[3],"Z p_{T} [GeV]",5,"Uncertainty (%)", 1);
  atributes(histoSyst[4],"Z p_{T} [GeV]",6,"Uncertainty (%)", 1);
  atributes(histoSyst[5],"Z p_{T} [GeV]",7,"Uncertainty (%)", 1);
  atributes(histoSyst[6],"Z p_{T} [GeV]",8,"Uncertainty (%)", 1);

  histoSyst[0]->SetMinimum(0.0);
  histoSyst[0]->Draw();
  histoSyst[1]->Draw("same,hist");
  histoSyst[2]->Draw("same,hist");
  histoSyst[3]->Draw("same,hist");
  histoSyst[4]->Draw("same,hist");
  histoSyst[5]->Draw("same,hist");
  histoSyst[6]->Draw("same,hist");

 TLatex * CMSLabel = new TLatex (0.10, 0.91, "#bf{CMS} Preliminary");
 CMSLabel->SetNDC ();
 CMSLabel->SetTextAlign (10);
 CMSLabel->SetTextFont (42);
 CMSLabel->SetTextSize (0.04);
 CMSLabel->Draw ("same") ;
 _lumiLabel = new TLatex (0.90, 0.91, "36.8 fb^{-1} (13 TeV)");
 _lumiLabel->SetNDC ();
 _lumiLabel->SetTextAlign (30);
 _lumiLabel->SetTextFont (42);
 _lumiLabel->SetTextSize (0.04);
 _lumiLabel->Draw ("same") ;
 TLegend* leg = new TLegend(0.20,0.60,0.50,0.80);                                                    
 leg ->SetFillStyle(0);
 leg ->SetFillColor(kWhite);
 leg ->SetBorderSize(0);
 leg->SetTextSize(0.035);                                                                         
 leg->AddEntry(histoSyst[0],"Total Unc.","l");
 leg->AddEntry(histoSyst[1],"Unfolding","l");
 leg->AddEntry(histoSyst[2],"Mom. Res.","l");
 leg->AddEntry(histoSyst[3],"Bkg.","l");
 leg->AddEntry(histoSyst[4],"Identification","l");
 leg->AddEntry(histoSyst[5],"Statistical","l");
 leg->AddEntry(histoSyst[6],"Luminosity","l");
 leg->Draw();

}
