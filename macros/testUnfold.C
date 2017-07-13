#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TUnfold.h"


// nsel = 0 (mm) 1 (ee)
// whichDY = 0 (MG), 1 (aMCNLO), 2 (POWHEG)
void helper_function(int nsel=0,int whichDY=0, int rebin=1, TString theHistName = "Pt", TString suffix = ""){
  TString theInputFolder = "inputs";

  TFile *_file0;
  TFile *_file1;
  TFile *_file2;
  double theXS = 2008.4*3;
  double theLumi = 35800.0;
  if     (whichDY == 0){
    _file0 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_LO.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),1,theHistName.Data()));
  }
  else if(whichDY == 1){
    _file0 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_NLO.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),2,theHistName.Data()));
  }
  else if(whichDY == 2 && nsel == 0){
    _file0 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToMM_POWHEG.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),1,theHistName.Data()));
    theXS = 1975.0;
  }
  else if(whichDY == 2 && nsel == 1){
    _file0 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToEE_POWHEG.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzll%sRecGen.root",theInputFolder.Data(),1,theHistName.Data()));
    theXS = 1975.0;
  }

TH1D* bdat = (TH1D*)_file0->Get(Form("histo%sRecDA_%d",theHistName.Data(),nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histo%sRecDY_%d",theHistName.Data(),nsel));

TH2D* Adet = (TH2D*)_file0->Get(Form("histo%sRecGen_%d",theHistName.Data(),nsel));

TString inputFileHist = Form("hDDil%sMM",theHistName.Data());
if(nsel == 1) inputFileHist = Form("hDDil%sEE",theHistName.Data());
TH1D* xini  = (TH1D*)_file1->Get(inputFileHist.Data());

TH1D* hDNEvt  = (TH1D*)_file1->Get("hDTotalMCWeight");

xini->Scale(theXS*theLumi/hDNEvt->GetSumOfWeights());

// add non-reconstructed fiducial events in gen bin 0
for(int i=1; i<=Adet->GetNbinsY(); i++) { // GEN
  double genBin = xini->GetBinContent(i);
  double recBin = 0;
  for(int j=1; j<=Adet->GetNbinsX(); j++) { // RECO
    recBin = recBin + Adet->GetBinContent(j,i);
  }
  if(genBin-recBin <= 0) printf("Bin %d totally efficient\n",i);
  Adet->SetBinContent(0,i,TMath::Max(genBin-recBin,0.0));
}
for(int i=1; i<=Adet->GetNbinsY(); i++) {
 Adet->SetBinContent(i,Adet->GetNbinsY()+1,0);
 Adet->SetBinContent(Adet->GetNbinsX()+1,i,0);
}

// subtract reconstructed non-fiducial events in data
for(int i=1; i<=Adet->GetNbinsX(); i++) { // RECO
  double allRecBin = bini->GetBinContent(i);
  double fidRecBin = 0;
  for(int j=1; j<=Adet->GetNbinsY(); j++) { // GEN
    fidRecBin = fidRecBin + Adet->GetBinContent(i,j);
  }
  double gendiff = allRecBin - fidRecBin;
  if(gendiff <= 0 || bdat->GetBinContent(i)-gendiff <= 0) printf("Bin %d totally efficient\n",i);
  //printf("dat(%2d) = %f -> ",i,bdat->GetBinContent(i));
  bdat->SetBinContent(i,TMath::Max(bdat->GetBinContent(i)-gendiff,0.0));
  //printf("%f\n",bdat->GetBinContent(i));
}

//bdat->Sumw2();
//bini->Sumw2();
//Adet->Sumw2();
//xini->Sumw2();
bdat->Rebin(rebin);
bini->Rebin(rebin);
Adet->Rebin2D(rebin,rebin);
xini->Rebin(rebin);

printf("data: %f\n",bdat->GetSumOfWeights());
printf("  DY: %f\n",bini->GetSumOfWeights());
printf("  2D: %f\n",Adet->GetSumOfWeights());
printf(" gen: %f\n",xini->GetSumOfWeights());

//Do Unfolding
TUnfold unfold(Adet,TUnfold::kHistMapOutputHoriz);
Double_t tau=1.E-4;
Double_t biasScale=0.0;
unfold.DoUnfold(tau,bdat,biasScale);

TH1D *unfoldedDistribution = (TH1D*)bdat->Clone();
unfoldedDistribution->Reset();
unfold.GetOutput(unfoldedDistribution);

// TH1D *x = (TH1D*)(unfold.GetOutput("x","myVariable",0,0));
// TH2D *rhoij=unfold.GetRhoIJ("correlation","myVariable");

unfoldedDistribution->Scale(1,"width");
bini->Scale(1,"width");

unfoldedDistribution->Draw();
bini->Draw("same");
bini->SetLineColor(kRed);
}

void testUnfold(TString theHistName = "Pt", TString suffix = "", int rebin=1){
  for (int i=0;i<=0;i++){
    for (int j=1;j<=1;j++){
      helper_function(i,j,rebin,theHistName.Data(),suffix.Data());

    }
  }
}
