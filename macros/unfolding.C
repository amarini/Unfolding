#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include "TUnfoldDensity.h"

// nsel = 0 (mm) 1 (ee)
// whichDY = 0 (MG), 1 (aMCNLO), 2 (POWHEG)
void helper_function(int nsel=0,int whichDY=0, int rebin=1, bool isPt2 = false, TString suffix = ""){
  TString theInputFolder = "inputs";
  TString theHistName = "";
  if(isPt2 == true){
    theInputFolder = "inputs";
    theHistName = "2";
  }
  TFile *_file0;
  TFile *_file1;
  TFile *_file2;
  double theXS = 2008.4*3;
  double theLumi = 35800.0;
  if     (whichDY == 0){
    _file0 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_LO.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),1,theHistName.Data()));
  }
  else if(whichDY == 1){
    _file0 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_NLO.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),2,theHistName.Data()));
  }
  else if(whichDY == 2 && nsel == 0){
    _file0 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToMM_POWHEG.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),1,theHistName.Data()));
    theXS = 1975.0;
  }
  else if(whichDY == 2 && nsel == 1){
    _file0 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),whichDY,theHistName.Data()));
    _file1 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToEE_POWHEG.root");
    _file2 = TFile::Open(Form("%s/histoDY%dzllPt%sRecGen.root",theInputFolder.Data(),1,theHistName.Data()));
    theXS = 1975.0;
  }

TH1D* bdat = (TH1D*)_file0->Get(Form("histoPt%sRecDA_%d",theHistName.Data(),nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histoPt%sRecDY_%d",theHistName.Data(),nsel));

TH2D* Adet = (TH2D*)_file0->Get(Form("histoPt%sRecGen_%d",theHistName.Data(),nsel));

TString inputFileHist = Form("hDDilPt%sMM",theHistName.Data());
if(nsel == 1) inputFileHist = Form("hDDilPt%sEE",theHistName.Data());
TH1D* xini  = (TH1D*)_file1->Get(inputFileHist.Data());

TH1D* hDNEvt  = (TH1D*)_file1->Get("hDTotalMCWeight");

xini->Scale(theXS*theLumi/hDNEvt->GetSumOfWeights());

bdat->Rebin(rebin);
bini->Rebin(rebin);
Adet->Rebin2D(rebin,rebin);
xini->Rebin(rebin);

printf("data: %f\n",bdat->GetSumOfWeights());
printf("  DY: %f\n",bini->GetSumOfWeights());
printf("  2D: %f\n",Adet->GetSumOfWeights());
printf(" gen: %f\n",xini->GetSumOfWeights());

RooUnfoldResponse theResponse(bini, xini, Adet);
theResponse.UseOverflow(kTRUE);
RooUnfoldBayes theRooUnfoldBayes(&theResponse, bdat, 5);
//RooUnfoldSvd theRooUnfoldBayes(&theResponse, bdat, 30);
//RooUnfoldTUnfold theRooUnfoldBayes(&theResponse, bdat);
theRooUnfoldBayes.SetVerbose(-1);
theRooUnfoldBayes.SetNToys(1000);

TH1D* hReco= (TH1D*) theRooUnfoldBayes.Hreco(RooUnfold::kCovToy);
theRooUnfoldBayes.PrintTable (cout, xini);
hReco->Draw();
hReco->Scale(1,"width");
//bdat->Draw("SAME");
xini->SetLineColor(2);
xini->Scale(1,"width");
xini->Draw("SAME");

char output[100];
if(isPt2 == false){
  sprintf(output,"histoUnfolding_nsel%d_dy%d_rebin%d%s.root",nsel,whichDY,rebin,suffix.Data());
} else {
  sprintf(output,"histoUnfoldingPt2_nsel%d_dy%d_rebin%d%s.root",nsel,whichDY,rebin,suffix.Data());
}
TFile* outFilePlots = new TFile(output,"recreate");
outFilePlots->cd();
hReco->Write();
xini->Write();
outFilePlots->Close();

}

void unfolding(TString suffix = "", int rebin=1){
  for (int i=0;i<=1;i++){
    for (int j=0;j<=2;j++){
      helper_function(i,j,rebin,true,suffix.Data());

    }
  }
}
