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
// LO = 0 (NLO) 1 (LO)
void helper_function(int nsel=0,int LO=1, int rebin=1, bool isFine = false, TString suffix = ""){
  TString theInputFolder = "inputs";
  TString theHistName = "";
  if(isFine == true){
    theInputFolder = "inputs_finerbinning";
    theHistName = "Low";
  }
  TFile *_file0;
  TFile *_file1;
  TFile *_file2;
  if (LO == 0){
    _file0 = TFile::Open(Form("%s/histozllPtRec_NLO.root",theInputFolder.Data()));
    _file1 = TFile::Open(Form("%s/histozllPtRecGen_NLO.root",theInputFolder.Data()));
    _file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_NLO.root");
  }
  else{
    _file0 = TFile::Open(Form("%s/histozllPtRec_LO.root",theInputFolder.Data()));
    _file1 = TFile::Open(Form("%s/histozllPtRecGen_LO.root",theInputFolder.Data()));
    _file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_LO.root");
  }

TH1D* bdat = (TH1D*)_file0->Get(Form("histoPtRecDA_%d",nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histoPtRecDY_%d",nsel));

TH2D* Adet = (TH2D*)_file1->Get(Form("histoPtRecGen_%d",nsel));

TString inputFileHist = Form("hDDil%sPtMM",theHistName.Data());
if(nsel == 1) inputFileHist = Form("hDDil%sPtEE",theHistName.Data());
TH1D* xini  = (TH1D*)_file2->Get(inputFileHist.Data());

TH1D* hDNEvt  = (TH1D*)_file2->Get("hDTotalMCWeight");

xini->Scale(2008.4*3*35800./hDNEvt->GetSumOfWeights());

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
//RooUnfoldBayes theRooUnfoldBayes(&theResponse, bdat, 5);
RooUnfoldSvd theRooUnfoldBayes(&theResponse, bdat, 30);
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
sprintf(output,"histoUnfolding_nsel%d_lo%d_rebin%d%s.root",nsel,LO,rebin,suffix.Data()); 
TFile* outFilePlots = new TFile(output,"recreate");
outFilePlots->cd();
hReco->Write();
xini->Write();
outFilePlots->Close();

}

void unfolding(int rebin=1, TString suffix = ""){
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      helper_function(i,j,rebin,false,suffix.Data());

    }
  }
}
