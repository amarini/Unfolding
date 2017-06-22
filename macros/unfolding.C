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
void helper_function(int nsel=0,int LO=1, int rebin=1){
  TFile *_file0;
  TFile *_file1;
  TFile *_file2;
  if (LO == 0){
    _file0 = TFile::Open("inputs_finerbinning/histozllPtRec_NLO.root");
    _file1 = TFile::Open("inputs_finerbinning/histozllPtRecGen_NLO.root");
    _file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_NLO.root");
  }
  else{
    _file0 = TFile::Open("inputs_finerbinning/histozllPtRec_LO.root");
    _file1 = TFile::Open("inputs_finerbinning/histozllPtRecGen_LO.root");
    _file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_LO.root");
  }

TH1D* bdat = (TH1D*)_file0->Get(Form("histoPtRecDA_%d",nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histoPtRecDY_%d",nsel));

TH2D* Adet = (TH2D*)_file1->Get(Form("histoPtRecGen_%d",nsel));

TString inputFileHist = "hDDilLowPtMM";
if(nsel == 1) inputFileHist = "hDDilLowPtEE";
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
RooUnfoldBayes theRooUnfoldBayes(&theResponse, bdat, 4);
//RooUnfoldSvd theRooUnfoldBayes(&theResponse, bdat, 20);
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
sprintf(output,"histoUnfolding_nsel%d_lo%d.root",nsel,LO); 
TFile* outFilePlots = new TFile(output,"recreate");
outFilePlots->cd();
hReco->Write();
xini->Write();
outFilePlots->Close();

}

void unfolding(int rebin=10){
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      helper_function(i,j,rebin);

    }
  }
}
