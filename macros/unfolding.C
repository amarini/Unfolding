#include <TFile.h>
#include <TSystem.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TSVDUnfold.h>
#include <TROOT.h>
// nsel = 0 (mm) 1 (ee)
//LO=0 (NLO) 1(LO)
TFile* outFilePlots = new TFile("histoUnfolding.root","recreate");

void helper_function(int nsel=0,int LO=1){
  TFile *_file0;
  TFile *_file1;
  TFile *_file2;
  if (LO == 0){
    _file0 = TFile::Open("inputs/histozllPtRec_NLO.root");
    _file1 = TFile::Open("inputs/histozllPtRecGen_NLO.root");
  _file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_NLO.root");
  }
  else{
    _file0 = TFile::Open("inputs/histozllPtRec_LO.root");
    _file1 = TFile::Open("inputs/histozllPtRecGen_LO.root");
    _file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_LO.root");
  }

TH1D* bdat = (TH1D*)_file0->Get(Form("histoPtRecDA_%d",nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histoPtRecDY_%d",nsel));

TH2D* Adet = (TH2D*)_file1->Get(Form("histoPtRecGen_%d",nsel));

TString inputFileHist = "hDDilPtMM";
if(nsel == 1) inputFileHist = "hDDilPtEE";
TH1D* xini  = (TH1D*)_file2->Get(inputFileHist.Data());

TH1D* hDNEvt  = (TH1D*)_file2->Get("hDTotalMCWeight");
xini->Scale(2008.4*3*35800./hDNEvt->GetSumOfWeights());

printf("data: %f\n",bdat->GetSumOfWeights());
printf("  DY: %f\n",bini->GetSumOfWeights());
printf("  2D: %f\n",Adet->GetSumOfWeights());
printf(" gen: %f\n",xini->GetSumOfWeights());

int nReg = 5;
RooUnfoldResponse theResponse(bini, xini, Adet);
RooUnfoldBayes theRooUnfoldBayes(&theResponse, bdat, nReg);
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
//CLONE
char reco_name[100];
TString input_name;

if (LO==0){
   input_name=inputFileHist+"_NLO";
   sprintf(reco_name,"reco_%d_NLO",nsel);
 }
else{
   input_name=inputFileHist+"_LO";
   sprintf(reco_name,"reco_%d_LO",nsel);

 }
TH1F *xini2 = (TH1F*)xini->Clone(input_name);
TH1F *hReco2 = (TH1F*)hReco->Clone(reco_name);

outFilePlots->cd();





hReco2->Write();
xini2->Write();
}


void unfolding(){
  
  for (int i=0;i<2;i++){
    for (int j=0;j<2;j++){
      helper_function(i,j);

    }
  }
  outFilePlots->Close();

}
