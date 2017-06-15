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

void unfolding(int nsel = 0){

TFile *_file0 = TFile::Open("inputs/histozllPtRec.root");

TH1D* bdat = (TH1D*)_file0->Get(Form("histoPtRecDA_%d",nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histoPtRecDY_%d",nsel));

TFile *_file1 = TFile::Open("inputs/histozllPtRecGen.root");

TH2D* Adet = (TH2D*)_file1->Get(Form("histoPtRecGen_%d",nsel));

TFile *_file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_LO.root");

TString inputFileHist = "hDDilPtMM";
if(nsel == 1) inputFileHist = "hDDilPtEE";
TH1D* xini  = (TH1D*)_file2->Get(inputFileHist.Data());
xini->Scale(2008.4*3*35800./49144272.);

printf("%f\n",bdat->GetSumOfWeights());
printf("%f\n",bini->GetSumOfWeights());
printf("%f\n",Adet->GetSumOfWeights());
printf("%f\n",xini->GetSumOfWeights());

int nReg = 5;
RooUnfoldResponse theResponse(bini, xini, Adet);
RooUnfoldBayes theRooUnfoldBayes(&theResponse, bdat, nReg);

TH1D* hReco= (TH1D*) theRooUnfoldBayes.Hreco();
theRooUnfoldBayes.PrintTable (cout, xini);
hReco->Draw();
//bdat->Draw("SAME");
xini->SetLineColor(8);
xini->Draw("SAME");

//TSVDUnfold *tsvdunf = new TSVDUnfold( bdat, bini, xini, Adet );
//TH1D* unfresult = tsvdunf->Unfold( 2 );
//unfresult->Draw();
//xini->Draw("e,same");
}
