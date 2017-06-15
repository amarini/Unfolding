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
//LO=0 (LO) 1(NLO)
void unfolding(int nsel = 0,int LO=0){
  if (LO==0){
    TFile *_file0 = TFile::Open("/afs/cern.ch/user/c/christiw/public/panda/MitAnalysisRunII/panda/macros/80x/LO/histozllPtRec.root");
    TFile *_file1 = TFile::Open("/afs/cern.ch/user/c/christiw/public/panda/MitAnalysisRunII/panda/macros/80x/NLO/histozllPtRecGen_NLO.root");

  }
  else{
    TFile *_file0 = TFile::Open("/afs/cern.ch/user/c/christiw/public/panda/MitAnalysisRunII/panda/macros/80x/NLO/histozllPtRec_NLO.root");
    TFile *_file1 = TFile::Open("/afs/cern.ch/user/c/christiw/public/panda/MitAnalysisRunII/panda/macros/80x/LO/histozllPtRecGen_LO.root");

  }

TH1D* bdat = (TH1D*)_file0->Get(Form("histoPtRecDA_%d",nsel));
TH1D* bini = (TH1D*)_file0->Get(Form("histoPtRecDY_%d",nsel));


TH2D* Adet = (TH2D*)_file1->Get(Form("histoPtRecGen_%d",nsel));

TFile *_file2 = TFile::Open("/afs/cern.ch/work/c/ceballos/public/samples/panda/v_004_0/DYJetsToLL_M-50_NLO.root");

TString inputFileHist = "hDDilPtMM";
if(nsel == 1) inputFileHist = "hDDilPtEE";
TH1D* xini  = (TH1D*)_file2->Get(inputFileHist.Data());

TH1D* hDNEvt  = (TH1D*)_file2->Get("hDTotalMCWeight");
xini->Scale(2008.4*3*35800./hDNEvt->GetSumOfWeights());

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
hReco->Scale(1,"width");
//bdat->Draw("SAME");
xini->SetLineColor(8);
xini->Scale(1,"width");
xini->Draw("SAME");

//TSVDUnfold *tsvdunf = new TSVDUnfold( bdat, bini, xini, Adet );
//TH1D* unfresult = tsvdunf->Unfold( 2 );
//unfresult->Draw();
//xini->Draw("e,same");
}
