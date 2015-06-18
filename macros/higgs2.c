#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"

  float E_min = 175, E_max = 700; // 150 ; 475; 14

bool keep_photon(photon &p)
{
    return p.tight && p.E > E_min && p.E < E_max /*&& p.convFlag <= 0*/ && abs(p.E-p.true_E) < 50 && abs(float(p.true_eta)) > 1.5f;
}

void higgs2()
{
//  TFile* f = TFile::Open("small.root");  // open the file
//    TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
  //  TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_ZH900_Pythia_NoMassCut.root");  // open the file
 // TTree* tree = (TTree*)f->Get("tree");

  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
 
  const int N = 6;
  TFile *files[N];
  TTree *trees[N];
  TH1D *histograms[N];
  int count[N];
  
 /*      p0           1.05852e+01     fixed    
   2  p1          -2.37115e-02     fixed    
   3  p2           1.80387e+02   2.98922e+01   9.44240e-02   4.40528e-05
   4  p3           1.26228e+02   3.40480e-01   1.07401e-03   2.40111e-03
   5  p4           1.80000e+00     fixed */
  
  double cross_sections[N] = { 22.3*2.3e-3, 8.19759 * 5.51e-5, 4.082711*1.25e-5, 3.191034*2.65e-6, 1.4391*3.11e-7, 0.61988*5e-8 };
  double masses[N] = { 125, 200, 300, 400, 500, 600 };
  
  for(int k = 0; k < N; ++k) cross_sections[k] *= 500 / (22.3*2.3e-3);
  
  
  TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
  c1->SetLogy();
 // c1->SetLogx();
  c1->SetGrid();
  
  TGraph *gr = new TGraph(N,masses,cross_sections);
  gr->SetTitle("N(m_{H});m_{H};\\log N");
  gr->SetMarkerStyle(3);
  gr->SetMarkerSize(2);
  gr->Draw("ACP");
  
  
   /*h->SetLineColor(kBlack);
   h->SetFillColorAlpha(kWhite, 0.2);
   h->Draw("same");*/
}
