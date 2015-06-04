#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "common.h"

void electrons()
{
  TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
  int true_zee = 0, true_zmm = 0;
  
  bind_attributes(tree, p1, p2, e1, e2, m1, m2);

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("mee", &me); 
  tree->SetBranchAddress("is_true_Zee", &true_zee); 
  tree->SetBranchAddress("is_true_Zmm", &true_zmm);

  TH1F* h = new TH1F("h", "m_{\\gamma \\gamma} \\mbox{ distribution};m_{\\gamma \\gamma} \\mbox{ (GeV)};\\mbox{Count}", 500, 0, 500); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "E", 500, 0, 25);
  //TH1F* h = new TH1F("h", "phi", 500, -5, 5);
  //TH1F* h2 = new TH1F("h2", "eta", 500, -5, 5);  
//TH1F* h2 = new TH1F("h2", "eta_true", 500, -5, 5);
  
  //h2->SetLineColor(kRed);

  int totalEntries = 0;
  int looseEntries = 0;
  int keptEntries = 0;

  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    //if(e1.true_E < 0) continue;
    tree->GetEntry(i);
    if(e1.p > 30 && e2.p > 30) h->Fill(me);

  }

  h->Sumw2();
  //h2->Sumw2();

  //TH1F *r = (TH1F *)h->Clone();
  //r->Divide(h, h2, 1., 1., "B");
  h->Draw();

  printf("distrib mgg: %.2f %.2f\n", avg_m, sigma_m);
  printf("loose ratio: %.6f\n", float(keptEntries)/float(totalEntries));
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
