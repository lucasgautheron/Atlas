#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"

void tests()
{
  TFile* f = TFile::Open("small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
  
  bind_attributes(tree, p1, p2, e1, e2, m1, m2);

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("mee", &me); 

  TH1F* h = new TH1F("h", "eta", 500, -5, 5); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "invariant mass", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //h2->SetLineColor(kRed);

  int totalEntries = 0;
  int looseEntries = 0;
  int keptEntries = 0;
  
  /*int n = 0, k = 0;
  
  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);
    if(p1.true_mother==25 && p2.true_mother==25)
    {
        if(p1.true_E > 0)
        {
            ++n;
        }
    }
  }
  
  float *x = new float[n+1], *y = new float[n+1];*/

  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);
    if(p1.true_mother==25 && p2.true_mother==25)
    {
        if(p1.true_E > 0)
        {
            float px = p1.true_E * TMath::Cos(p1.true_phi) + p2.true_E * TMath::Cos(p2.true_phi);
            float py = p1.true_E * TMath::Sin(p1.true_phi) + p2.true_E * TMath::Sin(p2.true_phi);
            float pz = p1.true_E * TMath::SinH(p1.true_eta) + p2.true_E * TMath::SinH(p2.true_eta) ;
            float p = TMath::Sqrt(px*px + py*py + pz*pz);
            float eta = TMath::ATanH(pz/p);
            //y[k] = invMassTrue(p1, p2);
            //x[k] = eta;
            //++k;
            h->Fill(eta);
            //h->Fill(p2.true_eta);
        }
    }
  }
  
  /*TGraph *gr = new TGraph(n,x,y);
  gr->Draw("ACP");*/

  printf("tight ratio: %.6f\n", float(keptEntries)/float(totalEntries));
  
    //h->Draw();
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
