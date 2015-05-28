#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"

void eta()
{
  TFile* f = TFile::Open("small.root");  // open the file
//    TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
//    TFile* f = TFile::Open("small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
  
  bind_attributes(tree, p1, p2, e1, e2, m1, m2);

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("mee", &me); 
  
  const int n = 75;
  float E_min = 0, E_max = 3.142;
  
  TF1 *fa = new TF1("fa1","(cos(x/2)/sin(x/2)  ) / ((cos(x/2)/sin(x/2)) *(cos(x/2)/sin(x/2)) + 1)",0,3.142);

  TH1F* h = new TH1F("h", "eta", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "invariant mass gamma gamma", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
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
  unsigned int count = tree->GetEntries();

  for (unsigned int i = 0; i < count; i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);

    //if(p1.mother==25 && p2.mother==25)
    {
        //if(p1.true_E > 0)
        {
            vec imp1, imp2;
            imp1.pr2c(p1.true_E, p1.true_phi, p1.true_eta);
            imp2.pr2c(p2.true_E, p2.true_phi, p2.true_eta);
            imp1.add(imp2);
            float p = imp1.norm();
            imp1.c2pr(imp1.x, imp1.y, imp1.z);
            /*float px = p1.true_E * TMath::Cos(p1.true_phi) + p2.true_E * TMath::Cos(p2.true_phi);
            float py = p1.true_E * TMath::Sin(p1.true_phi) + p2.true_E * TMath::Sin(p2.true_phi);
            float pz = p1.true_E * TMath::SinH(p1.true_eta) + p2.true_E * TMath::SinH(p2.true_eta) ;
            float p = TMath::Sqrt(px*px + py*py + pz*pz);
            float eta = TMath::ATanH(pz/p);*/
            //y[k] = invMassTrue(p1, p2);
            //x[k] = eta;
            //++k;
            //printf("%.2f\n", imp1.z - eta);
            //float theta = 2*atan(exp(-imp1.z));
            // todo: par unité d'angle solide en fct de theta ?
            h->Fill(2*atan(exp(-imp1.z))/*+TMath::Pi()/2.0f*/);

        }
    }
  }
  

  printf("tight ratio: %.6f\n", float(keptEntries)/float(totalEntries));
    h->Multiply(fa);
    h->Draw("E");
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
