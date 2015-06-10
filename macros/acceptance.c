#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"

void acceptance()
{
  TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_NoMassCut.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron el1, el2;
  muon m1, m2;
  int true_zee = 0, true_zmm = 0;
  
  bind_attributes(tree, p1, p2, el1, el2, m1, m2);

  float avg_m, sigma_m;

  //TH1F* h = new TH1F("h", "E", 50, 25, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //TH1F* h2 = new TH1F("h2", "E", 50, 25, 250);
  TH1F* h = new TH1F("h", "eta", 50, -3.5, 3.5);
  TH1F* h2 = new TH1F("h2", "eta", 75, -3.5, 3.5);  
//TH1F* h2 = new TH1F("h2", "eta_true", 50, -3, 3);
  
  //h2->SetLineColor(kRed);

  int totalEntries = 0;
  int looseEntries = 0;
  int keptEntries = 0;
 
 unsigned int count = tree->GetEntries();

 for (unsigned int i = 0; i < count; i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);

    //if(p1.mother==25 && p2.mother==25)
    {
        if(p1.tight && p2.tight)
        {
            vec imp1, imp2;
            imp1.pr2c(p1.E, p1.phi, p1.eta);
            imp2.pr2c(p2.E, p2.phi, p2.eta);
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
            
            if(p1.E > 0 && p2.E > 0)
            {
                float mass = invMass(p1, p2);
                if(mass > 122 && mass < 128) h->Fill(imp1.z);
                
            }
            //h->Fill(invMass(p1, p2));

        }
    }

  }

  h->Sumw2();
  /*h2->Sumw2();

  TH1F *r = (TH1F *)h->Clone();
  r->Divide(h, h2, 1., 1., "B");
  r->Draw();

  printf("distrib mgg: %.2f %.2f\n", avg_m, sigma_m);
  printf("loose ratio: %.6f\n", float(looseEntries)/float(totalEntries));*/
 h->Draw(); // plot the histogram
  //h2->Draw();
}
