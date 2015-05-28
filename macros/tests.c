#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"

void tests()
{
  TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_NoMassCut.root");  // open the file
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
  
  const int n = 100;
  float E_min = 100, E_max = 150;
  float x[n], y[n], p[n];

  TH1F* h = new TH1F("h", "invariant mass gamma gamma", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "invariant mass gamma gamma", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
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
            /*vec imp1, imp2;
            imp1.pr2c(p1.E, p1.phi, p1.eta);
            imp2.pr2c(p2.E, p2.phi, p2.eta);
            imp1.add(imp2);
            float p = imp1.norm();
            imp1.c2pr(imp1.x, imp1.y, imp1.z);*/
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
            
            h->Fill(invMass(p1,p2));
            //h->Fill(invMass(p1, p2));

        }
    }
  }
  
  /*TGraph *gr = new TGraph(n,x,y);
  gr->Draw("ACP");*/
  
  /*           h->Fit("expo", "0");
            TF1 *fit = h->GetFunction("expo");

  h->Sumw2();
  h->Add(fit, -1);
  h->GetFunction("expo")->Delete();
  
  float rms = h->GetStdDev();
  for(int i = 0; i < n; ++i)
  {
      x[i] = E_min + (E_max-E_min) * float(i)/float(n);
      y[i] = h->GetBinContent(i);
      p[i] = 0.5 * (1 + TMath::Erf(TMath::Abs(y[i])/(sqrt(2.0f) * rms)));
      
      printf("%.2f %.2f\n", x[i], p[i]);
  }
*/
  printf("tight ratio: %.6f\n", float(keptEntries)/float(totalEntries));
  
    TF1 *f1 = new TF1("f1","exp([0]+[1]*x)+[2] * exp((x-124.7)^2 / (2*1.733*1.733))",100,150);
    h->Fit("f1");
    TF1 *fit = h->GetFunction("f1");
    /*float amp = fit->GetParameter(0);
    float avg = fit->GetParameter(1);
    float sigma = fit->GetParameter(2);*/
    
   // printf("%.2f %.2f %.2f %.3f %.3f\n", amp, avg, sigma, fit->GetChisquare(), fit->GetChisquare()/fit->GetNDF());
    h->Draw("E");
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
