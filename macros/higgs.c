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

void higgs()
{
//  TFile* f = TFile::Open("small.root");  // open the file
//    TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
  //  TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_ZH900_Pythia_NoMassCut.root");  // open the file
 // TTree* tree = (TTree*)f->Get("tree");

  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
 
  const int N = 5;
  TFile *files[N];
  TTree *trees[N];
  TH1D *histograms[N];
  int count[N];
  
 /*      p0           1.05852e+01     fixed    
   2  p1          -2.37115e-02     fixed    
   3  p2           1.80387e+02   2.98922e+01   9.44240e-02   4.40528e-05
   4  p3           1.26228e+02   3.40480e-01   1.07401e-03   2.40111e-03
   5  p4           1.80000e+00     fixed */
  
  double cross_sections[N] = { 8.19759 * 5.51e-5, 4.082711*1.25e-5, 3.191034*2.65e-6, 1.4391*3.11e-7, 0.61988*5e-8 };
  
  const int n = 50;
  
  TH1D *h = new TH1D("h", "\\log{S/B};m_{\\gamma\\gamma}\\mbox{ (GeV)};\\log{S/B}", n, E_min, E_max);
  
  for(int k = 0; k < N; ++k)
  {
      char filename[128] = "";
      sprintf(filename, "Hgg_Moriond2013-Y2012_ZH%d00_Pythia_NoMassCut.root", k+2);
      files[k] = TFile::Open(filename);
      trees[k] = (TTree *)files[k]->Get("tree");
      count[k] = trees[k]->GetEntries();

      bind_attributes(trees[k], p1, p2, e1, e2, m1, m2);
      char h_name[5] = ""; sprintf(h_name, "h_%d", k);
      histograms[k] = new TH1D(h_name, "\\log{S/B};m_{\\gamma\\gamma}\\mbox{ (GeV)};\\log{S/B}", n, E_min, E_max);
      
      cross_sections[k] = 900 * cross_sections[k] / (22.3*2.3e-3);
  }
  
  char formula[64] = "";
  sprintf(formula, "exp(10.5852-0.0237115*x) * %f", 1 * float(E_max-E_min) / float(n));
//  sprintf(formula, "exp(9.89243-0.0237196*x) * %f", 1 * float(E_max-E_min) / float(n));
  printf("%s\n", formula);
  TF1 *background = new TF1("background", formula, E_min, E_max);
  
  h->Add(background);
  
  TCanvas *c = new TCanvas("c");
  c->SetLogy();
   
  for (unsigned a = 0; a < N; ++a)
  {
	  TTree *tree = trees[a];
	  for (unsigned int i = 0; i < count[a]; i++) {
		    tree->GetEntry(i);
            if(p1.E > 0 && p2.E > 0)
            {
                double mgg = (double)invMass(p1, p2);
                //if(!a) printf("%.5f\n", abs(mgg-(double(a+2)*100)));
                if(abs(mgg-(double(a+2)*100)) < 100) histograms[a]->Fill(mgg);
            }
	  }
   }
   
   //background->Draw();
   
  for (unsigned a = 0; a < N; ++a)
  {
  
      printf("factor: %e\n",cross_sections[a]/histograms[a]->Integral());
      histograms[a]->Scale( cross_sections[a]/histograms[a]->Integral());
      histograms[a]->Divide(h);
      printf("factor: %e\n",histograms[a]->Integral());
      
      //histograms[a]->Add(h);
      switch(a)
      {
          case 0: histograms[a]->SetFillColor(kRed); break;
          case 1: histograms[a]->SetFillColor(kGreen); break;
          case 2: histograms[a]->SetFillColor(kYellow); break;
          case 3: histograms[a]->SetFillColor(kOrange); break;
          case 4: histograms[a]->SetFillColor(kBlue); break;
          case 5: histograms[a]->SetFillColor(kPink); break;
      }
      
      histograms[a]->SetLineColor(kBlack);
      
      /*if(a)*/ histograms[a]->Draw("same");
      /*else *///if(!a) histograms[a]->Draw();
  }
  
   /*h->SetLineColor(kBlack);
   h->SetFillColorAlpha(kWhite, 0.2);
   h->Draw("same");*/
}
