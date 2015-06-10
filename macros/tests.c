#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"


bool reject = false;

Double_t expon(Double_t  *x, Double_t  *p)
{
    if(reject && x[0] > 120 && x[0] < 130)
    {
        TF1::RejectPoint();
        return 0;
    }
    return exp(p[0] + p[1] * x[0]);
}

void tests()
{
  double expected_mass = 125;
   TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_NoMassCut.root");  // open the file
 // TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_ZH600_Pythia_NoMassCut.root");
    //TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
  //  TFile* f = TFile::Open("small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
  
  bind_attributes(tree, p1, p2, e1, e2, m1, m2);

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("mee", &me); 
  
  const int n = 150;
  float E_min = 100, E_max = 250;
  TH1F* h = new TH1F("h", "m_{\\gamma\\gamma} \\mbox{ distribution};m_{\\gamma\\gamma} \\mbox{ (GeV)};\\mbox{Ev/GeV}", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //TH1F* h2 = new TH1F("h2", "invariant mass gamma gamma", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //h2->SetLineColor(kRed);
  
 // TF1 *fitgauss = new TF1("fitexpo", "gaus", 86, 96);

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
            
            if(p1.E > 0 && p2.E > 0) h->Fill(invMass(p1,p2));
            //h->Fill(invMass(p1, p2));

        }
    }
    /*if(m1.p > 0 && m2.p > 0)
    {
        m1.E = m1.p;
        m2.E = m2.p;
        h->Fill(invMass(m1, m2));

    }*/
    
    /*if(e1.p > 0 && e2.p > 0)
    {
        e1.E = e1.p;
        e2.E = e2.p;
        h->Fill(invMass(e1, e2));

    }*/
  }
  
  //h->Fit(fitgauss, "R");
  
  /*TGraph *gr = new TGraph(n,x,y);
  gr->Draw("ACP");*/
  
    TF1 *fitexpo = new TF1("fitexpo", expon, E_min, E_max, 2);
    fitexpo->SetLineStyle(2);
    reject = true;
   h->Fit(fitexpo);
   reject = false;
   fitexpo->Draw();
   
   printf("%.3f %.3f\n",  fitexpo->GetChisquare(), fitexpo->GetChisquare()/fitexpo->GetNDF());
   
   float a = fitexpo->GetParameter(0);
   float b = fitexpo->GetParameter(1);
   
  
    TF1 *f1 = new TF1("f1","exp([0]+[1]*x)+[2] * exp(-(x-[3])^2 / (2*[4]))",100,150);
    //f1->SetParLimits(3,120,130);
    f1->SetParameters(a,b,100,expected_mass,1.8);
    //f1->SetParLimits(4,0.5,2*2.0);
    f1->SetParLimits(4,2.0, 2.0);
    f1->SetParLimits(0, a, a);
    f1->SetParLimits(1, b, b);

    h->Fit("f1");
    TF1 *fit = h->GetFunction("f1");
    float c = fit->GetParameter(2), d = fit->GetParameter(3), e = fit->GetParameter(4);
    
   float chi_bg = 0, chi_tot = 0;
   for(int k = 1; k <= n; ++k)
   {
       float x = E_min + (E_max - E_min) * (float(k))/float(n);
       float y = fitexpo->Eval((float)x);
       printf("%.3f %.3f\n", y, (float)h->GetBinContent(k));
       chi_bg += (h->GetBinContent(k)-y)*(h->GetBinContent(k)-y) / y;
       y = fit->Eval((float)x);
       y = exp(a+b*x) + c * exp(-(x-d)*(x-d) / (2*e));
       chi_tot += (h->GetBinContent(k)-y)*(h->GetBinContent(k)-y) / y;
   }
   
   printf("%.3f %.3f %.3f %.3f\n", chi_bg, chi_bg/float(n-1), chi_tot, chi_tot/float(n-1));
    
      printf("%.3f %.3f\n",  fit->GetChisquare(), fit->GetChisquare()/fit->GetNDF());
       h->SetMarkerStyle(20);
       h->SetMarkerSize(1);
    h->Draw("E");
   fitexpo->Draw("same");
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
