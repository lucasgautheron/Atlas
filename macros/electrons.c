#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "common.h"

const int nx = 10, ny = 10;
float E_min = 40, E_max = 150;
float eta_min = 0, eta_max = 1.2;

void electrons()
{
  TFile* f = TFile::Open("small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron el1, el2;
  muon m1, m2;
  int true_zee = 0, true_zmm = 0;
  
  bind_attributes(tree, p1, p2, el1, el2, m1, m2);

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("mee", &me); 
  tree->SetBranchAddress("is_true_Zee", &true_zee); 
  tree->SetBranchAddress("is_true_Zmm", &true_zmm);

  TH1F* h = new TH1F("h", "\\gamma \\mbox{ efficiency};E_{T} \\mbox{ (GeV)};\\mbox{ratio}", nx, eta_min, eta_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "E",  nx, eta_min, eta_max);
  
    TH1F* e = new TH1F("h", "\\gamma \\mbox{ efficiency};E_{T} \\mbox{ (GeV)};\\mbox{ratio}", ny, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* e2 = new TH1F("h2", "E",  ny, E_min, E_max);
  
  TH2F* h2d = new TH2F("h2d", "\\mbox{Efficacite normalisee};|\\eta|;E_{T};\\mbox{ratio}", nx, eta_min, eta_max, ny, E_min, E_max); 
  TH2F* h2d2 = new TH2F("h2d2", "a;b;c;d", nx, eta_min, eta_max, ny, E_min, E_max); 
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
    
    h2d2->Fill(abs(p1.true_eta), p1.true_E);
    h2d2->Fill(abs(p2.true_eta), p2.true_E);
    h2->Fill(abs(p1.true_eta));
    h2->Fill(abs(p2.true_eta));
    e2->Fill(p1.true_E);
    e2->Fill(p2.true_E);
    
    totalEntries += 2;
    
    if(p1.tight && p1.E > 0)
    {
        h2d->Fill(abs(p1.true_eta), p1.true_E);
        h->Fill(abs(p1.true_eta));
        e->Fill(p1.true_E);
        keptEntries++;
    }
    
    if(p2.tight && p2.E > 0)
    {
        h2d->Fill(abs(p2.true_eta), p2.true_E);
        h->Fill(abs(p2.true_eta));
        e->Fill(p2.true_E);
        keptEntries++;
    }

  }

  h->Sumw2();
  //h2->Sumw2();
  
  h2d->Sumw2();
  h2d2->Sumw2();

  TH1F *r = (TH1F *)h->Clone();
  TH1F *r2d = (TH1F *)h2d->Clone();
  TH1F *re = (TH1F *)e->Clone();
  r->Divide(h, h2, 1., 1., "B");
  re->Divide(e, e2, 1., 1., "B");
  
  r2d->Divide(h2d, h2d2, 1., 1., "B");
  
  float k = float(keptEntries) / float(totalEntries);
  
  for(int x = 1; x <= nx; ++x)
  {
        float v = r->GetBinContent(x);
        for(int y = 1; y <= ny; ++y)
        {
            r2d->SetBinContent(x, y, r2d->GetBinContent(x, y)/v);
        }
  }
  
  for(int y = 1; y <= ny; ++y)
  {
        float v = re->GetBinContent(y);
        for(int x = 1; x <= nx; ++x)
        {
            r2d->SetBinContent(x, y, k*r2d->GetBinContent(x, y)/v);
        }
  }
  
  r2d->Draw("lego2");
  //h2d->Draw("box");

  printf("distrib mgg: %.2f %.2f\n", avg_m, sigma_m);
  printf("loose ratio: %.6f\n", float(keptEntries)/float(totalEntries));
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
