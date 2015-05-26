#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

struct params
{
    bool enable;
    float E[2];
    float eta[2];
};

struct part
{
   float E, eta, phi;
   float E_true, eta_true, phi_true;
   float true_E, true_eta, true_phi; // ??
   int tight, loose;
   int mother, true_mother;

   bool inBounds(params &p)
   {
       if(!p.enable) return true;
       return (E > p.E[0] && E < p.E[1])
 &&
              (abs(eta) > p.eta[0] && abs(eta) < p.eta[1]);
   }
   
};

struct photon : part {} ;
struct electron : part {} ;


float invMass(photon &e1, photon &e2)
{
    return TMath::Sqrt(2*e1.E*e2.E*(TMath::CosH(e1.eta-e2.eta)-TMath::Cos(e1.phi-e2.phi)));
}

float deltaR(photon &e1, photon &e2)
{
    return sqrt((e1.eta-e2.eta)*(e1.eta-e2.eta) + (e1.phi-e2.phi)*(e1.phi-e2.phi));
}

void electron()
{
  TFile* f = TFile::Open("/home/gautheron/Documents/git/Atlas/macros/small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m;
  electron e1, e2;

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("el1_Et", &e1.E);
  tree->SetBranchAddress("el2_Et", &e2.E);
  tree->SetBranchAddress("el1_eta", &e1.eta);
  tree->SetBranchAddress("el2_eta", &e2.eta);
  tree->SetBranchAddress("el1_phi", &e1.phi);
  tree->SetBranchAddress("el2_phi", &e2.phi);

  tree->SetBranchAddress("el1_mother", &e1.mother);
  tree->SetBranchAddress("el2_mother", &e2.mother);

  tree->SetBranchAddress("el1_Et_true", &e1.E_true);
  tree->SetBranchAddress("el2_Et_true", &e2.E_true);
  tree->SetBranchAddress("el1_eta_true", &e1.eta_true);
  tree->SetBranchAddress("el2_eta_true", &e2.eta_true);
  tree->SetBranchAddress("el1_phi_true", &e1.phi_true);
  tree->SetBranchAddress("el2_phi_true", &e2.phi_true);

  tree->SetBranchAddress("truePhoton1_Et", &e1.true_E);
  tree->SetBranchAddress("truePhoton2_Et", &e2.true_E);

  tree->SetBranchAddress("truePhoton1_eta", &e1.true_eta);
  tree->SetBranchAddress("truePhoton2_eta", &e2.true_eta);

  tree->SetBranchAddress("truePhoton1_phi", &e1.true_phi);
  tree->SetBranchAddress("truePhoton2_phi", &e2.true_phi);

  tree->SetBranchAddress("truePhoton1_motherType", &e1.true_mother);
  tree->SetBranchAddress("truePhoton2_motherType", &e2.true_mother);

  tree->SetBranchAddress("el1_isTight", &e1.tight);
  tree->SetBranchAddress("el2_isTight", &e2.tight);

  tree->SetBranchAddress("el1_isLoose", &e1.loose);
  tree->SetBranchAddress("el2_isLoose", &e2.loose);

  TH1F* h = new TH1F("h", "E", 500, 25, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "E", 500, 25, 250);
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

    if(e1.true_mother > 22)
    {
        h2->Fill(e1.true_E);
    }
    if(e2.true_mother > 22)
    {
        h2->Fill(e2.true_E);
    }
    if (!e1.loose || !e2.loose) continue;
    if(e1.mother > 22)
    {
        h->Fill(e1.true_E);
    }
    if(e2.mother > 22)
    {
        h->Fill(e2.true_E);
    }

  }

  h->Sumw2();
  h2->Sumw2();

  TH1F *r = (TH1F *)h->Clone();
  r->Divide(h, h2, 1., 1., "B");
  r->Draw();

  printf("distrib mgg: %.2f %.2f\n", avg_m, sigma_m);
  printf("loose ratio: %.6f\n", float(looseEntries)/float(totalEntries));
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
