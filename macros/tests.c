#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

struct params
{
    bool enable;
    float E[2];
    float eta[2];
};

struct photon
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

float invMass(photon &p1, photon &p2)
{
    return TMath::Sqrt(2*p1.E*p2.E*(TMath::CosH(p1.eta-p2.eta)-TMath::Cos(p1.phi-p2.phi)));
}

float invMassTrue(photon &p1, photon &p2)
{
    return TMath::Sqrt(2*p1.true_E*p2.true_E*(TMath::CosH(p1.true_eta-p2.true_eta)-TMath::Cos(p1.true_phi-p2.true_phi)));
}

float deltaR(photon &p1, photon &p2)
{
    return sqrt((p1.eta-p2.eta)*(p1.eta-p2.eta) + (p1.phi-p2.phi)*(p1.phi-p2.phi));
}

void tests()
{
  TFile* f = TFile::Open("/home/gautheron/Documents/git/Atlas/macros/small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m;
  photon p1, p2;

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("ph1_Et", &p1.E);
  tree->SetBranchAddress("ph2_Et", &p2.E);
  tree->SetBranchAddress("ph1_eta", &p1.eta);
  tree->SetBranchAddress("ph2_eta", &p2.eta);
  tree->SetBranchAddress("ph1_phi", &p1.phi);
  tree->SetBranchAddress("ph2_phi", &p2.phi);

  tree->SetBranchAddress("ph1_mother", &p1.mother);
  tree->SetBranchAddress("ph2_mother", &p2.mother);

  tree->SetBranchAddress("ph1_Et_true", &p1.E_true);
  tree->SetBranchAddress("ph2_Et_true", &p2.E_true);
  tree->SetBranchAddress("ph1_eta_true", &p1.eta_true);
  tree->SetBranchAddress("ph2_eta_true", &p2.eta_true);
  tree->SetBranchAddress("ph1_phi_true", &p1.phi_true);
  tree->SetBranchAddress("ph2_phi_true", &p2.phi_true);

  tree->SetBranchAddress("truePhoton1_Et", &p1.true_E);
  tree->SetBranchAddress("truePhoton2_Et", &p2.true_E);

  tree->SetBranchAddress("truePhoton1_eta", &p1.true_eta);
  tree->SetBranchAddress("truePhoton2_eta", &p2.true_eta);

  tree->SetBranchAddress("truePhoton1_phi", &p1.true_phi);
  tree->SetBranchAddress("truePhoton2_phi", &p2.true_phi);

  tree->SetBranchAddress("truePhoton1_motherType", &p1.true_mother);
  tree->SetBranchAddress("truePhoton2_motherType", &p2.true_mother);

  tree->SetBranchAddress("ph1_isTight", &p1.tight);
  tree->SetBranchAddress("ph2_isTight", &p2.tight);

  tree->SetBranchAddress("ph1_isLoose", &p1.loose);
  tree->SetBranchAddress("ph2_isLoose", &p2.loose);

  TH1F* h = new TH1F("h", "invariant mass", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "invariant mass", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //h2->SetLineColor(kRed);

  int totalEntries = 0;
  int looseEntries = 0;
  int keptEntries = 0;

  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);
    if(p1.mother==25 && p2.mother==25) h->Fill(invMassTrue(p1,p2));
    if(p1.loose && p2.loose) h2->Fill(invMass(p1,p2));
  }

  h->Draw();

  printf("distrib mgg: %.2f %.2f\n", avg_m, sigma_m);
  printf("loose ratio: %.6f\n", float(looseEntries)/float(totalEntries));
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
