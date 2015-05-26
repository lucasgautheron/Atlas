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

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("ph1_Et", &p1.E);
  tree->SetBranchAddress("ph2_Et", &p2.E);
  tree->SetBranchAddress("ph1_eta", &p1.eta);
  tree->SetBranchAddress("ph2_eta", &p2.eta);
  tree->SetBranchAddress("ph1_phi", &p1.phi);
  tree->SetBranchAddress("ph2_phi", &p2.phi);

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

  tree->SetBranchAddress("ph1_isTight", &p1.tight);
  tree->SetBranchAddress("ph2_isTight", &p2.tight);

  tree->SetBranchAddress("ph1_isLoose", &p1.loose);
  tree->SetBranchAddress("ph2_isLoose", &p2.loose);

  // TH1F* h = new TH1F("h", "invariant mass", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h = new TH1F("h", "eta", 500, -5, 5);
  //TH1F* h2 = new TH1F("h2", "eta_true", 500, -5, 5);
  
  //h2->SetLineColor(kRed);

  int totalEntries = tree->GetEntries();
  int looseEntries = 0;

  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (p1.loose && p2.loose) ++looseEntries;
    if (!p1.tight || !p2.tight) continue;
    //if (p1.E < 50 || p2.E < 50) continue; 
    //printf("%.5f %.5f\n", m, inv_mass(p1, p2));
    //if(deltaR(p1, p2) < 1) printf("%.5f %s %s\n", deltaR(p1, p2), p1.tight ? "t" : "n", p2.tight ? "t" : "n");
    h->Fill(p1.eta_true); // Enter the mass value into the histogram, for the events passing the cuts above
    //h2->Fill(p1.eta);
  }
  printf("%.2f", float(looseEntries)/float(totalEntries));
  h->Draw(); // plot the histogram
  //h2->Draw();
}
