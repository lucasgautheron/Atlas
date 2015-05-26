#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

struct photon
{
   float E, eta, phi;
   float E_true, eta_true, phi_true;
   int tight;
   
};

float inv_mass(photon &p1, photon &p2)
{
    return TMath::Sqrt(2*p1.E*p2.E*(TMath::CosH(p1.eta-p2.eta)-TMath::Cos(p1.phi-p2.phi)));
}

void tests()
{
  TFile* f = TFile::Open("/home/gautheron/Documents/git/Atlas/macros/small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  photon p1, p2;

  float m, ET1, ET2;
  int isTight1, isTight2;
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

  tree->SetBranchAddress("ph1_isTight", &p1.tight);
  tree->SetBranchAddress("ph2_isTight", &p2.tight);

  //TH1F* h = new TH1F("h", "invariant mass", 500, 100, 600); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h = new TH1F("h", "eta", 500, -5, 5);

  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (!p1.tight || !p2.tight) continue;
    if (p1.E < 50 || p2.E < 50) continue; 
    //printf("%.5f %.5f\n", m, inv_mass(p1, p2));
    //if(i%2) printf("%.5f %.5f\n", p1.eta_true, p1.eta);
    h->Fill(p1.eta); // Enter the mass value into the histogram, for the events passing the cuts above
    h->Fill(p2.eta);
  }
  h->Draw(); // plot the histogram
}
