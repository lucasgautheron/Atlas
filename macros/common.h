
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
   float p;

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
struct muon : part {} ;


float invMass(part &e1, part &e2)
{
    return TMath::Sqrt(2*e1.E*e2.E*(TMath::CosH(e1.eta-e2.eta)-TMath::Cos(e1.phi-e2.phi)));
}

float invMassTrue(part &p1, part &p2)
{
    return TMath::Sqrt(2*p1.true_E*p2.true_E*(TMath::CosH(p1.true_eta-p2.true_eta)-TMath::Cos(p1.true_phi-p2.true_phi)));
}

float deltaR(part &e1, part &e2)
{
    return sqrt((e1.eta-e2.eta)*(e1.eta-e2.eta) + (e1.phi-e2.phi)*(e1.phi-e2.phi));
}

void bind_attributes(TTree *tree, photon &p1, photon &p2, electron &e1, electron &e2, muon &mu1, muon &mu2)
{
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
    
    tree->SetBranchAddress("el1_pt", &e1.p);
    tree->SetBranchAddress("el2_pt", &e2.p);
    tree->SetBranchAddress("el1_eta", &e1.eta);
    tree->SetBranchAddress("el2_eta", &e2.eta);
    tree->SetBranchAddress("el1_phi", &e1.phi);
    tree->SetBranchAddress("el2_phi", &e2.phi);

    tree->SetBranchAddress("el1_mother", &e1.mother);
    tree->SetBranchAddress("el2_mother", &e2.mother);
    
    tree->SetBranchAddress("mu1_pt", &mu1.p);
    tree->SetBranchAddress("mu2_pt", &mu2.p);
    tree->SetBranchAddress("mu1_eta", &mu1.eta);
    tree->SetBranchAddress("mu2_eta", &mu2.eta);
    tree->SetBranchAddress("mu1_phi", &mu1.phi);
    tree->SetBranchAddress("mu2_phi", &mu2.phi);

    tree->SetBranchAddress("mu1_mother", &mu1.mother);
    tree->SetBranchAddress("mu2_mother", &mu2.mother);
    
    tree->SetBranchAddress("mu1_isTight", &mu1.tight);
    tree->SetBranchAddress("mu2_isTight", &mu2.tight);
    
    tree->SetBranchAddress("mu1_isLoose", &mu1.loose);
    tree->SetBranchAddress("mu2_isLoose", &mu2.loose);
}

