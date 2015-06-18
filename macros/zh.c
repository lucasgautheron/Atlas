#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"


bool reject = false;

Double_t expon(Double_t  *x, Double_t  *p)
{
    if(reject && x[0] > 80 && x[0] < 100)
    {
        TF1::RejectPoint();
        return 0;
    }
    return p[0];
}

TH1F * DrawOverflow(TH1F *h)
{
      // This function paint the histogram h with an extra bin for overflows
   UInt_t nx    = h->GetNbinsX()+1;
   Double_t *xbins= new Double_t[nx+1];
   for (UInt_t i=0;i<nx;i++)
     xbins[i]=h->GetBinLowEdge(i+1);
   xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
   char *tempName= new char[strlen(h->GetName())+10];
   sprintf(tempName,"%swtOverFlow",h->GetName());
   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
   // Reset the axis labels
   htmp->SetXTitle(h->GetXaxis()->GetTitle());
   htmp->SetYTitle(h->GetYaxis()->GetTitle());
   // Fill the new hitogram including the extra bin for overflows
   for (UInt_t i=1; i<=nx; i++)
     htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
   // Fill the underflows
   htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   // Restore the number of entries
   htmp->SetEntries(h->GetEntries());
   // FillStyle and color
   htmp->SetFillStyle(h->GetFillStyle());
   htmp->SetFillColor(h->GetFillColor());
   return htmp;
}

void zh()
{
  double expected_mass = 91.2;
   //TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_NoMassCut.root");  // open the file
 // TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_ZH600_Pythia_NoMassCut.root");
    TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
   // TFile* f = TFile::Open("small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
  electron th;
  bool dielectron = 0, dimuon = 0;
  int zee, zmm;
  
  bind_attributes(tree, p1, p2, e1, e2, m1, m2);

  float avg_m, sigma_m;

  tree->SetBranchAddress("mgg", &m); 
  tree->SetBranchAddress("mee", &me); 
  tree->SetBranchAddress("isDielectron", &dielectron);
  tree->SetBranchAddress("isMuon", &dimuon);
  tree->SetBranchAddress("trueHiggs_pt", &th.p);
  tree->SetBranchAddress("is_true_Zee", &zee);
  tree->SetBranchAddress("is_true_Zmm", &zmm);
  
  
  const int n = 50;
  float E_min = 150, E_max = 1250;
  TH1F* h = new TH1F("h", "m_{ZH} \\mbox{ distribution};m_{ZH} \\mbox{ (GeV)};\\mbox{Ev/10GeV}", n, E_min, E_max); // create a histogram : 500 bins
  TH1F* h2 = new TH1F("h2", "m_{ZH} \\mbox{ distribution};m_{ZH} \\mbox{ (GeV)};\\mbox{Ev/5GeV}", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //TH1F* h2 = new TH1F("h2", "invariant mass gamma gamma", 500, 50, 250); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //h2->SetLineColor(kRed);
  
 //TF1 *fitgauss = new TF1("fitexpo", "gaus", 70, 110);

  int totalEntries = 0;
  int looseEntries = 0;
  int keptEntries = 0;
  
  unsigned int count = tree->GetEntries();
  
  double max_m = 0;
  
  int found_zee = 0, found_zmm = 0, true_zee = 0, true_zmm = 0;
  int found_hgg = 0;

  for (unsigned int i = 0; i < count; i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);
    
    if(p1.tight && p2.tight && p1.E > 0 && p2.E > 0) found_hgg++;

    //if(p1.mother==25 && p2.mother==25)
    {
        //if(p1.tight && p2.tight)
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
            
            double z_mass = -1, z_energy, z_momentum;
            electron z;
            vec imp1, imp2;
            vec a, b, c, d;
            float ae, be, ce, de;
            
            bool has_z_and_h = true;
            
            if(zee) true_zee++;
            else if(zmm) true_zmm++;
            
            if(m1.p > 0 && m2.p > 0 /*&& dimuon*/)
            {
                m1.E = m1.p;
                m2.E = m2.p;
                z_mass = invMass(m1,m2);
                
                imp1.pr2c(m1.E, m1.phi, m1.eta);
                imp2.pr2c(m2.E, m2.phi, m2.eta);
                imp1.add(imp2);
                float p = imp1.norm();
                z_momentum = p;
                z_energy =  m1.E * cosh(m1.eta) + m2.E * cosh(m2.eta);
                imp1.c2pr(imp1.x, imp1.y, imp1.z);
                z.E = p;
                z.phi = imp1.y;
                z.eta = imp1.z;
                found_zmm++;
                
                a.pr2c(m1.E, m1.phi, m1.eta);
                b.pr2c(m2.E, m2.phi, m2.eta);
                ae = m1.E * cosh(m1.eta);
                be = m2.E * cosh(m2.eta);
            }
            else if(e1.p > 0 && e2.p > 0)
            {
                e1.E = e1.p;
                e2.E = e2.p;
                z_mass = invMass(e1,e2);
                
                imp1.pr2c(e1.E, e1.phi, e1.eta);
                imp2.pr2c(e2.E, e2.phi, e2.eta);
                imp1.add(imp2);
                float p = imp1.norm();
                z_momentum = p;
                z_energy = e1.E * cosh(e1.eta) + e2.E * cosh(e2.eta);
                imp1.c2pr(imp1.x, imp1.y, imp1.z);
                z.E = p;
                z.phi = imp1.y;
                z.eta = imp1.z;
                found_zee++;
                
                a.pr2c(e1.E, e1.phi, e1.eta);
                b.pr2c(e2.E, e2.phi, e2.eta);
                ae = e1.E * cosh(e1.eta);
                be = e2.E * cosh(e2.eta);
            }
            
            if(z_mass < 0) continue;
            has_z_and_h &= (z_mass < 100 && z_mass > 80);
           
            
            if(p1.tight && p2.tight && p1.E > 0 && p2.E > 0)
            {
                double h_mass = -1, h_energy, h_momentum;
                h_mass = invMass(p1,p2);
                imp1.pr2c(p1.E, p1.phi, p1.eta);
                imp2.pr2c(p2.E, p2.phi, p2.eta);
                
                c.pr2c(p1.E, p1.phi, p1.eta);
                d.pr2c(p2.E, p2.phi, p2.eta);
                
                ce =  p1.E * cosh(p1.eta);
                de =  p2.E * cosh(p2.eta);
                
                electron H;
                imp1.add(imp2);
                float p = imp1.norm();
                h_momentum = p;
                h_energy = p1.E * cosh(p1.eta) + p2.E * cosh(p2.eta);
                imp1.c2pr(imp1.x, imp1.y, imp1.z);
                H.E = p;
                H.phi = imp1.y;
                H.eta = imp1.z;
                
                imp1.pr2c(z_momentum/cosh(z.eta), z.phi, z.eta);
                imp2.pr2c(h_momentum/cosh(H.eta), H.phi, H.eta);
                imp1.add(imp2);
                float p2 = imp1.norm();
                float m = sqrt((z_energy+h_energy)*(z_energy+h_energy)-p2*p2);
                
                a.add(b);
                a.add(c);
                a.add(d);
                
                m = sqrt((ae+be+ce+de)*(ae+be+ce+de) - a.norm() * a.norm());
                
                
                has_z_and_h &= (h_mass < 130 && h_mass > 120);
                
                //has_z_and_h &= (zee || zmm);
                
                if(m > 0)
                {
                    if(has_z_and_h) printf("%.3f %.3f\n", z_mass+h_mass, m);
                	h->Fill(m);
	                if(!has_z_and_h) h2->Fill(m);
	                //if(has_z_and_h) h->Fill(z_mass+h_mass);
	            }
	            
	            if(m > max_m /*&& has_z_and_h*/) max_m = m;
            }

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
  
  printf("Masse max: %.3f; Int: %.3f\n, ee %.4f, mm %.4f, %d, %d, %d, %d", max_m, h->Integral(), float(found_zee)/float(true_zee), float(found_zmm)/float(true_zmm), true_zee, true_zmm, found_hgg, (int)count);
  
  //h->Fit(fitgauss, "R");
  
  /*TGraph *gr = new TGraph(n,x,y);
  gr->Draw("ACP");*/
  
  /* TF1 *fitexpo = new TF1("fitexpo", expon, E_min, E_max, 2);
    fitexpo->SetLineStyle(2);
    reject = true;
   h->Fit(fitexpo);
   reject = false;
   fitexpo->Draw();*/
   
 /*  printf("%.3f %.3f\n",  fitexpo->GetChisquare(), fitexpo->GetChisquare()/fitexpo->GetNDF());
   
   float a = fitexpo->GetParameter(0);
   float b = fitexpo->GetParameter(1);
   
  
    TF1 *f1 = new TF1("f1","[0]+[1] * exp(-(x-[2])^2 / (2*[3]))",100,150);
    //f1->SetParLimits(3,120,130);
    f1->SetParameters(a,50,expected_mass,2);
    //f1->SetParLimits(4,0.5,2*2.0);
    f1->SetParLimits(3,3.0, 10.0);
    f1->SetParLimits(0, a, a);

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
    
      printf("%.3f %.3f\n",  fit->GetChisquare(), fit->GetChisquare()/fit->GetNDF());*/
       h->SetMarkerStyle(20);
       h->SetMarkerSize(1);
       h2->SetFillColor(kRed);
       
           h->Draw("P");
           h2->Draw("same");
   //fitexpo->Draw("same");
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
