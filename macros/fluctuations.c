#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"


bool reject = false;

float min_mass = 0, max_mass = 0;
  const int n = 1400;
  float E_min = 100, E_max = 450;

TH1F* h;

double total_chi = 0;
unsigned int n_d_f = 0;

Double_t expon(Double_t  *x, Double_t  *p)
{
    if(reject && x[0] > min_mass && x[0] < max_mass)
    {
        TF1::RejectPoint();
        return 0;
    }
    return exp(p[0] + p[1] * x[0]);
}

double poly(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e, Double_t x)
{
    return a + b*x + c*x*x + d*x*x*x + e*x*x*x*x;
}

Double_t polyn(Double_t  *x, Double_t  *p)
{
    if(reject && x[0] > min_mass && x[0] < max_mass)
    {
        TF1::RejectPoint();
        return 0;
    }
    return poly(p[0],p[1],p[2],p[3],p[4],x[0]);
}

const int essais = 103;
double deltachi[essais];
double deltachi_poly[essais];
double masses[essais];
TF1 *fitexpo[essais];

double scan_window_poly(int k, float mass, float width)
{

   min_mass = mass-2*width;
   max_mass = mass+2*width;
   
   float min_window = mass-8*width;
   float max_window = mass+8*width;
   
   float expected_ndf = float(n) * ((max_window-min_window) - (max_mass-min_mass))/(E_max-E_min);
   
   fitexpo[k] = new TF1("fitexpo", /*expon*/polyn, min_window, max_window, /*2*/5);
    reject = true;
   TFitResultPtr r = h->Fit(fitexpo[k], "QRS");
   reject = false;
   
   //printf("%.3f %.3f\n",  fitexpo->GetChisquare(), fitexpo->GetChisquare()/fitexpo->GetNDF());
   
   float a = fitexpo[k]->GetParameter(0);
   float b = fitexpo[k]->GetParameter(1);
   
   float c = fitexpo[k]->GetParameter(2);
   float d = fitexpo[k]->GetParameter(3);
   float e = fitexpo[k]->GetParameter(4);
   
    char fname[5] = "";
    sprintf(fname, "f%d", k);
   
  
//    TF1 *f1 = new TF1(fname,"exp([0]+[1]*x)+[2] * exp(-(x-[3])^2 / (2*[4]*[4]))", mass-4*width, mass+4*width);
    TF1 *f1 = new TF1(fname,"[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x+[5] * exp(-(x-[6])^2 / (2*[7]*[7]))", mass-4*width, mass+4*width);
    //f1->SetParLimits(3,120,130);
//    f1->SetParameters(a,b,1,mass,1.8);
    f1->SetParameters(a,b,c,d,e,1,mass,1.8);
    f1->SetParLimits(2+3,0,5000);
    f1->SetParLimits(3+3,mass,mass);
    f1->SetParLimits(4+3,0, 3*width);
    f1->SetParLimits(0, a, a);
    f1->SetParLimits(1, b, b);
    f1->SetParLimits(2, c, c);
    f1->SetParLimits(3, d, d);
    f1->SetParLimits(4, e, e);
    

    h->Fit(fname, "QR");
    TF1 *fit = h->GetFunction(fname);
    float ac = fit->GetParameter(2), ad = fit->GetParameter(3), ae = fit->GetParameter(4);
    
   float chi_bg = 0, chi_tot = 0;
   int NDF = 0;
   for(int j = 1; j <= n; ++j)
   {
   
       float x = E_min + (E_max - E_min) * (float(j))/float(n);
       if(x < mass-2*width || x > mass+2*width) continue;
       ++NDF;
       //float y = exp(a+b*x);
       float y = poly(a,b,c,d,e,x);
       //printf("%.3f, %.3f %.3f\n", x, y, (float)h->GetBinContent(j));
       chi_bg += (h->GetBinContent(j)-y)*(h->GetBinContent(j)-y) / y;
       total_chi += (h->GetBinContent(j)-y)*(h->GetBinContent(j)-y) / (y*y);
       y = exp(a+b*x) + ac * exp(-(x-ad)*(x-ad) / (2*ae*ae));
       chi_tot += (h->GetBinContent(j)-y)*(h->GetBinContent(j)-y) / y;
   }
   
   if(abs(mass-126) > 3)
   {
       n_d_f += r->Ndf();
   }
   
   printf("%.7f %.3f %.3f %.3f %d %.3f %.3f \n", (float)ROOT::Math::chisquared_cdf_c(chi_bg, double((int)NDF -1)), (float)chi_bg, (float)fit->GetNDF(), (float)r->Chi2(), r->Ndf(), expected_ndf, float(NDF));
   
   return ROOT::Math::chisquared_cdf_c(chi_bg, double((int)NDF -1));
}

double scan_window(int k, float mass, float width)
{

   min_mass = mass-2*width;
   max_mass = mass+2*width;
   
   float min_window = mass-8*width;
   float max_window = mass+8*width;
   
   float expected_ndf = float(n) * ((max_window-min_window) - (max_mass-min_mass))/(E_max-E_min);
   
   fitexpo[k] = new TF1("fitexpo", expon, min_window, max_window, 2);
    reject = true;
   TFitResultPtr r = h->Fit(fitexpo[k], "QRS");
   reject = false;
   
   //printf("%.3f %.3f\n",  fitexpo->GetChisquare(), fitexpo->GetChisquare()/fitexpo->GetNDF());
   
   float a = fitexpo[k]->GetParameter(0);
   float b = fitexpo[k]->GetParameter(1);
   
   /*float c = fitexpo[k]->GetParameter(2);
   float d = fitexpo[k]->GetParameter(3);
   float e = fitexpo[k]->GetParameter(4);*/
   
    char fname[5] = "";
    sprintf(fname, "f%d", k);
   
  
    TF1 *f1 = new TF1(fname,"exp([0]+[1]*x)+[2] * exp(-(x-[3])^2 / (2*[4]*[4]))", mass-4*width, mass+4*width);
    //f1->SetParLimits(3,120,130);
//    f1->SetParameters(a,b,1,mass,1.8);
    f1->SetParameters(a,b,1,mass,1.8);
    f1->SetParLimits(2,0,5000);
    f1->SetParLimits(3,mass,mass);
    f1->SetParLimits(4,0, 3*width);
    f1->SetParLimits(0, a, a);
    f1->SetParLimits(1, b, b);
    

    h->Fit(fname, "QR");
    TF1 *fit = h->GetFunction(fname);
    float ac = fit->GetParameter(2), ad = fit->GetParameter(3), ae = fit->GetParameter(4);
    
   float chi_bg = 0, chi_tot = 0;
   int NDF = 0;
   for(int j = 1; j <= n; ++j)
   {
   
       float x = E_min + (E_max - E_min) * (float(j))/float(n);
       if(x < mass-2*width || x > mass+2*width) continue;
       ++NDF;
       float y = exp(a+b*x);
       //printf("%.3f, %.3f %.3f\n", x, y, (float)h->GetBinContent(j));
       chi_bg += (h->GetBinContent(j)-y)*(h->GetBinContent(j)-y) / y;
       total_chi += (h->GetBinContent(j)-y)*(h->GetBinContent(j)-y) / (y*y);
       y = exp(a+b*x) + ac * exp(-(x-ad)*(x-ad) / (2*ae*ae));
       chi_tot += (h->GetBinContent(j)-y)*(h->GetBinContent(j)-y) / y;
   }
   
   if(abs(mass-126) > 3)
   {
       n_d_f += r->Ndf();
   }
   
   printf("%.7f %.3f %.3f %.3f %d %.3f %.3f \n", (float)ROOT::Math::chisquared_cdf_c(chi_bg, double((int)NDF -1)), (float)chi_bg, (float)fit->GetNDF(), (float)r->Chi2(), r->Ndf(), expected_ndf, float(NDF));
   
   return ROOT::Math::chisquared_cdf_c(chi_bg, double((int)NDF -1));
}


void fluctuations()
{
  double expected_mass = 125;
   TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_NoMassCut.root");  // open the file
 // TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_ZH600_Pythia_NoMassCut.root");
    //TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
  //  TFile* f = TFile::Open("small.root");  // open the file
  TTree* tree = (TTree*)f->Get("tree");

  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
  
  bind_attributes(tree, p1, p2, e1, e2, m1, m2);

  float avg_m, sigma_m;

 // tree->SetBranchAddress("mgg", &m); 
 // tree->SetBranchAddress("mee", &me); 
  

  h = new TH1F("h", "m_{\\gamma\\gamma} \\mbox{ distribution};m_{\\gamma\\gamma} \\mbox{ (GeV)};\\mbox{Ev/GeV}", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
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
  
  for(int j = 0; j < essais; ++j)
  {
      //double m = (350-109.2) * float(j)/float(essais) + 109.2;
      double m = (350-110) * float(j)/float(essais) + 110;
      deltachi[j] = scan_window(j, m, 1.4);
      deltachi_poly[j] = scan_window_poly(j, m, 1.4);
      masses[j] = m;
      printf("%.3f %.3f %.3f\n", float(m), float(deltachi[j]), float(deltachi_poly[j]));
  }
  
   TCanvas *c1 = new TCanvas("c1","",200,10,700,500);
  c1->SetLogy();
  c1->SetGrid();
  
  TGraph *gr = new TGraph(essais,masses,deltachi);
  gr->SetTitle("\\log \\mbox{ p-value};m_{\\gamma\\gamma}\\mbox{ (GeV)};\\log p_{0}");
  gr->SetMarkerStyle(2);
  gr->Draw("ACP");
  
        TGraph *gr_poly = new TGraph(essais,masses,deltachi_poly);
  gr_poly->SetTitle("\\log \\mbox{ p-value};m_{\\gamma\\gamma}\\mbox{ (GeV)};\\log p_{0}");
  gr_poly->SetMarkerStyle(2);
  gr_poly->SetMarkerColor(kBlue);
  gr_poly->SetLineStyle(2);
  gr_poly->SetLineColor(kBlue);
  gr_poly->Draw("CP");
 
  c1->Update();
  
printf("%.3f (%.3f %d)", (float)(total_chi/double(n_d_f)), float(total_chi), n_d_f);
     /*  h->SetMarkerStyle(20);
       h->SetMarkerSize(1);
    h->Draw("E");*/
  
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
