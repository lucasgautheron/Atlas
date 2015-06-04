#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TGraph.h"
#include "common.h"
#include "maths.h"

  float E_min = 150, E_max = 500;

bool keep_photon(photon &p)
{
    return p.tight && p.E > E_min && p.E < E_max && p.convFlag <= 0 && abs(p.E-p.true_E) < 50;
}

void eta()
{
//  TFile* f = TFile::Open("small.root");  // open the file
//    TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_merge_200804_216432_Presel_1lepton.root");
  //  TFile* f = TFile::Open("Hgg_Moriond2013-Y2012_ZH900_Pythia_NoMassCut.root");  // open the file
 // TTree* tree = (TTree*)f->Get("tree");

  float m, me;
  photon p1, p2;
  electron e1, e2;
  muon m1, m2;
 

  float avg_m, sigma_m;
  const int N = 8;
  TFile *files[N];
  TTree *trees[N];
  int count[N];
  
  for(int k = 0; k < N; ++k)
  {
      char filename[128] = "";
      sprintf(filename, "Hgg_Moriond2013-Y2012_ZH%d00_Pythia_NoMassCut.root", k+2);
      files[k] = TFile::Open(filename);
      trees[k] = (TTree *)files[k]->Get("tree");
      count[k] = trees[k]->GetEntries();
      bind_attributes(trees[k], p1, p2, e1, e2, m1, m2);
  }
  
  const int n = 15;
  
  TF1 *fa = new TF1("fa1","(cos(x/2)/sin(x/2)  ) / ((cos(x/2)/sin(x/2)) *(cos(x/2)/sin(x/2)) + 1)",0,3.142);

  TH1F* h = new TH1F("h", "\\mbox{ erreur energie };\\Delta E_{T}/E_{T};\\mbox{count}", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  TH1F* h2 = new TH1F("h2", "invariant mass gamma gamma", n, E_min, E_max); // create a histogram : 500 bins ranging from 100 to 600 GeV.
  //h2->SetLineColor(kRed);


  int totalEntries = 0;
  int looseEntries = 0;
  int keptEntries = 0;
  
  /*int n = 0, k = 0;
  
  for (unsigned int i = 0; i < tree->GetEntries(); i++) {
    //if(p1.true_E < 0) continue;
    tree->GetEntry(i);
    if(p1.true_mother==25 && p2.true_mother==25)
    {
        if(p1.true_E > 0)
        {
            ++n;
        }
    }
  }
  
  float *x = new float[n+1], *y = new float[n+1];*/
 
   //TF1 *f1 = new TF1("f1","gaus",122, 128);
   
   double *err[n];
   int elems[n], index[n];
   double en[n], den[n], xerr[n], yerr[n];
  for (unsigned a = 0; a < N; ++a)
  {
	  TTree *tree = trees[a];
	  for (unsigned int i = 0; i < count[a]; i++) {
		//if(p1.true_E < 0) continue;
		tree->GetEntry(i);

		//if(p1.mother==25 && p2.mother==25)
		{
		    //if(p1.true_E > 0)
		    {
		        vec imp1, imp2;
		        imp1.pr2c(p1.true_E, p1.true_phi, p1.true_eta);
		        imp2.pr2c(p2.true_E, p2.true_phi, p2.true_eta);
		        imp1.add(imp2);
		        float p = imp1.norm();
		        imp1.c2pr(imp1.x, imp1.y, imp1.z);
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
		        // todo: par unité d'angle solide en fct de theta ?
		        //h->Fill(2*atan(exp(-imp1.z))/*+TMath::Pi()/2.0f*/);
		        /*float i1 = p1.true_E * cosh(p1.true_eta);
		        float i2 = p2.true_E * cosh(p2.true_eta);*/
		        
		        
		        /*if(p1.tight && p1.true_E > E_min && p1.true_E < E_max)
		        {
		            h->Fill(p1.true_E);   
		        }


		        if(p2.tight && p2.true_E > E_min && p2.true_E < E_max)
		        {
		            h->Fill(p1.true_E);   
		        }*/
		        
		        p1.E *= cosh(p1.true_eta);
		        p2.E *= cosh(p2.true_eta);
		        p1.true_E *= cosh(p1.true_eta);
		        p2.true_E *= cosh(p2.true_eta);
		        
		        if(keep_photon(p1))
		        {
		           h->Fill(p1.E);
		           if(abs(p1.E-p1.true_E) > 50) printf("%.3f %.3f %.3f %.3f\n", p1.true_E, p1.E, p1.eta, p1.isolTrack);
		        }
		        
		        if(keep_photon(p2))
		        {
		           h->Fill(p2.E);
		           if(abs(p2.E-p2.true_E) > 50) printf("%.3f %.3f %.3f %.3f\n", p2.true_E, p2.E, p2.eta, p2.isolTrack);
		        }
		    }
		}
	  }
   }
  
  for(int k = 0; k < n; ++k)
  {
      err[k] = NULL;
      elems[k] = 0;
  }
  
  for (unsigned a = 0; a < N; ++a)
  {
	  TTree *tree = trees[a];
	  for(unsigned int i = 0; i < count[a]; ++i)
	  {
		  tree->GetEntry(i);
		  TAxis *xaxis = h->GetXaxis();
		  int bin1 = xaxis->FindBin(p1.true_E)-1;
		  int bin2 = xaxis->FindBin(p2.true_E)-1;

		  if(keep_photon(p1) && bin1 >= 0 && bin1<n)
		  {
		      elems[bin1]++;
		  }
		  
		  if(keep_photon(p2) && bin2 >= 0 && bin2<n)
		  {
		      elems[bin2]++;
		  }
	  }
   }
  
  for(int k = 0; k < n; ++k)
  {
      err[k] = new double[elems[k]];
      index[k] = 0;
  }
  
  
  for (unsigned a = 0; a < N; ++a)
  {
	  TTree *tree = trees[a];
	  for(unsigned int i = 0; i < count[a]; ++i)
	  {
		  tree->GetEntry(i);
		  TAxis *xaxis = h->GetXaxis();
		  int bin1 = xaxis->FindBin(p1.true_E)-1;
		  int bin2 = xaxis->FindBin(p2.true_E)-1;
		  
		  

		  if(keep_photon(p1) && bin1 >= 0 &&  bin1 < n)
		  {
		      err[bin1][index[bin1]++] = p1.E-p1.true_E;
		  }
		  
		  if(keep_photon(p2) && bin2 >= 0 &&  bin2 < n)
		  {
		      err[bin2][index[bin2]++] = p2.E-p2.true_E;
		  }
	  }
   }
  
  for(int k = 0; k < n; ++k)
  {
      double mean = 0;
      for(int j = 0; j < elems[k]; ++j)
      {
         mean += err[k][j];
      }
      mean /= double(elems[k]);
      
      double tot = 0;
      //mean = 0;
      for(int j = 0; j < elems[k]; ++j)
      {
         //if(k == 2) printf("err: %.3f \n", err[k][j]);
         tot += (err[k][j]-mean) * (err[k][j]-mean);
      }
      tot /= double(elems[k]);
      double E = float(E_min + float(k)/float(n) * (E_max-E_min));
      /*printf("%.2f %.6f %.6f %.6f %% %.7f %% %.7f %d\n", E, float(mean), float(sqrt(tot)), 100 * float(sqrt(tot))/float(E), 100 * float(sqrt(tot))/float(E)/sqrt(index[k]-1), float(1)/float(n) * (E_max-E_min)/2.0, index[k]);*/
      printf("%.2f %.6f %% %.7f %% %.7f %d\n", E, 100 * float(sqrt(tot))/float(E), 100 * float(sqrt(tot))/float(E)/sqrt(index[k]-1), float(1)/float(n) * (E_max-E_min)/2.0, index[k]);
      en[k] = E;
      den[k] = 100 * float(sqrt(tot))/float(E);
      yerr[k] = 100 * float(sqrt(tot))/float(E)/sqrt(index[k]-1);
      xerr[k] = float(1)/float(n) * (E_max-E_min)/2.0;
  }
  
   TF1 *energy = new TF1("energy", " [0]/sqrt(x) + [1]", 100, 1000);
   TGraphErrors *gr = new TGraphErrors(n,en,den,xerr,yerr);
   gr->Fit("energy");
   gr->SetTitle("\\mbox{Resolution en energie};E\\mbox{ (GeV) };\\sigma_{E} / \\Delta E \\mbox{ \% }");
   gr->SetLineStyle(0);
   gr->SetLineWidth(0);
   /* h->Sumw2();
  h2->Sumw2();
  
  h->Fit("f1");

  TH1F *r = (TH1F *)h->Clone();
  r->Divide(h, h2, 1., 1., "B");
  //r->Draw();
  h->Draw();*/
  gr->Draw("AP");
  

  printf("tight ratio: %.6f\n", float(keptEntries)/float(totalEntries));
    /*h->Multiply(fa);
    h->Draw("E");*/
 // h->Draw(); // plot the histogram
  //h2->Draw();
}
