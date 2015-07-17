#include <TMath.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include <TVector.h>
#include "AtlasStyle.C"

#define cos(x) (TMath::Cos((x)))
#define sin(x) (TMath::Sin((x)))
#define cosh(x) (TMath::CosH((x)))
#define sinh(x) (TMath::SinH((x)))
#define atanh(x) (TMath::ATanH((x)))
#define sqrt(x) (TMath::Sqrt((x)))
#define atan(x) (TMath::ATan((x)))
#define acos(x) (TMath::ACos((x)))
#define atan2(x,y) (TMath::ATan2((x),(y)))
#define exp(x) (TMath::Exp((x)))
#define pi TMath::Pi
#define ABS(x) ((x) > 0 ? (x) : -(x))

#define SIGMA 5.67e-8

struct Rectangle
{
    TVector3 pos;
    TVector3 A, B, C, D;
    TVector3 normal;
    double a,b,c,d;
    double width, length;
    TMatrix transform;
    static const int nx = 40, ny = 20;
    
    
    double temperature;
    double emissivity;
    double power;
};

struct Cavity
{
    double radius, length;
    static const int theta_segments = 150, r_segments = 150, z_segments = 150;
    static const int rays = 200;
    
    double emissivity;
    double temperature;
    double total_power;
};

void calculate_power(Cavity *c, Rectangle *r)
{
    r->power = 0;
    double total_power = 0;
    
    r->normal = r->A.Cross(r->B);
    r->normal *= 1/(r->normal.Mag());
    c->total_power = 0;
            
    // "barrel"
    for(int i = 0; i < Cavity::theta_segments; ++i)
    {
        double theta = 2*pi() * double(i) / double(Cavity::theta_segments);
        for(int j = 0; j < Cavity::z_segments; ++j)
        {
            double z = -0.5 * c->length + c->length * double(j) / double(Cavity::z_segments);
            TVector3 origin;
            origin.SetX(c->radius * cos(theta));
            origin.SetY(c->radius * sin(theta));
            origin.SetZ(z);
            
            /*TVector3 OA, OB, OC, OD;
            OA = r->A - origin;
            OB = r->B - origin;
            OC = r->C - origin;
            OD = r->D - origin;
            
            double alpha = OC.Angle(OB);
            double beta = OC.Angle(OA);
            double gamma = OA.Angle(OB);
            
            double beta_p = OC.Angle(OD);
            double gamma_p = OD.Angle(OB);
            
            double A = acos((cos(alpha) - cos(beta) * cos(gamma))/(sin(beta) * sin(gamma)));
            double B = acos((cos(beta) - cos(gamma) * cos(alpha))/(sin(gamma) * sin(alpha)));
            double C = acos((cos(gamma) - cos(beta) * cos(alpha))/ (sin(alpha) * sin(beta)));
            
            if(abs(sin(beta) * sin(gamma) * sin(alpha)) < 0.00001) A = B = C = pi()/3.0;
            
            double A_p = acos((cos(alpha) - cos(beta_p) * cos(gamma_p))/(sin(beta_p) * sin(gamma_p)));
            double B_p = acos((cos(beta_p) - cos(gamma_p) * cos(alpha))/(sin(gamma_p) * sin(alpha)));
            double C_p = acos((cos(gamma_p) - cos(beta_p) * cos(alpha))/ (sin(alpha) * sin(beta_p)));
            
            if(abs(sin(beta_p) * sin(gamma_p) * sin(alpha)) < 0.00001) A_p = B_p = C_p = pi()/3.0;*/
            
            double surface = c->radius * (c->length / double(Cavity::z_segments)) * 2*pi() / double(Cavity::theta_segments);
            double solid = 0;
            
            double dx = r->width / double(Rectangle::nx);
            double dy = r->length / double(Rectangle::ny);
            
            for(int k = 0; k < Rectangle::nx; ++k)
            {
                for(int l = 0; l < Rectangle::ny; ++l)
                {
                    TVector3 pos;
                    pos.SetX(-r->width/2.0 + dx * double(k));
                    pos.SetY(-r->length/2.0 + dy * double(l));
                    pos.SetZ(0);
                    TVector3 diff = pos - origin;
                    TVector3 n1, n2;
                    n2.SetZ(1); n2.SetX(0); n2.SetY(0);
                    n1.SetZ(0); n1.SetX(-cos(theta)); n1.SetY(-sin(theta));
                    
                    double cos1 = cos(n1.Angle(diff)), cos2 = cos(n2.Angle(diff));
                    
                    solid += abs(cos1 * cos2 * surface* dx*dy / (pi() * diff.Mag() * diff.Mag()));
                }
            }
           
            
            //double solid = ABS(A+B+C-pi()) + ABS(A_p+B_p+C_p-pi());
           // printf("%f %f %f %f %f\n", solid, ABS(A+B+C-pi()) + ABS(A_p+B_p+C_p-pi()), origin.X(), origin.Y(), origin.Z());
            //double solid = 4*pi() * hits/(double(Cavity::rays) * double(Cavity::rays));
            
            //printf("%f %f %f\n", (ABS(A+B+C-pi()) + ABS(A_p+B_p+C_p-pi())), solid, r->width * r->length / (origin.Mag() * origin.Mag()));
            
            double power_emitted = SIGMA * pow(c->temperature, 4) * c->emissivity;
            
              //      if(ABS(r->power) < 100)    printf("%f %f %f %f %f %f %f %f\n",  alpha, beta, gamma, r->width * r->length  / (4*pi() * origin.Mag() * origin.Mag()), origin.X(), origin.Y(), origin.Z(), r->power);
            
            r->power += power_emitted * solid;
            
            c->total_power += power_emitted;
        }
    }
    
    // "disks"
    
    for(int k = 0; k < 2; ++k)
    for(int i = 0; i < Cavity::theta_segments; ++i)
    {
        double theta = 2*pi() * double(i) / double(Cavity::theta_segments);
        for(int j = 0; j < Cavity::r_segments; ++j)
        {
            double rad = c->radius* double(j) / double(Cavity::r_segments);
            TVector3 origin;
            origin.SetX(rad * cos(theta));
            origin.SetY(rad * sin(theta));
            origin.SetZ(k ? -c->length/2.0 : +c->length/2.0);
            
           /* TVector3 OA, OB, OC, OD;
            OA = r->A - origin;
            OB = r->B - origin;
            OC = r->C - origin;
            OD = r->D - origin;
            
            double alpha = OC.Angle(OB);
            double beta = OC.Angle(OA);
            double gamma = OA.Angle(OB);
            
            double beta_p = OC.Angle(OD);
            double gamma_p = OD.Angle(OB);
            
            double A = acos((cos(alpha) - cos(beta) * cos(gamma))/(sin(beta) * sin(gamma)));
            double B = acos((cos(beta) - cos(gamma) * cos(alpha))/(sin(gamma) * sin(alpha)));
            double C = acos((cos(gamma) - cos(beta) * cos(alpha))/ (sin(alpha) * sin(beta)));
            
            if(abs(sin(beta) * sin(gamma) * sin(alpha)) < 0.00001) A = B = C = pi()/3.0;
            
            double A_p = acos((cos(alpha) - cos(beta_p) * cos(gamma_p))/(sin(beta_p) * sin(gamma_p)));
            double B_p = acos((cos(beta_p) - cos(gamma_p) * cos(alpha))/(sin(gamma_p) * sin(alpha)));
            double C_p = acos((cos(gamma_p) - cos(beta_p) * cos(alpha))/ (sin(alpha) * sin(beta_p)));
            
            if(abs(sin(beta_p) * sin(gamma_p) * sin(alpha)) < 0.00001) A_p = B_p = C_p = pi()/3.0;*/
            
           double solid = 0;
            
            double dx = r->width / double(Rectangle::nx);
            double dy = r->length / double(Rectangle::ny);
            
            double surface = rad * (c->radius / double(Cavity::r_segments)) * 2*pi() / double(Cavity::theta_segments);
            
            for(int k = 0; k < Rectangle::nx; ++k)
            {
                for(int l = 0; l < Rectangle::ny; ++l)
                {
                    TVector3 pos;
                    pos.SetX(-r->width/2.0 + dx * double(k));
                    pos.SetY(-r->length/2.0 + dy * double(l));
                    pos.SetZ(0);
                    TVector3 diff = pos - origin;
                    TVector3 n1, n2;
                    n2.SetZ(1); n2.SetX(0); n2.SetY(0);
                    n1.SetZ(1); n1.SetX(0); n1.SetY(0);
                    
                    double cos1 = cos(n1.Angle(diff)), cos2 = cos(n2.Angle(diff));
                    
                    solid += abs(cos1 * cos2 * surface* dx*dy / (pi() * diff.Mag() * diff.Mag()));
                }
            }
            //printf("%f %f %f %f %f\n", solid, 4*pi() * hits/2500.0, origin.X(), origin.Y(), origin.Z());
            //double solid = 4*pi() * hits/(double(Cavity::rays) * double(Cavity::rays));

            double power_emitted = SIGMA * pow(c->temperature, 4) * c->emissivity;
            
                   //     printf("%f %f %f %f %f %f\n",  ABS(A+B+C-pi())+ABS(A_p+B_p+C_p-pi()), r->width * r->length  / (4*pi() * origin.Mag() * origin.Mag()), origin.X(), origin.Y(), origin.Z(), r->power);
            
            r->power += power_emitted * solid;
            
            c->total_power += power_emitted;
        }
    }
    //printf("Total power: %f W\n", c->total_power / 1e6);
}

void radiation_2()
{
    Rectangle *r = new Rectangle();
    r->width = 40;
    r->length = 20;
    r->power = 0;
    r->A = *new TVector3(-r->width/2.0, -r->length/2.0, 0);
    r->B = *new TVector3(r->width/2.0, -r->length/2.0, 0);
    r->C = *new TVector3(r->width/2.0, r->length/2.0, 0);
    r->D = *new TVector3(-r->width/2.0, +r->length/2.0, 0);
    r->temperature = 273.15 - 30;
    
    Rectangle *support = new Rectangle();
    support->width = 30;
    support->length = 150;
    support->power = 0;
    support->A = *new TVector3(-r->width/2.0, 0, -r->length/2.0);
    support->B = *new TVector3(r->width/2.0, 0, -r->length/2.0);
    support->C = *new TVector3(r->width/2.0, 0, r->length/2.0);
    support->D = *new TVector3(-r->width/2.0, 0, r->length/2.0);
    support->temperature = 273.15 - 35;
    
    Cavity *c = new Cavity();
    c->radius = 75;
    c->length = 500;
    c->emissivity = 1;
    c->temperature = 273.15 + 26;
    //calculate_power(c, r);
    //calculate_power(c, support);
    printf("recv = %e W, approx = %e W, loss = %e W\n", r->power/1e6, 2 * SIGMA * pow(c->temperature, 4) * r->width * r->length / 1e6, 2 * SIGMA * pow(r->temperature, 4) * r->width * r->length / 1e6);
    support->power = 0;
    
    const int num_temp = 10;
    float min_temp = -40, max_temp = 30;
    float x[num_temp], y[num_temp], z[num_temp];
    
    for(int k = 0; k < num_temp; ++k)
    {
        x[k] = min_temp + (max_temp - min_temp) * float(k) / float(num_temp);
        c->temperature = 273.15 + x[k];
        if(!k)
        {
            calculate_power(c, r);
            printf("%.3f\n", r->power/1e6);
        }
        else r->power *= pow((273.15 + x[k])/(x[k-1] + 273.15 ), 4);
        //calculate_power(c, support);
        y[k] = (r->power)/1e6;
        z[k] = 2 * SIGMA * pow(r->temperature, 4) * r->width * r->length / 1e6;
    }
    
    SetAtlasStyle();
    
     TCanvas *c1 = new TCanvas("c1","...",200,10,700,500);
     
   //calculate_temperature(r, 0.6);

   //c1->SetFillColor(42);
   c1->SetGrid();

   TGraph *gr = new TGraph(num_temp,x,y);
   gr->SetLineColor(2);
   gr->SetMarkerStyle(21);
   gr->SetTitle("Radiative power");
   gr->GetXaxis()->SetTitle("T (C)");
   gr->GetYaxis()->SetTitle("P (W)");
   gr->SetMinimum(0);
   gr->Draw("ACP");    
   
   TGraph *gr2 = new TGraph(num_temp,x,z);
   gr2->SetLineColor(kBlue);
   gr2->SetMarkerStyle(0);
   gr2->SetTitle("Radiative power");
   gr2->GetXaxis()->SetTitle("T (C)");
   gr2->GetYaxis()->SetTitle("P (W)");
   gr2->Draw("sameCP");
   
   c1->Update();    
}
