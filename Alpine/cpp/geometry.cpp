#include <TMath.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TH1.h>
#include <TH2.h>
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include <TVector.h>
#include <vector>
#include "AtlasStyle.C"

#define LAYERS 4

#define cos(x) (TMath::Cos((x)))
#define sin(x) (TMath::Sin((x)))
#define cosh(x) (TMath::CosH((x)))
#define sinh(x) (TMath::SinH((x)))
#define atanh(x) (TMath::ATanH((x)))
#define sqrt(x) (TMath::Sqrt((x)))
#define atan(x) (TMath::ATan((x)))
#define atan2(x,y) (TMath::ATan2((x),(y)))
#define exp(x) (TMath::Exp((x)))
#define pi TMath::Pi

TH2 *h_hits_r_z;
TH2 *h_hits_x_y;

void MatrixToEuler(const TMatrix &mat, TVector3 &vec)
{
    // x = theta
    // y = psi
    // z = phi
    if(abs(mat[2][0]) != 1)
    {
        vec.SetX(-asin(mat[2][0]));
        vec.SetY(atan2(mat[2][1]/cos(vec.X()), mat[2][2]/cos(vec.X())));
        vec.SetZ(atan2(mat[1][0]/cos(vec.X()), mat[0][0]/cos(vec.X())) );
    }
    else
    {
        vec.SetZ(0);
        if(abs(mat[2][0]) > 0)
        {
            vec.SetX(pi()/2.0);
            vec.SetY(atan2(mat[0][1], mat[0][2]));
        }
        else
        {
            vec.SetX(-pi()/2.0);
            vec.SetY(atan2(-mat[0][1], -mat[0][2]));
        }
    }
}

struct Module
{
    TVector3 pos;
    TVector3 normal; // normal to the surface
    
    // center pos
    float px, py, pz;
    
    // transformation matrix
    float a11, a12, a13, a21, a22, a23, a31, a32, a33;
    TMatrix transform, inverse;
    
    // surface cartesian equation coefficients (ax + by + cz + d = 0)
    float a, b, c, d;
    
    // dimensions 
    float width, length;
    
    // properties
    int type, layer, disk, p, e;
};

struct Ray
{
    TVector3 dir;
    TVector3 a, b;
    float eta, phi, z0;
    
    int hits[LAYERS], totalhits;
    
    vector<double *> hitvec;
};

vector<Module *> modules;
vector<Ray *> rays;

void load_modules()
{
FILE *fp = fopen("geometry_Alpine_XML.dat", "r");
    char line[128] = "";
    int mod_line = 0;
    Module *m = new Module();
    while(fgets(line, 128, fp) != NULL)
    {
        if(!mod_line && sscanf(line, "P %d %d %d %d %d", &m->type, &m->layer, &m->disk, &m->p, &m->e))
        {
            mod_line++;
            continue;
        }
        else
        if(mod_line == 1 && sscanf(line, "width x length %f %f", &m->width, &m->length))
        {
            mod_line++;
            continue;
        }
        else
        if(mod_line == 2 && sscanf(line, "pos : %f %f %f", &m->px, &m->py, &m->pz))
        {
            m->pos.SetX(m->px); m->pos.SetY(m->py); m->pos.SetZ(m->pz);
            mod_line++;
            continue;
        }
        else
        if(mod_line == 3 && sscanf(line, "%f %f %f", &m->a11, &m->a12, &m->a13))
        {
            mod_line++;
            continue;
        }
        else
        if(mod_line == 4 && sscanf(line, "%f %f %f", &m->a21, &m->a22, &m->a23))
        {
            mod_line++;
            continue;
        }
        else
        if(mod_line == 5 && sscanf(line, "%f %f %f", &m->a31, &m->a32, &m->a33))
        {
            m->transform.ResizeTo(3, 3);
            m->inverse.ResizeTo(3, 3);
            m->transform[0][0] = m->a11; m->transform[0][1] = m->a12; m->transform[0][2]  = m->a13;
            m->transform[1][0] = m->a21; m->transform[1][1] = m->a22; m->transform[1][2]  = m->a23;
            m->transform[2][0] = m->a31; m->transform[2][1] = m->a32; m->transform[2][2]  = m->a33;
            
            m->inverse = m->transform;
            m->inverse.Invert();
            
            TVector3 a, b;
            a[0] = 1; a[1] = 0; a[2] = 0;
            b[0] = 0; b[1] = 1; b[2] = 0;
            a = m->transform * a;
            b = m->transform * b;
            m->normal = a.Cross(b);
            m->normal *= 1/(m->normal.Mag());
            
            m->a = m->normal.X();
            m->b = m->normal.Y();
            m->c = m->normal.Z();
            m->d = -(m->a*m->px + m->b*m->py + m->c*m->pz);
           
            modules.push_back(m);
            mod_line = 0;
            printf("Module %.2f %.2f at %.3f %.3f %.3f (%.3f %.3f %.3f)\n", m->width, m->length, m->px, m->py, m->pz, m->normal.X(), m->normal.Y(), m->normal.Z());
            m = new Module();
        }
    }
    fclose(fp);
}

void create_rays(const int amount = 1000)
{
    TRandom *rand = new TRandom();
    for(int i = 0; i < amount; ++i)
    {
        Ray *r = new Ray();
        r->eta = rand->Uniform(-2.5, 2.5);
        r->phi = rand->Uniform(-pi(), +pi());
        r->z0 = rand->Uniform(-150, 150);
        float theta = 2*atan(exp(-r->eta));
        float d = 5000;
        r->a.SetX(0); r->a.SetY(0); r->a.SetZ(r->z0);
        r->b.SetX(d*sin(theta)*cos(r->phi)); r->b.SetY(d*sin(theta)*sin(r->phi)); r->b.SetZ(r->z0+d*cos(theta));
        r->dir.SetX(sin(theta)*cos(r->phi)); r->dir.SetY(sin(theta)*sin(r->phi)); r->dir.SetZ(cos(theta));
        rays.push_back(r);
        //printf("%.3f\n", 180*theta/pi());
        
        r->totalhits = 0;
        for(int j = 0; j < LAYERS; ++j) r->hits[j] = 0;
    }
}

void trace_rays()
{
    for(unsigned int i = 0; i < rays.size(); ++i)
    {
        Ray *r = rays[i];
        for(unsigned int j = 0; j < modules.size(); ++j)
        {
            // foreach ray and surface, find the intersection between the ray and the associated surface plan
            // then, check whether the intersection is inside the rectangle
            Module *m = modules[j];
            float t = -(m->d + m->c * r->z0) / (m->a * r->dir.X() + m->b * r->dir.Y() + m->c * r->dir.Z());
            if(t <= 0) continue;
            TVector3 inter, inter_atlas;
            inter.SetX(r->a.X() + r->dir.X() * t);
            inter.SetY(r->a.Y() + r->dir.Y() * t);
            inter.SetZ(r->a.Z() + r->dir.Z() * t);
            inter_atlas = inter;
            inter = inter - m->pos;
            inter = m->inverse * inter;
            
            if(abs(inter.X()) < m->width/2.0 && abs(inter.Y()) < m->length/2.0)
            {
                double *hr = new double[3];
                hr[0] = inter_atlas.X(); hr[1] = inter_atlas.Y(); hr[2] = inter_atlas.Z();
                //r->hitvec.push_back(hr);
                
                if(hr[2] >= 0) h_hits_r_z->Fill( hr[2], sqrt(hr[0]*hr[0]+hr[1]*hr[1]));
                
                h_hits_x_y->Fill( hr[0], hr[1]);
                
                r->hits[m->layer]++;
                r->totalhits++;
                //printf("%d %.3f %.3f %.3f\n", m->layer, inter.X(), inter.Y(), inter.Z());   
            }
        }
    }
}

void draw_ray(Ray *r)
{
    TPolyLine3D *pl = new TPolyLine3D(2);
    TVector3 dir(0, 0, 0), end(0, 0, 0);
    dir = r->b - r->a; dir *= 3000 / (dir.Mag());
    end = r->a + dir;
    pl->SetPoint(0, r->a.X(), r->a.Y(), r->a.Z());
    pl->SetPoint(1, end.X(), end.Y(), end.Z());
    pl->SetLineColor(kRed);
    pl->Draw();
    
    /*TPolyMarker3D *m = TPolyMarker3D(r->hitvec.size());
    m->SetMarkerColor(kGreen);
    for(unsigned int i = 0; i < r->hitvec.size(); ++i)
    {
        double *h = r->hitvec[i];
        m->SetPoint(i, h[0], h[1], h[2]);
    }*/
}

void render()
{
    //TRandom *rando = new TRandom();
    TGeoManager * geom = new TGeoManager("alpine", "Alpine Layout");
    TGeoMaterial *vacuumMatter = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *sensorMatter = new TGeoMaterial("Sensor"    , .938, 1., 10000.);
    TGeoMaterial *rayMatter = new TGeoMaterial("Ray", .938, 1., 10000);
    sensorMatter->SetTransparency(75);
    vacuumMatter->SetTransparency(100);
    
    TGeoMedium *vacuum = new TGeoMedium("vacuum", 1, vacuumMatter);
    TGeoMedium *sensor = new TGeoMedium("sensor", 1, sensorMatter);
    TGeoMedium *ray = new TGeoMedium("ray", 1, rayMatter);
    
    Double_t worldx = 600.;
    Double_t worldy = 600.;
    Double_t worldz = 2000.;
    
    TGeoVolume *top = geom->MakeBox("WORLD", vacuum, worldx, worldy, worldz);
    geom->SetTopVolume(top);
    
    TGeoVolume *normal_sensor  = geom->MakeBox("normal_sensor", sensor,  33.6, 19.8, 0.01);
    TGeoVolume *transition_sensor  = geom->MakeBox("transition_sensor", sensor,  16.8, 19.8, 0.01);
    TGeoVolume *ray_line = geom->MakeTube("ray_line", ray, 0, 0.01, 1);
    
    for(unsigned int i = 0; i < modules.size(); ++i)
    {
       Module *m = modules[i];
       TVector3 angles;
       MatrixToEuler(m->transform, angles);
       //angles.SetX(rando->Uniform(-2, 2));
       //angles.SetY(rando->Uniform(-2, 2));
       TGeoTranslation trans(m->px, m->py, m->pz);
       angles *= 180.0/pi();
       printf("%.3f %.3f %.3f\n", angles.X(), angles.Y(), angles.Z());
       TGeoRotation rot("euler", angles.Z(), angles.X(), angles.Y());
       
       if(m->width > 30) top->AddNode(normal_sensor, i, new TGeoCombiTrans(trans, rot));
       else top->AddNode(transition_sensor, i, new TGeoCombiTrans(trans, rot));
    }
    
   geom->CloseGeometry();
   geom->SetVisLevel(4);
   top->Draw("ogl");
   
   for(unsigned int i = 0; i < 1000; ++i) if(abs(rays[i]->eta) > 2.0)
   {
       draw_ray(rays[i]);
       break;
   }
   
    
}

void geometry()
{
    load_modules();
    printf("Loaded %d modules\n", (int)modules.size());
    create_rays(10000);
    
    gStyle->SetOptStat(0);
    SetAtlasStyle();
    gStyle->SetOptStat(0);
    
    gStyle->SetMarkerSize(0.1);
    
    h_hits_r_z = new TH2F("h_hits_r_z", "Hits (r,z)", 1000, 500, 0, 1200, 0, 300);
    h_hits_x_y = new TH2F("h_hits_x_y", "Hits (x,y)", 1000, -300, 300, 1000, -300, 300);
    
    
    trace_rays();
    
    TProfile *p_hits_eta = new TProfile("p_hits_eta", "p_hits_eta", 20, -3, 3, -0.5, 15.5);
    
    
    for(unsigned int i = 0; i < rays.size(); ++i)
    {
        Ray *r = rays[i];
        if(r->totalhits < 4)
        {
            printf("Total hits: %d (%d %d %d %d) (eta = %.2f, phi = %.2f, z0 = %.2f)\n", r->totalhits, r->hits[0], r->hits[1], r->hits[2], r->hits[3], r->eta, r->phi, r->z0);
        }
        p_hits_eta->Fill(r->eta, r->totalhits, 1);
    }
    
    //p_hits_eta->Draw();
    //h_hits_r_z->Draw();
    h_hits_x_y->Draw();
    //render();
    
}
