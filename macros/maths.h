#define cos(x) (TMath::Cos((x)))
#define sin(x) (TMath::Sin((x)))
#define cosh(x) (TMath::CosH((x)))
#define sinh(x) (TMath::SinH((x)))
#define atanh(x) (TMath::ATanH((x)))
#define sqrt(x) (TMath::Sqrt((x)))
#define atan(x) (TMath::ATan((x)))
#define exp(x) (TMath::Exp((x)))
#define pi TMath::Pi

struct vec
{
    float x, y, z;
    
    vec() { x = y = z = 0; }
    
    vec(float xn, float yn, float zn)
    {
        x = xn; y = yn; z = zn;
    }
    
    vec &add(const vec v)
    {
        this->x += v.x;
        this->y += v.y;
        this->z += v.z;
        return *this;
    }
    
    vec &sub(const vec v)
    {
        this->x -= v.x;
        this->y -= v.y;
        this->z -= v.z;
        return *this;
    }
    
    const float norm() { return sqrt(x*x + y*y + z*z); }
    
    vec &pr2c(float pt, float phi, float eta)
    {
        x = pt * cos(phi);
        y = pt * sin(phi);
        z = pt * sinh(eta);
        return *this;
    }
    
    vec &c2pr(float xn, float yn, float zn)
    {
        float p = sqrt(xn*xn + yn*yn + zn*zn);
        x = zn;
        z = atanh(zn/p);
        return *this;
    }
};
