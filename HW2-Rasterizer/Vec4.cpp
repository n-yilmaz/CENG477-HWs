
#include "Vec4.h"
#include <iomanip>
#include "Vec3.h"

using namespace std;

Vec4::Vec4()
{
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->w = 0.0;
    this->colorId = -1;
}

Vec4::Vec4(double x, double y, double z, double w, int colorId)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
    this->colorId = colorId;
}

Vec4::Vec4(Vec3& v, int cId){
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
    this->w = 1;
    this->colorId = cId;
}

Vec4::Vec4(const Vec4 &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->w = other.w;
    this->colorId = other.colorId;
}


ostream& operator<<(ostream& os, const Vec4& v) {
    
    os << fixed << setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "]";

    return os;
}
