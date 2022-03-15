#include "Vec3.h"
#include <iomanip>
#include <cmath>

using namespace std;

Vec3::Vec3()
{
    this->x = 0.0;
    this->y = 0.0;
    this->z = 0.0;
    this->colorId = -1;
}

Vec3::Vec3(double x, double y, double z, int colorId)
{
    this->x = x;
    this->y = y;
    this->z = z;
    this->colorId = colorId;
}

Vec3::Vec3(const Vec3 &other)
{
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->colorId = other.colorId;
}

Vec3 Vec3::operator+(Vec3& rhs){

    return {x+rhs.x, y+rhs.y, z+rhs.z, -1};
}

void Vec3::operator+=(Vec3& rhs){
    x += rhs.x;
    y += rhs.y; 
    z += rhs.z;
}

Vec3 Vec3::operator*(double rhs){

    return {x*rhs, y*rhs, z*rhs, colorId};
}

void Vec3::operator*=(double rhs){
    x *= rhs;
    y *= rhs;
    z *= rhs;
}

Vec3 Vec3::operator-(Vec3& rhs){

    return {x-rhs.x, y-rhs.y, z-rhs.z, -1};
}

Vec3 Vec3::operator-(){

    return {-x, -y, -z, colorId};
}

void Vec3::operator-=(Vec3& rhs){
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
}



void Vec3::normalize(){

    double len = std::sqrt(x*x + y*y + z*z);
    x /= len;
    y /= len;
    z /= len;

}

ostream& operator<<(ostream& os, const Vec3& v) {
    
    os << fixed << setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << "]";

    return os;
}
