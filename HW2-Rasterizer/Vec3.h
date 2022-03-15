#ifndef __VEC3_H__
#define __VEC3_H__

#include <iostream>
using namespace std;

class Vec3
{
public:
    double x, y, z;
    int colorId;

    Vec3();
    Vec3(double x, double y, double z, int colorId);
    Vec3(const Vec3 &other);

    void normalize();


    Vec3 operator*(double scalar);
    Vec3 operator/(double scalar){
        return {x/scalar, y/scalar, z/scalar, colorId};
    }
    friend Vec3 operator*(double scalar, Vec3& rhs){
        return rhs*scalar;
    }

    friend Vec3 operator*(double scalar, Vec3&& rhs){
        return rhs*scalar;
    }
    Vec3 operator+(Vec3&& rhs){
        return {x+rhs.y, y+rhs.y, z+rhs.z, -1};
    }
    Vec3 operator+(Vec3& rhs);
    void operator+=(Vec3& rhs);
    void operator+=(Vec3&& rhs){
        x += rhs.x;
        y += rhs.y; 
        z += rhs.z;
    }
    Vec3 operator-();
    Vec3 operator-(Vec3& rhs);
    Vec3 operator-(Vec3&& rhs){
        return {x-rhs.x, y-rhs.y, z-rhs.z, -1};
    }
    void operator-=(Vec3& rhs);
    void operator-=(Vec3&& rhs){
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
    }
    void operator*=(double scalar);

    
    friend std::ostream& operator<<(std::ostream& os, const Vec3& v);
};

#endif
