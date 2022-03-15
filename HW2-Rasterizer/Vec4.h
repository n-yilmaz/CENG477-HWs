#ifndef __VEC4_H__
#define __VEC4_H__

#include <iostream>
#include "Vec3.h"

using namespace std;

class Vec4
{
public:
    double x, y, z, w;
    int colorId;


    Vec4();
    Vec4(double x, double y, double z, double w, int colorId);
    Vec4(const Vec4 &other);
    Vec4(Vec3& v, int cId);

    Vec4 operator*(double t){
    return {x*t, y*t, z*t, w*t, colorId};
    }
  
    Vec4 operator-(Vec4&& rhs){
    return {x-rhs.x, y-rhs.y, z-rhs.z, w-rhs.w, -1};
    }

    Vec4 operator-(Vec4& rhs){
    return {x-rhs.x, y-rhs.y, z-rhs.z, w-rhs.w, -1};
    }

    Vec4 operator+(Vec4& rhs){
    return {x+rhs.x, y+rhs.y, z+rhs.z, w+rhs.w, -1};
    }
    
    Vec4 operator+(Vec4&& rhs){
    return {x+rhs.x, y+rhs.y, z+rhs.z, w+rhs.w, -1};
    }

    Vec4 operator-(){
    return {-x, -y, -z, -w, colorId};
    }


    friend std::ostream& operator<<(std::ostream& os, const Vec4& v);
};

#endif
