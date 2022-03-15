#ifndef __COLOR_H__
#define __COLOR_H__

#include <iostream>

class Color
{
public:
    double r, g, b;

    Color();
    Color(double r, double g, double b);
    Color(const Color &other);

    Color operator*(double val);
    Color operator+(Color& rhs);
    Color operator-(Color& rhs);
    Color operator/(double val);
    void operator+=(Color& rhs);

    Color operator+(Color&& rhs){
        return {r+rhs.r, g+rhs.g, b+rhs.b};
    }
    friend Color operator*(double val, Color& rhs){
    return rhs*val;
    }
    friend Color operator*(double val, Color&& rhs){
        return rhs*val;
    }
    friend std::ostream& operator<<(std::ostream& os, const Color& c);
};

#endif
