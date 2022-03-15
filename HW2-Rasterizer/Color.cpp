#include "Color.h"
#include <iostream>
#include <iomanip>

using namespace std;

Color::Color() {}

Color::Color(double r, double g, double b)
{
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

Color Color::operator*(double val){
    return {r*val, g*val, b*val};
}

Color Color::operator+(Color& rhs){
    return {r+rhs.r, g+rhs.g, b+rhs.b};
}

Color Color::operator-(Color& rhs){
    return {r-rhs.r, g-rhs.g, b-rhs.b};
}

Color Color::operator/(double rhs){
    return {r/rhs, g/rhs, b/rhs};
}

void Color::operator+=(Color& rhs){
    this->r += rhs.r;
    this->g += rhs.g;
    this->b += rhs.b;
}

ostream& operator<<(ostream& os, const Color& c)
{
    os << fixed << setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}
