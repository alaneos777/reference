#include <bits/stdc++.h>

using namespace std;

typedef long long int lli;

struct fraccion{
    lli num, den;
    fraccion(){
        num = 0, den = 1;
    }
    fraccion(lli x, lli y){
        if(y < 0)
            x *= -1, y *=-1;
        lli d = __gcd(abs(x), abs(y));
        num = x/d, den = y/d;
    }
    fraccion(lli v){
        num = v;
        den = 1;
    }
    fraccion operator+(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return fraccion(num*(f.den/d) + f.num*(den/d), den*(f.den/d));
    }
    fraccion operator-() const{
        return fraccion(-num, den);
    }
    fraccion operator-(const fraccion& f) const{
        return *this + (-f);
    }
    fraccion operator*(const fraccion& f) const{
        return fraccion(num*f.num, den*f.den);
    }
    fraccion operator/(const fraccion& f) const{
        return fraccion(num*f.den, den*f.num);
    }
    fraccion operator+=(const fraccion& f){
        *this = *this + f;
        return *this;
    }
    fraccion operator-=(const fraccion& f){
        *this = *this - f;
        return *this;
    }
    fraccion operator++(int xd){
        *this = *this + 1;
        return *this;
    }
    fraccion operator--(int xd){
        *this = *this - 1;
        return *this;
    }
    fraccion operator*=(const fraccion& f){
        *this = *this * f;
        return *this;
    }
    fraccion operator/=(const fraccion& f){
        *this = *this / f;
        return *this;
    }
    bool operator==(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return (num*(f.den/d) == (den/d)*f.num);
    }
    bool operator!=(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return (num*(f.den/d) != (den/d)*f.num);
    }
    bool operator >(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return (num*(f.den/d) > (den/d)*f.num);
    }
    bool operator <(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return (num*(f.den/d) < (den/d)*f.num);
    }
    bool operator >=(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return (num*(f.den/d) >= (den/d)*f.num);
    }
    bool operator <=(const fraccion& f) const{
        lli d = __gcd(den, f.den);
        return (num*(f.den/d) <= (den/d)*f.num);
    }
    fraccion inverso() const{
        return fraccion(den, num);
    }
    fraccion fabs() const{
        fraccion nueva;
        nueva.num = abs(num);
        nueva.den = den;
        return nueva;
    }
    double value() const{
    	return (double)num / (double)den;
    }
    string str() const{
        stringstream ss;
        ss << num;
        if(den != 1) ss << "/" << den;
        return ss.str();
    }
};

ostream &operator<<(ostream &os, const fraccion & f) { 
    return os << f.str();
}

istream &operator>>(istream &is, fraccion & f){
    lli num = 0, den = 1;
    string str;
    is >> str;
    size_t pos = str.find("/");
    if(pos == string::npos){
        istringstream(str) >> num;
    }else{
        istringstream(str.substr(0, pos)) >> num;
        istringstream(str.substr(pos + 1)) >> den;
    }
    f = fraccion(num, den);
    return is;
}