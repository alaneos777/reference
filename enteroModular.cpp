#include <bits/stdc++.h>

using namespace std;

typedef long long int lli;

lli mod(lli a, lli b){
    lli ans = a % b;
    if(ans < 0) ans += b;
    return ans;
}

struct enteroModular{
    lli a, n;
    enteroModular(lli x, lli y){
    	if(y == 0) a = x;
        else a = mod(x, y);
        n = y;
    }
    enteroModular(lli x){
        a = x, n = 0;
    }
    enteroModular(){
        a = 0, n = 0;
    }
    enteroModular operator+(const enteroModular & e) const{
        return enteroModular(a + e.a, max(n, e.n));
    }
    enteroModular operator-() const{
        return enteroModular(-a, n);
    }
    enteroModular operator-(const enteroModular & e) const{
        return *this + (-e);
    }
    enteroModular operator*(const enteroModular & e) const{
        return enteroModular(a * e.a, max(n, e.n));
    }
    enteroModular inverso() const{
        lli r0 = a, r1 = n, ri, s0 = 1, s1 = 0, si;
        while(r1 != 0){
            ri = r0 % r1;
            si = s0 - s1 * (r0 / r1);
            r0 = r1, r1 = ri;
            s0 = s1, s1 = si;
        }
        return enteroModular(s0, n);
    }
    enteroModular operator/(const enteroModular & e) const{
    	enteroModular tmp(e.a, max(n, e.n));
        return enteroModular(a * tmp.inverso().a, max(n, e.n));
    }
    enteroModular operator+=(const enteroModular & e){
        *this = *this + e;
        return *this;
    }
    enteroModular operator-=(const enteroModular & e){
        *this = *this - e;
        return *this;
    }
    enteroModular operator*=(const enteroModular & e){
        *this = *this * e;
        return *this;
    }
    enteroModular operator/=(const enteroModular & e){
        *this = *this / e;
        return *this;
    }
    enteroModular operator++(int xd){
        *this = *this + enteroModular(1, n);
        return * this;
    }
    enteroModular operator--(int xd){
        *this = *this - enteroModular(1, n);
        return * this;
    }
    bool operator==(const enteroModular & e) const{
    	if(max(n, e.n) == 0)
    		return a == e.a;
        return mod(a, max(n, e.n)) == mod(e.a, max(n, e.n));
    }
    bool operator!=(const enteroModular & e) const{
    	if(max(n, e.n) == 0)
    		return a != e.a;
        return mod(a, max(n, e.n)) != mod(e.a, max(n, e.n));
    }
    string str() const{
    	stringstream ss;
    	ss << a << "," << n;
    	return ss.str();
    }
};

ostream &operator<<(ostream &os, const enteroModular & e) { 
    return os << e.a;
}