#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli mod = 1e7 + 19;

//Solves a linear recurrence relation of degree d of the form
//F(n) = a(d-1)*F(n-1) + a(d-2)*F(n-2) + ... + a(1)*F(n-(d-1)) + a(0)*F(n-d)
//with initial values F(0), F(1), ..., F(d-1)
//It finds the nth term of the recurrence, F(n)
//The values of a[0,...,d) are in the array P[]
lli solveRecurrence(lli *P, lli *init, int deg, lli n){
	lli *ans = new lli[deg]();
	lli *R = new lli[2*deg]();
	ans[0] = 1;
	lli p = 1;
	for(lli v = n; v >>= 1; p <<= 1);
	auto mult = [&](int d){
		fill(R, R + 2*deg, 0);
		for(int i = 0; i < deg; i++)
			for(int j = 0; j < deg; j++)
				R[i + j + d] += ans[i] * ans[j];

		for(int i = 0; i < 2*deg; ++i) R[i] %= mod;
		for(int i = deg-1; i >= 0; i--){
			R[i + deg] %= mod;
			for(int j = 0; j < deg; j++)
				R[i + j] += R[i + deg] * P[j];
		}

		for(int i = 0; i < deg; i++) R[i] %= mod;
		copy(R, R + deg, ans);
	};
	while(p){
		mult((n & p)!=0);
		p >>= 1;
	}
	lli nValue = 0;
	for(int i = 0; i < deg; i++)
		nValue += ans[i] * init[i];
	return nValue % mod;
}

int main(){
	int deg;
	cin >> deg;
	lli *P = new lli[deg];
	lli *init = new lli[deg];
	for(int i = 0; i < deg; i++)
		cin >> P[i];
	for(int i = 0; i < deg; i++)
		cin >> init[i];
	lli n;
	cin >> n;
	lli F_n = solveRecurrence(P, init, deg, n);
	cout << F_n;
	return 0;
}