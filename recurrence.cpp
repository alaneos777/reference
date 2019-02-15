#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli mod = 1e7 + 19;

//Solves a linear homogeneous recurrence relation of degree "deg" of the form
//F(n) = a(d-1)*F(n-1) + a(d-2)*F(n-2) + ... + a(1)*F(n-(d-1)) + a(0)*F(n-d)
//with initial values F(0), F(1), ..., F(d-1)
//It finds the nth term of the recurrence, F(n)
//The values of a[0,...,d) are in the array P[]
lli solveRecurrence(const vector<lli> & P, const vector<lli> & init, lli n){
	int deg = P.size();
	vector<lli> ans(deg), R(2*deg);
	ans[0] = 1;
	lli p = 1;
	for(lli v = n; v >>= 1; p <<= 1);
	do{
		int d = (n & p) != 0;
		fill(R.begin(), R.end(), 0);
		//only if deg(mod-1)^2 overflows, do mod in all the multiplications
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
		copy(R.begin(), R.begin() + deg, ans.begin());
	}while(p >>= 1);
	lli nValue = 0;
	for(int i = 0; i < deg; i++)
		nValue += ans[i] * init[i];
	return nValue % mod;
}

int main(){
	int deg;
	cin >> deg;
	vector<lli> P(deg), init(deg);
	for(int i = 0; i < deg; i++)
		cin >> P[i];
	for(int i = 0; i < deg; i++)
		cin >> init[i];
	lli n;
	cin >> n;
	lli F_n = solveRecurrence(P, init, n);
	cout << F_n;
	return 0;
}