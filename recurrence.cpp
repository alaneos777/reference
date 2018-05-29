#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli mod = 1e7 + 19;

void multByOne(lli *polynomial, lli *original, int degree){
	lli first = polynomial[degree - 1];
	for(int i = degree - 1; i >= 0; --i){
		polynomial[i] = first * original[i];
		if(i > 0) polynomial[i] += polynomial[i - 1];
	}
	for(int i = 0; i < degree; ++i) polynomial[i] %= mod;
}

lli *mult(lli *P, lli *Q, lli **residues, int degree){
	lli *R = new lli[degree]();
	lli *S = new lli[degree - 1]();
	for(int i = 0; i < degree; i++){
		for(int j = 0; j < degree; j++){
			if(i + j < degree) R[i + j] += P[i] * Q[j];
			else S[i + j - degree] += P[i] * Q[j];
		}
	}
	for(int i = 0; i < degree - 1; i++) S[i] %= mod;
	for(int i = 0; i < degree - 1; i++){
		for(int j = 0; j < degree; j++)
			R[j] += S[i] * residues[i][j];
	}
	for(int i = 0; i < degree; i++) R[i] %= mod;
	return R;
}

lli **getResidues(lli *charPoly, int degree){
	lli **residues = new lli*[degree - 1];
	lli *current = new lli[degree];
	copy(charPoly, charPoly + degree, current);
	for(int i = 0; i < degree - 1; i++){
		residues[i] = new lli[degree];
		copy(current, current + degree, residues[i]);
		if(i != degree - 2) multByOne(current, charPoly, degree);
	}
	return residues;
}

//Solves a linear recurrence relation of degree d of the form
//F(n) = a(d-1)*F(n-1) + a(d-2)*F(n-2) + ... + a(1)*F(n-(d-1)) + a(0)*F(n-d)
//with initial values F(0), F(1), ..., F(d-1)
//It finds the nth term of the recurrence, F(n)
//The values of a[0,...,d) are in the array charPoly[]
lli solveRecurrence(lli *charPoly, lli *initValues, int degree, lli n){
	lli **residues = getResidues(charPoly, degree);
	lli *tmp = new lli[degree]();
	lli *ans = new lli[degree]();
	ans[0] = 1;
	if(degree > 1) tmp[1] = 1;
	else tmp[0] = charPoly[0];
	while(n){
		if(n & 1) ans = mult(ans, tmp, residues, degree);
		n >>= 1;
		if(n) tmp = mult(tmp, tmp, residues, degree);
	}
	lli nValue = 0;
	for(int i = 0; i < degree; i++) nValue += ans[i] * initValues[i];
	return nValue % mod;
}

int main(){
	int degree;
	cin >> degree;
	lli *charPoly = new lli[degree];
	lli *initValues = new lli[degree];
	for(int i = 0; i < degree; i++)
		cin >> charPoly[i];
	for(int i = 0; i < degree; i++)
		cin >> initValues[i];
	lli n;
	cin >> n;
	lli F_n = solveRecurrence(charPoly, initValues, degree, n);
	cout << F_n;
	return 0;
}