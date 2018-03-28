#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli mod = 1e7 + 19;

void multByOne(vector<lli> & polynomial, vector<lli> & original){
	lli first = polynomial.back();
	for(int i = polynomial.size() - 1; i >= 0; i--){
		polynomial[i] = (first * original[i]) % mod;
		if(i > 0){
			polynomial[i] += polynomial[i - 1];
			polynomial[i] %= mod;
		}
	}
}

vector<lli> mult(vector<lli> & P, vector<lli> & Q, vector<vector<lli>> & residues){
	int degree = P.size();
	vector<lli> R(degree), S(degree - 1);
	for(int i = 0; i < degree; i++){
		for(int j = 0; j < degree; j++){
			int pos = i + j;
			if(pos < degree){
				R[pos] += (P[i] * Q[j]) % mod;
				R[pos] %= mod;
			}else{
				pos -= degree;
				S[pos] += (P[i] * Q[j]) % mod;
				S[pos] %= mod;
			}
		}
	}
	for(int i = 0; i < degree - 1; i++){
		for(int j = 0; j < degree; j++){
			R[j] += (S[i] * residues[i][j]) % mod;
			R[j] %= mod;
		}
	}
	return R;
}

lli solveRecurrence(vector<lli> & charPoly, vector<lli> & initValues, vector<vector<lli>> & residues, lli n){
	int degree = charPoly.size();
	vector<lli> tmp(degree), ans(degree);
	ans[0] = 1;
	if(degree > 1){
		tmp[1] = 1;
	}else{
		tmp[0] = charPoly[0];
	}
	while(n){
		if(n & 1) ans = mult(ans, tmp, residues);
		n >>= 1;
		if(n) tmp = mult(tmp, tmp, residues);
	}
	lli nValue = 0;
	for(int i = 0; i < degree; i++){
		nValue += (ans[i] * initValues[i]) % mod;
		nValue %= mod;
	}
	if(nValue < 0) nValue += mod;
	return nValue;
}

int main(){
	int degree;
	cin >> degree;
	vector<lli> charPoly(degree), initValues(degree);
	for(int i = 0; i < degree; i++){
		cin >> charPoly[i];
	}
	for(int i = 0; i < degree; i++){
		cin >> initValues[i];
	}
	vector<vector<lli>> residues;
	vector<lli> current = charPoly;
	for(int i = 0; i < degree - 1; i++){
		residues.push_back(current);
		if(i != degree - 2) multByOne(current, charPoly);
	}

	lli n;
	cin >> n;
	lli F_n = solveRecurrence(charPoly, initValues, residues, n);
	cout << F_n;
	return 0;
}