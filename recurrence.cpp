#include <bits/stdc++.h>
using namespace std;
using lli = long long int;

const lli mod = 1e7 + 19;

//Solves a linear homogeneous recurrence relation of degree "deg"
//of the form F(n) = a(d-1)*F(n-1) + a(d-2)*F(n-2) + ... + a(1)*F(n-(d-1)) + a(0)*F(n-d)
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
		for(int i = 0; i < deg; i++)
			for(int j = 0; j < deg; j++)
				(R[i + j + d] += ans[i] * ans[j]) %= mod;
		for(int i = deg-1; i >= 0; i--)
			for(int j = 0; j < deg; j++)
				(R[i + j] += R[i + deg] * P[j]) %= mod;
		copy(R.begin(), R.begin() + deg, ans.begin());
	}while(p >>= 1);
	lli nValue = 0;
	for(int i = 0; i < deg; i++)
		(nValue += ans[i] * init[i]) %= mod;
	return nValue;
}

lli inverse(lli a, lli n){
	lli r0 = a, r1 = n, ri, s0 = 1, s1 = 0, si;
	while(r1){
		si = s0 - s1 * (r0 / r1), s0 = s1, s1 = si;
		ri = r0 % r1, r0 = r1, r1 = ri;
	}
	assert(r0 == 1);
	if(s0 < 0) s0 += n;
	return s0;
}

//Finds the shortest linear recurrence relation for the
//given init values. Only works for prime modulo.
vector<lli> BerlekampMassey(const vector<lli> & init){
	vector<lli> cur, ls;
	lli ld;
	for(int i = 0, m; i < init.size(); ++i){
		lli eval = 0;
		for(int j = 0; j < cur.size(); ++j)
			eval = (eval + init[i-j-1] * cur[j]) % mod;
		eval -= init[i];
		if(eval < 0) eval += mod;
		if(eval == 0) continue;
		if(cur.empty()){
			cur.resize(i + 1);
			m = i;
			ld = eval;
		}else{
			lli k = eval * inverse(ld, mod) % mod;
			vector<lli> c(i - m - 1);
			c.push_back(k);
			for(int j = 0; j < ls.size(); ++j)
				c.push_back((mod-ls[j]) * k % mod);
			if(c.size() < cur.size()) c.resize(cur.size());
			for(int j = 0; j < cur.size(); ++j){
				c[j] += cur[j];
				if(c[j] >= mod) c[j] -= mod;
			}
			if(i - m + ls.size() >= cur.size())
				ls = cur, m = i, ld = eval;
			cur = c;
		}
	}
	if(cur.empty()) cur.push_back(0);
	reverse(cur.begin(), cur.end());
	return cur;
}

int main(){
	/*int deg;
	cin >> deg;
	vector<lli> P(deg), init(deg);
	for(int i = 0; i < deg; i++)
		cin >> P[i];
	for(int i = 0; i < deg; i++)
		cin >> init[i];
	lli n;
	cin >> n;
	lli F_n = solveRecurrence(P, init, n);
	cout << F_n;*/
	vector<lli> init;
	lli vi, n;
	string str;
	getline(cin, str);
	stringstream ss(str);
	while(ss >> vi) init.push_back(vi);
	auto P = BerlekampMassey(init);
	init.resize(P.size());
	for(int pi : P) cout << pi << " "; cout << "\n";
	getline(cin, str);
	ss = stringstream(str);
	while(ss >> n) cout << solveRecurrence(P, init, n) << " ";
	return 0;
}