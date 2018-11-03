#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli Mod = 1e9 + 7;

lli modularInverse(lli a, lli m){
	lli r0 = a, r1 = m, ri, s0 = 1, s1 = 0, si;
	while(r1){
		si = s0 - s1 * (r0 / r1), s0 = s1, s1 = si;
		ri = r0 % r1, r0 = r1, r1 = ri;
	}
	if(s0 < 0) s0 += m;
	return s0;
}

lli powMod(lli a, int n, lli m){
	lli ans = 1;
	while(n){
		if(n & 1) ans = ans * a % m;
		n >>= 1;
		a = a * a % m;
	}
	return ans;
}

const lli inv_2 = modularInverse(2, Mod);
const lli inv_6 = modularInverse(6, Mod);
const lli inv_30 = modularInverse(30, Mod);

lli sum(lli n, int k){
	n %= Mod;
	if(k == 0) return n;
	if(k == 1) return n * (n + 1) % Mod * inv_2 % Mod;
	if(k == 2) return n * (n + 1) % Mod * (2*n + 1) % Mod * inv_6 % Mod;
	if(k == 3) return powMod(n * (n + 1) % Mod * inv_2 % Mod, 2, Mod);
	if(k == 4) return n * (n + 1) % Mod * (2*n + 1) % Mod * (3*n*(n+1)%Mod -1) % Mod * inv_30 % Mod;
	return 1;
}

//finds the sum of the kth powers of the primes
//less than or equal to n (0<=k<=4, add more if you need)
lli SumPrimePi(lli n, int k){
	lli v = sqrt(n), p, temp, q, j, end, i, d;
	vector<lli> lo(v+2), hi(v+2);
	vector<bool> used(v+2);
	for(p = 1; p <= v; p++){
		lo[p] = sum(p, k) - 1;
		hi[p] = sum(n/p, k) - 1;
	}
	for(p = 2; p <= v; p++){
		if(lo[p] == lo[p-1]) continue;
		temp = lo[p-1];
		q = p * p;
		hi[1] -= (hi[p] - temp) * powMod(p, k, Mod) % Mod;
		if(hi[1] < 0) hi[1] += Mod;
		j = 1 + (p & 1);
		end = (v <= n/q) ? v : n/q;
		for(i = p + j; i <= 1 + end; i += j){
			if(used[i]) continue;
			d = i * p;
			if(d <= v)
				hi[i] -= (hi[d] - temp) * powMod(p, k, Mod) % Mod;
			else
				hi[i] -= (lo[n/d] - temp) * powMod(p, k, Mod) % Mod;
			if(hi[i] < 0) hi[i] += Mod;
		}
		if(q <= v)
			for(i = q; i <= end; i += p*j)
				used[i] = true;
		for(i = v; i >= q; i--){
			lo[i] -= (lo[i/p] - temp) * powMod(p, k, Mod) % Mod;
			if(lo[i] < 0) lo[i] += Mod;
		}
	}
	return hi[1] % Mod;
}

int main(){
	lli n;
	int k;
	cin >> n >> k;
	clock_t start = clock();
	lli ans = SumPrimePi(n, k);
	clock_t end = clock();
	cout << ans << "\n" << (double)(end - start) / (double)CLOCKS_PER_SEC << "s\n";
	return 0;
}