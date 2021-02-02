#include <bits/stdc++.h>
using namespace std;
using lli = long long int;

const lli Mod = 1e9 + 7;

template<typename T>
struct SumPrimePi{
	int v, k;
	lli n;
	vector<T> lo, hi;
	vector<int> primes;

	SumPrimePi(lli n, int k = 0): n(n), v(sqrt(n)), k(k){
		lo.resize(v+2), hi.resize(v+2);
	}

	T power(T a, lli b){
		T ans = 1;
		while(b){
			if(b & 1) ans *= a;
			b >>= 1;
			a *= a;
		}
		return ans;
	}

	T powerSum(T n, int k){
		if(k == 0) return n;
		if(k == 1) return n * (n + 1) / 2;
		return 0;
	}

	void build(){
		lli p, q, j, end, i, d;
		T temp;
		for(p = 1; p <= v; p++){
			lo[p] = powerSum(p, k) - 1;
			hi[p] = powerSum(n/p, k) - 1;
		}
		for(p = 2; p <= v; p++){
			T pk = power(p, k);
			if(lo[p] == lo[p-1]) continue;
			primes.push_back(p);
			temp = lo[p-1];
			q = p * p;
			end = (v <= n/q) ? v : n/q;
			for(i = 1; i <= end; ++i){
				d = i * p;
				if(d <= v)
					hi[i] -= (hi[d] - temp) * pk;
				else
					hi[i] -= (lo[n/d] - temp) * pk;
			}
			for(i = v; i >= q; i--){
				lo[i] -= (lo[i/p] - temp) * pk;
			}
		}
	}

	T get(lli i) const{
		if(i <= v) return lo[i];
		else return hi[n/i];
	}
};

template<typename T>
struct MultiplicativeSum{
	int v;
	lli n;
	vector<T> lo, hi, smallFP;
	vector<int> primes;

	MultiplicativeSum(lli n, const vector<int> & primes): n(n), v(sqrt(n)), primes(primes){
		lo.resize(v+2), hi.resize(v+2), smallFP.resize(v+2);
	}

	void add(T coef, const auto & pi){
		assert(pi.n == n);
		for(int i = 1; i <= v; ++i){
			smallFP[i] += coef * pi.get(i);
			hi[i] += coef * (pi.get(n/i) - pi.get(v));
		}
	}

	T getAdded(lli i, lli p){
		if(i <= v){
			return lo[i] + smallFP[max(i, p)] - smallFP[p];
		}else{
			return hi[n/i] + smallFP[v] - smallFP[p];
		}
	}

	void build(function<T(lli, int)> g){
		for(int i = 1; i <= v; ++i){
			lo[i] += 1;
			hi[i] += 1;
		}
		for(int r = (int)primes.size()-1; r >= 0; --r){
			lli p = primes[r];
			vector<lli> p_power(1, 1);
			vector<T> gs(1, T(1));
			lli p_pow = p;
			for(int e = 1; ; ++e){
				p_power.push_back(p_pow);
				gs.push_back(g(p, e));
				if(p_pow > n/p) break;
				p_pow *= p;
			}
			for(int i = 1; i <= v; ++i){
				lli next = n / i;
				if(next < p*p) break;
				for(int e = 1; e < p_power.size() && p_power[e] <= next; ++e){
					hi[i] += gs[e] * getAdded(next / p_power[e], p);
				}
				hi[i] -= gs[1];
			}
			for(int i = v; i >= 1; --i){
				if(i < p*p) break;
				for(int e = 1; e <= p_power.size() && p_power[e] <= i; ++e){
					lo[i] += gs[e] * getAdded(i / p_power[e], p);
				}
				lo[i] -= gs[1];
			}
		}
		for(int i = 1; i <= v; ++i){
			lo[i] += smallFP[i];
			hi[i] += smallFP[v];
		}
	}

	T get(lli i) const{
		if(i <= v) return lo[i];
		else return hi[n/i];
	}
};

template<typename T>
T f(const SumPrimePi<T> & pi, function<T(lli, int)> g, lli n, int idx = 0){
	// sum of g(p, 1) for primes p such that primes[idx] <= p <= n
	int lo = idx ? pi.primes[idx-1] : 0;
	T ans = pi.get(n) - pi.get(lo);
	if(idx == 0) ans++;
	for(int i = idx; i < pi.primes.size(); ++i){
		lli p = pi.primes[i];
		if(p * p > n) break;
		int e = 1;
		lli curr = n / p;
		while(curr >= p){
			ans += g(p, e) * f(pi, g, curr, i+1) + g(p, e+1);
			curr /= p;
			++e;
		}
	}
	return ans;
}

int main(){
	lli n;
	int k;
	cin >> n >> k;
	clock_t start = clock();
	SumPrimePi<lli> pi(n, k);
	pi.build();
	lli ans = pi.get(n);
	clock_t end = clock();
	cout << ans << "\n" << (double)(end - start) / (double)CLOCKS_PER_SEC << "s\n";
	return 0;
}