#include <bits/stdc++.h>
using namespace std;
typedef long long int ull;

ull potencia(ull n, ull p){
	ull div = p, pot = 0;
	while(div <= n){
		pot += n / div;
		div *= p;
	}
	return pot;
}

void criba(ull n, vector<ull> & primos){
	vector<bool> es_primo(n + 1, true);
	primos.push_back(2);
	for(ull i = 3; i <= n; i += 2){
		if(es_primo[i]){
			primos.push_back(i);
			for(ull j = i * i; j <= n; j += 2 * i){
				es_primo[j] = false;
			}
		}
	}
}

void factorizar_factorial(ull n, map<ull, ull> & f, vector<ull> & primos){
	for(ull & p : primos){
		if(p > n) break;
		ull pot = potencia(n, p);
		if(pot > 0){
			f[p] = pot;
		}
	}
}

void factorizar_dividir(ull n, map<ull, ull> & f, vector<ull> & primos){
	for(ull & p : primos){
		if(p * p > n) break;
		if(n % p == 0){
			ull pot = 0;
			while(n % p == 0){
				n /= p;
				++pot;
			}
			f[p] -= pot;
		}
	}
	if(n > 1){
		f[n] -= 1;
	}
}

int main(){
	ios_base::sync_with_stdio(0);
	vector<ull> primos;
	criba(100, primos);
	ull n, d;
	while(cin >> n >> d && !(n == 0 && d == 0)){
		d = abs(d);
		map<ull, ull> f;
		factorizar_factorial(n, f, primos);
		factorizar_dividir(d, f, primos);
		ull ans = 1;
		for(pair<const ull, ull> & par : f){
			if(par.second < 0){
				ans = 0;
				break;
			}
			ans *= par.second + 1;
		}
		cout << ans << "\n";
	}
	return 0;
}