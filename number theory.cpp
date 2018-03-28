#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

lli piso(lli a, lli b){
	if((a >= 0 && b > 0) || (a < 0 && b < 0)){
		return a / b;
	}else{
		if(a % b == 0) return a / b;
		else return a / b - 1;
	}
}

lli techo(lli a, lli b){
	if((a >= 0 && b > 0) || (a < 0 && b < 0)){
		if(a % b == 0) return a / b;
		else return a / b + 1;
	}else{
		return a / b;
	}
}

lli pow(lli b, lli e){
	lli ans = 1;
	while(e){
		if(e & 1) ans *= b;
		e >>= 1;
		b *= b;
	}
	return ans;
}

lli powMod(lli b, lli e, lli m){
	lli ans = 1;
	while(e){
		if(e & 1) ans = (ans * b) % m;
		e >>= 1;
		b = (b * b) % m;
	}
	return ans;
}

lli multMod(lli a, lli b, lli n){
	lli ans = 0;
	a %= n, b %= n;
	if(abs(b) > abs(a)) swap(a, b);
	if(b < 0){
		a *= -1, b *= -1;
	}
	while(b){
		if(b & 1) ans = (ans + a) % n;
		b >>= 1;
		a = (a + a) % n;
	}
	return ans;
}

lli gcd(lli a, lli b){
	lli r;
	while(b != 0) r = a % b, a = b, b = r;
	return a;
}

lli lcm(lli a, lli b){
	return b * (a / gcd(a, b));
}

lli gcd(vector<lli> & nums){
	lli ans = 0;
	for(lli & num : nums) ans = gcd(ans, num);
	return ans;
}

lli lcm(vector<lli> & nums){
	lli ans = 1;
	for(lli & num : nums) ans = lcm(ans, num);
	return ans;
}

lli extendedGcd(lli a, lli b, lli & s, lli & t){
	lli q, r0 = a, r1 = b, ri, s0 = 1, s1 = 0, si, t0 = 0, t1 = 1, ti;
	while(r1){
		q = r0 / r1;
		ri = r0 % r1, r0 = r1, r1 = ri;
		si = s0 - s1 * q, s0 = s1, s1 = si;
		ti = t0 - t1 * q, t0 = t1, t1 = ti;
	} 
	s = s0, t = t0;
	return r0;
}

lli modularInverse(lli a, lli m){
	lli r0 = a, r1 = m, ri, s0 = 1, s1 = 0, si;
	while(r1){
		si = s0 - s1 * (r0 / r1), s0 = s1, s1 = si;
		ri = r0 % r1, r0 = r1, r1 = ri;
	} 
	if(s0 < 0) s0 += m;
	return s0;
}

pair<lli, lli> chinese(vector<lli> & a, vector<lli> & n){
	lli prod = 1, p, ans = 0;
	for(lli & ni : n) prod *= ni;
	for(int i = 0; i < a.size(); i++){
		p = prod / n[i];
		ans = (ans + (a[i] % n[i]) * modularInverse(p, n[i]) % prod * p) % prod;
	}
	if(ans < 0) ans += prod;
	return make_pair(ans, prod);
}

vector<lli> divisorsSum;
vector<vector<lli>> divisors;
void divisorsSieve(lli n){
	divisorsSum.resize(n + 1, 0);
	divisors.resize(n + 1, vector<lli>());
	for(lli i = 1; i <= n; i++){
		for(lli j = i; j <= n; j += i){
			divisorsSum[j] += i;
			divisors[j].push_back(i);
		}
	}
}

vector<lli> primes;
vector<bool> isPrime;
void primesSieve(lli n){
	isPrime.resize(n + 1, true);
	isPrime[0] = isPrime[1] = false;
	primes.push_back(2);
	for(lli i = 4; i <= n; i += 2){
		isPrime[i] = false;
	}
	for(lli i = 3; i <= n; i += 2){
		if(isPrime[i]){
			primes.push_back(i);
			for(lli j = i * i; j <= n; j += 2 * i){
				isPrime[j] = false;
			}
		}
	}
}

vector<lli> lowestPrime;
void lowestPrimeSieve(lli n){
	lowestPrime.resize(n + 1, 1);
	lowestPrime[0] = lowestPrime[1] = 0;
	for(lli i = 2; i <= n; i++) lowestPrime[i] = (i & 1 ? i : 2);
	lli limit = sqrt(n);
	for(lli i = 3; i <= limit; i += 2){
		if(lowestPrime[i] == i){
			for(lli j = i * i; j <= n; j += 2 * i){
				if(lowestPrime[j] == j) lowestPrime[j] = i;
			}
		}
	}
}

vector<vector<lli>> primeFactors;
void primeFactorsSieve(lli n){
	primeFactors.resize(n + 1, vector<lli>());
	for(int i = 0; i < primes.size(); i++){
		lli p = primes[i];
		for(lli j = p; j <= n; j += p){
			primeFactors[j].push_back(p);
		}
	}
}

vector<lli> Phi;
void phiSieve(lli n){
	Phi.resize(n + 1);
	for(lli i = 1; i <= n; i++) Phi[i] = i;
	for(lli i = 2; i <= n; i ++){
		if(Phi[i] == i){
			for(lli j = i; j <= n; j += i){
				Phi[j] -= Phi[j] / i;
			}
		}
	}
}

vector<vector<lli>> Ncr;
void ncrSieve(lli n){
	Ncr.resize(n + 1, vector<lli>());
	Ncr[0] = {1};
	for(lli i = 1; i <= n; i++){
		Ncr[i].resize(i + 1);
		Ncr[i][0] = Ncr[i][i] = 1;
		for(lli j = 1; j <= i / 2; j++){
			Ncr[i][i - j] = Ncr[i][j] = Ncr[i - 1][j - 1] + Ncr[i - 1][j];
		}
	}
}

vector<pair<lli, int>> factorize(lli n){
	vector<pair<lli, int>> f;
	for(lli & p : primes){
		if(p * p > n) break;
		int pot = 0;
		while(n % p == 0){
			pot++;
			n /= p;
		}
		if(pot) f.push_back(make_pair(p, pot));
	}
	if(n > 1) f.push_back(make_pair(n, 1));
	return f;
}

lli phi(lli n){
	lli ans = n;
	vector<pair<lli, int>> f = factorize(n);
	for(auto & factor : f){
		ans -= ans / factor.first;
	}
	return ans;
}

lli carmichaelLambda(lli n){
	lli ans = 1;
	vector<pair<lli, int>> f = factorize(n);
	for(auto & factor : f){
		lli p = factor.first;
		int a = factor.second;
		lli tmp = pow(p, a);
		tmp -= tmp / p;
		if(a <= 2 || p >= 3) ans = lcm(ans, tmp);
		else ans = lcm(ans, tmp >> 1);
	}
	return ans;
}

int mu(lli n){
	int ans = 1;
	vector<pair<lli, int>> f = factorize(n);
	for(auto & factor : f){
		if(factor.second > 1) return 0;
		ans *= -1;
	}
	return ans;
}

bool testPrimitiveRoot(lli x, lli m){
	if(gcd(x, m) != 1) return false;
	lli order = phi(m);
	vector<pair<lli, int>> f = factorize(order);
	for(auto & factor : f){
		lli p = factor.first;
		if(powMod(x, order / p, m) == 1) return false;
	}
	return true;
}

lli multiplicativeOrder(lli x, lli m){
	if(gcd(x, m) != 1) return -1;
	lli order = phi(m);
	vector<pair<lli, int>> f = factorize(order);
	for(auto & factor : f){
		lli p = factor.first;
		int a = factor.second;
		order /= pow(p, a);
		lli tmp = powMod(x, order, m);
		while(tmp != 1){
			tmp = powMod(tmp, p, m);
			order *= p;
		}
	}
	return order;
}

lli findFirstGenerator(lli m){
	lli order = phi(m);
	vector<pair<lli, int>> f = factorize(order);
	for(lli x = 1; x < m; x++){
		if(gcd(x, m) != 1) continue;
		bool test = true;
		for(auto & factor : f){
			lli p = factor.first;
			if(powMod(x, order / p, m) == 1){
				test = false;
				break;
			}
		}
		if(test) return x;
	}
	return -1;
}

vector<lli> getAllGenerators(lli m){
	vector<lli> ans;
}

// a^x = b mod m, a and m coprime
pair<lli, lli> discreteLogarithm(lli a, lli b, lli m){
	if(gcd(a, m) != 1) return make_pair(-1, 0);
	lli order = multiplicativeOrder(a, m);
	lli n = sqrt(order) + 1;
	lli a_n = powMod(a, n, m);
	lli ans = 0;
	unordered_map<lli, lli> firstHalf;
	set<lli> solutions;
	lli current = a_n;
	for(lli p = 1; p <= n; p++){
		firstHalf[current] = p;
		current = (current * a_n) % m;
	}
	current = b;
	for(lli q = 0; q <= n; q++){
		if(firstHalf.count(current)){
			lli p = firstHalf[current];
			lli x = n * p - q;
			return make_pair(x, order);
		}
		current = (current * a) % m;
	}
	return make_pair(-1, 0);
}

// x^k = b mod m, m prime
vector<lli> discreteRoot(lli k, lli b, lli m){
	if(b == 0) return {0};
	lli g = findFirstGenerator(m);
	lli power = powMod(g, k, m);
	pair<lli, lli> y0 = discreteLogarithm(power, b, m);
	if(y0.first == -1) return {};
	lli phi_m = phi(m);
	lli d = gcd(k, phi_m);
	vector<lli> x(d);
	x[0] = powMod(g, y0.first, m);
	lli inc = powMod(g, phi_m / d, m);
	for(lli i = 1; i < d; i++){
		x[i] = x[i - 1] * inc % m;
	}
	sort(x.begin(), x.end());
	return x;
}

lli potInFactorial(lli n, lli p){
	lli ans = 0;
	lli div = p;
	while(div <= n){
		ans += n / div;
		div *= p;
	}
	return ans;
}

vector<pair<lli, lli>> factorizeFactorial(lli n){
	vector<pair<lli, lli>> f;
	for(lli & p : primes){
		if(p > n) break;
		f.push_back(make_pair(p, potInFactorial(n, p)));
	}
	return f;
}

lli ncr(lli n, lli r){
	if(r < 0 || r > n) return 0;
	r = min(r, n - r);
	lli ans = 1;
	for(lli den = 1, num = n; den <= r; den++, num--){
		ans = ans * num / den;
	}
	return ans;
}

string decimalToBaseB(lli n, lli b){
	string ans = "";
	lli digito;
	do{
		digito = n % b;
		if(0 <= digito && digito <= 9){
			ans = (char)(48 + digito) + ans;
		}else if(10 <= digito && digito <= 35){
			ans = (char)(55 + digito) + ans;
		}
		n /= b;
	}while(n != 0);
	return ans;
}

lli baseBtoDecimal(string & n, lli b){
	lli ans = 0;
	for(char & digito : n){
		if(48 <= digito && digito <= 57){
			ans = ans * b + (digito - 48);
		}else if(65 <= digito && digito <= 90){
			ans = ans * b + (digito - 55);
		}else if(97 <= digito && digito <= 122){
			ans = ans * b + (digito - 87);
		}
	}
	return ans;
}

string decimalToRoman(int n){
	int digito, base = 0;
	string ans = "";
	vector< vector<char> > datos = {{'I', 'V'}, {'X', 'L'}, {'C', 'D'}, {'M', '\0'}};
	int miles = n / 1000;
	do{
		string tmp = "";
		digito = n % 10;
		n /= 10;
		if(base < 3){
			if(0 <= digito && digito <= 3){
				tmp.append(digito, datos[base][0]);
			}else if(digito == 4){
				tmp += datos[base][0];
				tmp += datos[base][1];
			}else if(5 <= digito && digito <= 8){
				tmp += datos[base][1];
				tmp.append(digito - 5, datos[base][0]);
			}else if(digito == 9){
				tmp += datos[base][0];
				tmp += datos[base + 1][0];
			}
		}else{
			tmp.append(miles, 'M');
			ans = tmp + ans;
			break;
		}
		ans = tmp + ans;
		base++;
	}while(n != 0);
	return ans;
}

int romanToDecimal(string n){
	int ans = 0;
	char actual, anterior;
	bool f = false;
	map<char, int> datos = {{'I', 1}, {'V', 5}, {'X', 10}, {'L', 50}, {'C', 100}, {'D', 500}, {'M', 1000}};
	for(int i = n.size() - 1; i >= 0; i--){
		actual = n[i];
		if(i > 0) anterior = n[i - 1];
		if(actual == 'V' && anterior == 'I') ans += 4, f = true;
		else if(actual == 'X' && anterior == 'I') ans += 9, f = true;
		else if(actual == 'L' && anterior == 'X') ans += 40, f = true;
		else if(actual == 'C' && anterior == 'X') ans += 90, f = true;
		else if(actual == 'D' && anterior == 'C') ans += 400, f = true;
		else if(actual == 'M' && anterior == 'C') ans += 900, f = true;
		else{
			if(!f) ans += datos[actual];
			f = false;
		}
	}
	return ans;
}

int mod = 1e9 + 7;

vector<int> P;
int partitionsP(int n){
	if(n < 0) return 0;
	if(P[n]) return P[n];
	int ans = 0, tmp, pos1 = 1, pos2 = 2, inc1 = 4, inc2 = 5;
	for(int k = 1; k <= n; k++){
		tmp = (n >= pos1 ? P[n - pos1] : 0) + (n >= pos2 ? P[n - pos2] : 0);
		if(tmp >= mod) tmp -= mod;
		if(k & 1){
			ans += tmp;
			if(ans >= mod) ans -= mod;
		}else{
			ans -= tmp;
			if(ans < 0) ans += mod;
		}
		if(n < pos2) break;
		pos1 += inc1, pos2 += inc2;
		inc1 += 3, inc2 += 3;
	}
	return ans;
}

void calculateFunctionP(int n){
	P.resize(n + 1);
	P[0] = 1;
	for(int i = 1; i <= n; i++){
		P[i] = partitionsP(i);
	}
}

vector<int> Q;

bool isPerfectSquare(int n){
	int r = sqrt(n);
	return r * r == n;
}

int s(int n){
	int r = 1 + 24 * n;
	if(isPerfectSquare(r)){
		int j;
		r = sqrt(r);
		if((r + 1) % 6 == 0) j = (r + 1) / 6;
		else j = (r - 1) / 6;
		if(j & 1) return -1;
		else return 1;
	}else{
		return 0;
	}
}

int partitionsQ(int n){
	if(n < 0) return 0;
	if(Q[n]) return Q[n];
	int ans = 0, pos = 1, inc = 3;
	int limit = sqrt(n);
	for(int k = 1; k <= limit; k++){
		if(k & 1){
			ans += Q[n - pos];
			if(ans >= mod) ans -= mod;
		}else{
			ans -= Q[n - pos];
			if(ans < 0) ans += mod;
		}
		pos += inc;
		inc += 2;
	}
	ans <<= 1;
	ans += s(n);
	if(ans >= mod) ans -= mod;
	if(ans < 0) ans += mod;
	return ans;
}

void calculateFunctionQ(int n){
	Q.resize(n + 1);
	Q[0] = 1;
	for(int i = 1; i <= n; i++){
		Q[i] = partitionsQ(i);
	}
}

int main(){
	primesSieve(140);
	/*int N = 54;
	for(int i = 1; i < N; i++){
		cout << i << " " << multiplicativeOrder(i, N) << " " << testPrimitiveRoot(i, N) << "\n";
	}*/

	/*lli a, b, m;
	cin >> a >> b >> m;
	pair<lli, lli> ans = discreteLogarithm(a, b, m);
	cout << ans.first << " + " << ans.second << "t";*/

	lli k, b, m;
	cin >> k >> b >> m;
	vector<lli> roots = discreteRoot(k, b, m);
	for(lli & root : roots){
		cout << root << " ";
	}

	/*calculateFunctionP(1e5);
	for(int i = 99900; i <= 100000; i++){
		cout << "P(" << i << ") = " << P[i] << "\n";
	}*/

	/*calculateFunctionQ(1e5);
	for(int i = 99900; i <= 100000; i++){
		cout << "Q(" << i << ") = " << Q[i] << "\n";
	}*/
	return 0;
}