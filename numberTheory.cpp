#include <bits/stdc++.h>
using namespace std;
using lli = __int128_t;

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

lli power(lli b, lli e){
	lli ans = 1;
	while(e){
		if(e & 1) ans *= b;
		e >>= 1;
		b *= b;
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
	if(r0 < 0) s0 *= -1;
	if(s0 < 0) s0 += m;
	return s0;
}

lli powerMod(lli b, lli e, lli m){
	lli ans = 1;
	b %= m;
	if(e < 0){
		b = modularInverse(b, m);
		e *= -1;
	}
	while(e){
		if(e & 1) ans = (ans * b) % m;
		e >>= 1;
		b = (b * b) % m;
	}
	return ans;
}

pair<lli, lli> chinese(vector<lli> & a, vector<lli> & m){
	lli prod = 1, p, ans = 0;
	for(lli & ni : m) prod *= ni;
	for(int i = 0; i < a.size(); ++i){
		p = prod / m[i];
		ans += (a[i] % m[i]) * modularInverse(p, m[i]) % prod * p % prod;
		while(ans >= prod) ans -= prod; while(ans < 0) ans += prod;
	}
	return {ans, prod};
}

vector<lli> divisorsSum;
vector<vector<int>> divisors;
void divisorsSieve(int n){
	divisorsSum.resize(n + 1, 0);
	divisors.resize(n + 1);
	for(int i = 1; i <= n; ++i){
		for(int j = i; j <= n; j += i){
			divisorsSum[j] += i;
			divisors[j].push_back(i);
		}
	}
}

vector<int> primes;
vector<bool> isPrime;
void primesSieve(int n){
	isPrime.resize(n + 1, true);
	isPrime[0] = isPrime[1] = false;
	primes.push_back(2);
	for(int i = 4; i <= n; i += 2) isPrime[i] = false;
	int limit = sqrt(n);
	for(int i = 3; i <= n; i += 2){
		if(isPrime[i]){
			primes.push_back(i);
			if(i <= limit)
				for(int j = i * i; j <= n; j += 2 * i)
					isPrime[j] = false;
		}
	}
}

vector<int> lowestPrime;
void lowestPrimeSieve(int n){
	lowestPrime.resize(n + 1, 1);
	lowestPrime[0] = lowestPrime[1] = 0;
	for(int i = 2; i <= n; ++i) lowestPrime[i] = (i & 1 ? i : 2);
	int limit = sqrt(n);
	for(int i = 3; i <= limit; i += 2)
		if(lowestPrime[i] == i)
			for(int j = i * i; j <= n; j += 2 * i)
				if(lowestPrime[j] == j) lowestPrime[j] = i;
}

vector<vector<int>> primeFactors;
void primeFactorsSieve(lli n){
	primeFactors.resize(n + 1);
	for(int i = 0; i < primes.size(); ++i){
		int p = primes[i];
		for(int j = p; j <= n; j += p)
			primeFactors[j].push_back(p);
	}
}

vector<int> Phi;
void phiSieve(int n){
	Phi.resize(n + 1);
	for(int i = 1; i <= n; ++i) Phi[i] = i;
	for(int i = 2; i <= n; ++i)
		if(Phi[i] == i)
			for(int j = i; j <= n; j += i)
				Phi[j] -= Phi[j] / i;
}

vector<vector<lli>> Ncr;
void ncrSieve(lli n){
	Ncr.resize(n + 1);
	Ncr[0] = {1};
	for(lli i = 1; i <= n; ++i){
		Ncr[i].resize(i + 1);
		Ncr[i][0] = Ncr[i][i] = 1;
		for(lli j = 1; j <= i / 2; j++)
			Ncr[i][i - j] = Ncr[i][j] = Ncr[i - 1][j - 1] + Ncr[i - 1][j];
	}
}

vector<pair<lli, int>> factorize(lli n){
	vector<pair<lli, int>> f;
	for(lli p : primes){
		if(p * p > n) break;
		int pot = 0;
		while(n % p == 0){
			pot++;
			n /= p;
		}
		if(pot) f.emplace_back(p, pot);
	}
	if(n > 1) f.emplace_back(n, 1);
	return f;
}

//divisor power sum of n
//if pot=0 we get the number of divisors
//if pot=1 we get the sum of divisors
lli sigma(lli n, lli pot){
	lli ans = 1;
	auto f = factorize(n);
	for(auto & factor : f){
		lli p = factor.first;
		int a = factor.second;
		if(pot){
			lli p_pot = power(p, pot);
			ans *= (power(p_pot, a + 1) - 1) / (p_pot - 1);
		}else{
			ans *= a + 1;
		}
	}
	return ans;
}

//number of total primes with multiplicity dividing n
int Omega(lli n){
	int ans = 0;
	auto f = factorize(n);
	for(auto & factor : f)
		ans += factor.second;
	return ans;
}

//number of distinct primes dividing n
int omega(lli n){
	int ans = 0;
	auto f = factorize(n);
	for(auto & factor : f)
		++ans;
	return ans;
}

int liouvilleLambda(lli n){
	int exponent = Omega(n);
	return (exponent & 1) ? -1 : 1;
}

//number of coprimes with n less than n
lli phi(lli n){
	lli ans = n;
	auto f = factorize(n);
	for(auto & factor : f)
		ans -= ans / factor.first;
	return ans;
}

//the smallest positive integer k such that for
//every coprime x with n, x^k=1 mod n
lli carmichaelLambda(lli n){
	lli ans = 1;
	auto f = factorize(n);
	for(auto & factor : f){
		lli p = factor.first;
		int a = factor.second;
		lli tmp = power(p, a);
		tmp -= tmp / p;
		if(a <= 2 || p >= 3) ans = lcm(ans, tmp);
		else ans = lcm(ans, tmp >> 1);
	}
	return ans;
}

//1 if n is square-free with an even number of prime factors
//-1 if n is square-free with an odd number of prime factors
//0 is n has a square prime factor
int mu(lli n){
	int ans = 1;
	auto f = factorize(n);
	for(auto & factor : f){
		if(factor.second > 1) return 0;
		ans *= -1;
	}
	return ans;
}

// the smallest positive integer k such that x^k = 1 mod m
lli multiplicativeOrder(lli x, lli m){
	if(gcd(x, m) != 1) return 0;
	lli order = phi(m);
	auto f = factorize(order);
	for(auto & factor : f){
		lli p = factor.first;
		int a = factor.second;
		order /= power(p, a);
		lli tmp = powerMod(x, order, m);
		while(tmp != 1){
			tmp = powerMod(tmp, p, m);
			order *= p;
		}
	}
	return order;
}

//number of generators modulo m
lli numberOfGenerators(lli m){
	lli phi_m = phi(m);
	lli lambda_m = carmichaelLambda(m);
	if(phi_m == lambda_m) return phi(phi_m);
	else return 0;
}

//test if order(x, m) = phi(m), i.e., x is a generator for Z/mZ
bool testPrimitiveRoot(lli x, lli m){
	if(gcd(x, m) != 1) return false;
	lli order = phi(m);
	auto f = factorize(order);
	for(auto & factor : f){
		lli p = factor.first;
		if(powerMod(x, order / p, m) == 1) return false;
	}
	return true;
}

//test if x^k = 1 mod m and k is the smallest for such x, i.e., x^(k/p) != 1 for every prime divisor of k
bool testPrimitiveKthRootUnity(lli x, lli k, lli m){
	if(powerMod(x, k, m) != 1) return false;
	auto f = factorize(k);
	for(auto & factor : f){
		lli p = factor.first;
		if(powerMod(x, k / p, m) == 1) return false;
	}
	return true;
}

lli findFirstGenerator(lli m){
	lli order = phi(m);
	if(order != carmichaelLambda(m)) return -1; //just an optimization, not required
	auto f = factorize(order);
	for(lli x = 1; x < m; x++){
		if(gcd(x, m) != 1) continue;
		bool test = true;
		for(auto & factor : f){
			lli p = factor.first;
			if(powerMod(x, order / p, m) == 1){
				test = false;
				break;
			}
		}
		if(test) return x;
	}
	return -1; //not found
}

lli findFirstPrimitiveKthRootUnity(lli k, lli m){
	if(carmichaelLambda(m) % k != 0) return -1; //just an optimization, not required
	auto f = factorize(k);
	for(lli x = 1; x < m; x++){
		if(powerMod(x, k, m) != 1) continue;
		bool test = true;
		for(auto & factor : f){
			lli p = factor.first;
			if(powerMod(x, k / p, m) == 1){
				test = false;
				break;
			}
		}
		if(test) return x;
	}
	return -1; //not found
}

// Solves for x in the equation a^x = b mod m
pair<lli, lli> discreteLogarithm(lli a, lli b, lli m){
	lli m1 = m, pw = 1, d, x, y, nonRep = 0;
	for(; (d = gcd(a, m1)) > 1; ++nonRep, m1 /= d, pw = pw * a % m){
		if(pw == b) return {nonRep, 0}; //aperiodic solution found
	}
	d = extendedGcd(pw, m, x, y);
	if(b % d > 0) return {-1, 0}; //solution not found
	b = x * (b / d) % m;
	if(b < 0) b += m;
	lli order = multiplicativeOrder(a, m1);
	lli n = sqrt(order) + 1;
	lli a_n = powerMod(a, n, m1);
	unordered_map<lli, lli> firstHalf;
	pw = a_n;
	for(lli p = 1; p <= n; ++p, pw = pw * a_n % m1){
		firstHalf[pw] = p;
	}
	pw = b % m1;
	for(lli q = 0; q <= n; ++q, pw = pw * a % m1){
		if(firstHalf.count(pw)) return {nonRep + (n * firstHalf[pw] - q) % order, order}; //periodic solution found
	}
	return {-1, 0}; //solution not found
}

// x^k = b mod m, m has at least one generator
vector<lli> discreteRoot(lli k, lli b, lli m){
	if(b % m == 0) return {0};
	lli g = findFirstGenerator(m);
	lli power = powerMod(g, k, m);
	auto y0 = discreteLogarithm(power, b, m);
	if(y0.first == -1) return {};
	lli phi_m = phi(m);
	lli d = gcd(k, phi_m);
	vector<lli> x(d);
	x[0] = powerMod(g, y0.first, m);
	lli inc = powerMod(g, phi_m / d, m);
	for(lli i = 1; i < d; i++)
		x[i] = x[i - 1] * inc % m;
	sort(x.begin(), x.end());
	return x;
}

lli potInFactorial(lli n, lli p){
	lli ans = 0, div = n;
	while(div /= p) ans += div;
	return ans;
}

vector<pair<lli, lli>> factorizeFactorial(lli n){
	vector<pair<lli, lli>> f;
	for(lli p : primes){
		if(p > n) break;
		f.emplace_back(p, potInFactorial(n, p));
	}
	return f;
}

lli ncr(lli n, lli r){
	if(r < 0 || r > n) return 0;
	r = min(r, n - r);
	lli ans = 1;
	for(lli den = 1, num = n; den <= r; den++, num--)
		ans = ans * num / den;
	return ans;
}

string decimalToBaseB(lli n, lli b){
	string ans = "";
	lli d;
	do{
		d = n % b;
		if(0 <= d && d <= 9) ans = (char)(48 + d) + ans;
		else if(10 <= d && d <= 35) ans = (char)(55 + d) + ans;
		n /= b;
	}while(n != 0);
	return ans;
}

lli baseBtoDecimal(const string & n, lli b){
	lli ans = 0;
	for(const char & d : n){
		if(48 <= d && d <= 57) ans = ans * b + (d - 48);
		else if(65 <= d && d <= 90) ans = ans * b + (d - 55);
		else if(97 <= d && d <= 122) ans = ans * b + (d - 87);
	}
	return ans;
}

string decimalToRoman(int n){
	int d, b = 0;
	string ans = "";
	vector<vector<char>> datos = {{'I', 'V'}, {'X', 'L'}, {'C', 'D'}, {'M', '\0'}};
	int miles = n / 1000;
	do{
		string tmp = "";
		d = n % 10;
		n /= 10;
		if(b < 3){
			if(0 <= d && d <= 3){
				tmp.append(d, datos[b][0]);
			}else if(d == 4){
				tmp += datos[b][0];
				tmp += datos[b][1];
			}else if(5 <= d && d <= 8){
				tmp += datos[b][1];
				tmp.append(d - 5, datos[b][0]);
			}else if(d == 9){
				tmp += datos[b][0];
				tmp += datos[b + 1][0];
			}
		}else{
			tmp.append(miles, 'M');
			ans = tmp + ans;
			break;
		}
		ans = tmp + ans;
		b++;
	}while(n != 0);
	return ans;
}

int romanToDecimal(string n){
	int ans = 0;
	char curr, prev;
	bool f = false;
	map<char, int> datos = {{'I', 1}, {'V', 5}, {'X', 10}, {'L', 50}, {'C', 100}, {'D', 500}, {'M', 1000}};
	for(int i = n.size() - 1; i >= 0; i--){
		curr = n[i];
		if(i > 0) prev = n[i - 1];
		if(curr == 'V' && prev == 'I') ans += 4, f = true;
		else if(curr == 'X' && prev == 'I') ans += 9, f = true;
		else if(curr == 'L' && prev == 'X') ans += 40, f = true;
		else if(curr == 'C' && prev == 'X') ans += 90, f = true;
		else if(curr == 'D' && prev == 'C') ans += 400, f = true;
		else if(curr == 'M' && prev == 'C') ans += 900, f = true;
		else{
			if(!f) ans += datos[curr];
			f = false;
		}
	}
	return ans;
}

lli mod = 1e9 + 7;

vector<lli> P;

//number of ways to write n as a sum of positive integers
lli partitionsP(int n){
	if(n < 0) return 0;
	if(P[n]) return P[n];
	int pos1 = 1, pos2 = 2, inc1 = 4, inc2 = 5;
	lli ans = 0;
	for(int k = 1; k <= n; k++){
		lli tmp = (n >= pos1 ? P[n - pos1] : 0) + (n >= pos2 ? P[n - pos2] : 0);
		if(k & 1) ans += tmp;
		else ans -= tmp;
		if(n < pos2) break;
		pos1 += inc1, pos2 += inc2;
		inc1 += 3, inc2 += 3;
	}
	ans %= mod;
	if(ans < 0) ans += mod;
	return ans;
}

void calculateFunctionP(int n){
	P.resize(n + 1);
	P[0] = 1;
	for(int i = 1; i <= n; i++)
		P[i] = partitionsP(i);
}

vector<lli> Q;

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

//number of ways to write n as a sum of distinct positive integers
//number of ways to write n as a sum of odd positive integers
lli partitionsQ(int n){
	if(n < 0) return 0;
	if(Q[n]) return Q[n];
	int pos = 1, inc = 3;
	lli ans = 0;
	int limit = sqrt(n);
	for(int k = 1; k <= limit; k++){
		if(k & 1) ans += Q[n - pos];
		else ans -= Q[n - pos];
		pos += inc;
		inc += 2;
	}
	ans <<= 1;
	ans += s(n);
	ans %= mod;
	if(ans < 0) ans += mod;
	return ans;
}

void calculateFunctionQ(int n){
	Q.resize(n + 1);
	Q[0] = 1;
	for(int i = 1; i <= n; i++)
		Q[i] = partitionsQ(i);
}

//continued fraction of (p+sqrt(n))/q, where p,n,q are positive integers
//returns a vector of terms and the length of the period,
//the periodic part is taken from the right of the array
pair<vector<lli>, int> ContinuedFraction(lli p, lli n, lli q){
	vector<lli> coef;
	lli r = sqrt(n);
	//Skip this if you know that n is not a perfect square
	if(r * r == n){
		lli num = p + r;
		lli den = q;
		lli residue;
		while(den){
			residue = num % den;
			coef.push_back(num / den);
			num = den;
			den = residue;
		}
		return {coef, 0};
	}
	if((n - p * p) % q != 0){
		n *= q * q;
		p *= q;
		q *= q;
		r = sqrt(n);
	}
	lli a = (r + p) / q;
	coef.push_back(a);
	int period = 0;
	map<pair<lli, lli>, int> pairs;
	while(true){
		p = a * q - p;
		q = (n - p * p) / q;
		a = (r + p) / q;
		//if p=0 and q=1, we can just ask if q==1 after inserting a
		if(pairs.count({p, q})){
			period -= pairs[{p, q}];
			break;
		}
		coef.push_back(a);
		pairs[{p, q}] = period++;
	}
	return {coef, period};
}

//first solution (x, y) to the equation x^2-ny^2=1, n IS NOT a perfect aquare
pair<lli, lli> PellEquation(lli n){
	vector<lli> cf = ContinuedFraction(0, n, 1).first;
	lli num = 0, den = 1;
	int k = cf.size() - 1;
	for(int i = ((k & 1) ? (2 * k - 1) : (k - 1)); i >= 0; i--){
		lli tmp = den;
		int pos = i % k;
		if(pos == 0 && i != 0) pos = k;
		den = num + cf[pos] * den;
		num = tmp;
	}
	return {den, num};
}

bool isPrimeMillerRabin(lli n){
	if(n < 2) return false;
	if(!(n & 1)) return n == 2;
	lli d = n - 1, s = 0;
	for(; !(d & 1); d >>= 1, ++s);
	for(int a : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}){
		if(n == a) return true;
		lli m = powerMod(a, d, n);
		if(m == 1 || m == n - 1) continue;
		int k = 0;
		for(; k < s; ++k){
			m = m * m % n;
			if(m == n - 1) break;
		}
		if(k == s) return false;
	}
	return true;
}

mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
lli aleatorio(lli a, lli b){
	std::uniform_int_distribution<lli> dist(a, b);
	return dist(rng);
}
lli getFactor(lli n){
	lli a = aleatorio(1, n - 1), b = aleatorio(1, n - 1);
	lli x = 2, y = 2, d = 1;
	while(d == 1){
		x = x * (x + b) % n + a;
		y = y * (y + b) % n + a;
		y = y * (y + b) % n + a;
		d = gcd(abs(x - y), n);
	}
	return d;
}

map<lli, int> fact;
void factorizePollardRho(lli n, bool clean = true){
	if(clean) fact.clear();
	while(n > 1 && !isPrimeMillerRabin(n)){
		lli f = n;
		for(; f == n; f = getFactor(n));
		n /= f;
		factorizePollardRho(f, false);
		for(auto & it : fact){
			while(n % it.first == 0){
				n /= it.first;
				++it.second;
			}
		}
	}
	if(n > 1) ++fact[n];
}

//find all inverses (from 1 to p-1) modulo p
vector<lli> allInverses(lli p){
	vector<lli> ans(p);
	ans[1] = 1;
	for(lli i = 2; i < p; ++i)
		ans[i] = p - (p / i) * ans[p % i] % p;
	return ans;
}

//very fast fibonacci
inline void modula(lli & n, lli mod){
	while(n >= mod) n -= mod;
}

lli fibo(lli n, lli mod){
	array<lli, 2> F = {1, 0};
	lli p = 1;
	for(lli v = n; v >>= 1; p <<= 1);
	array<lli, 4> C;
	do{
		int d = (n & p) != 0;
		C[0] = C[3] = 0;
		C[d] = F[0] * F[0] % mod;
		C[d+1] = (F[0] * F[1] << 1) % mod;
		C[d+2] = F[1] * F[1] % mod;
		F[0] = C[0] + C[2] + C[3];
		F[1] = C[1] + C[2] + (C[3] << 1);
		modula(F[0], mod), modula(F[1], mod);
	}while(p >>= 1);
	return F[1];
}

//number of ordered factorizations of n
lli orderedFactorizations(lli n){
	//skip the factorization if you already know the powers
	auto fact = factorize(n);
	int k = 0, q = 0;
	vector<int> powers(fact.size() + 1);
	for(auto & f : fact){
		powers[k + 1] = f.second;
		q += f.second;
		++k;
	}
	vector<lli> prod(q + 1, 1);
	//we need Ncr until the max_power+Omega(n) row
	//module if needed
	for(int i = 0; i <= q; i++){
		for(int j = 1; j <= k; j++){
			prod[i] = prod[i] * Ncr[powers[j] + i][powers[j]];
		}
	}
	lli ans = 0;
	for(int j = 1; j <= q; j++){
		int alt = 1;
		for(int i = 0; i < j; i++){
			ans = ans + alt * Ncr[j][i] * prod[j - i - 1];
			alt *= -1;
		}
	}
	return ans;
}

//Number of unordered factorizations of n with
//largest part at most m
//Call unorderedFactorizations(n, n) to get all of them
//Add this to the main to speed up the map:
//mem.reserve(1024); mem.max_load_factor(0.25);
struct HASH{
  size_t operator()(const pair<int,int>&x)const{
    return hash<long long>()(((long long)x.first)^(((long long)x.second)<<32));
  }
};
unordered_map<pair<int, int>, lli, HASH> mem;
lli unorderedFactorizations(int m, int n){
	if(m == 1 && n == 1) return 1;
	if(m == 1) return 0;
	if(n == 1) return 1;
	if(mem.count({m, n})) return mem[{m, n}];
	lli ans = 0;
	int l = sqrt(n);
	for(int i = 1; i <= l; ++i){
		if(n % i == 0){
			int a = i, b = n / i;
			if(a <= m) ans += unorderedFactorizations(a, b);
			if(a != b && b <= m) ans += unorderedFactorizations(b, a);
		}
	}
	return mem[{m, n}] = ans;
}

vector<int> Mu;
void muSieve(int n){
	Mu.resize(n + 1, -1);
	Mu[0] = 0, Mu[1] = 1;
	for(int i = 2; i <= n; ++i)
		if(Mu[i])
			for(int j = 2*i; j <= n; j += i)
				Mu[j] -= Mu[i];
}

vector<int> linearPrimeSieve(int n){
	vector<int> primes;
	vector<bool> isPrime(n+1, true);
	for(int i = 2; i <= n; ++i){
		if(isPrime[i])
			primes.push_back(i);
		for(int p : primes){
			int d = i * p;
			if(d > n) break;
			isPrime[d] = false;
			if(i % p == 0) break;
		}
	}
	return primes;
}

//suppose f(n) is a multiplicative function and
//we want to find f(1), f(2), ..., f(n)
//we have f(pq) = f(p)f(q) if gcd(p, q) = 1
//and f(p^a) = g(p, a), where p is prime and a>0
vector<int> generalSieve(int n, function<int(int, int)> g){
	vector<int> f(n+1, 1), cnt(n+1), acum(n+1), primes;
	vector<bool> isPrime(n+1, true);
	for(int i = 2; i <= n; ++i){
		if(isPrime[i]){ //case base: f(p)
			f[i] = g(i, 1);
			primes.push_back(i);
			cnt[i] = 1;
			acum[i] = i;
		}
		for(int p : primes){
			int d = i * p;
			if(d > n) break;
			isPrime[d] = false;
			if(i % p == 0){ //gcd(i, p) != 1
				f[d] = f[i / acum[i]] * g(p, cnt[i] + 1);
				cnt[d] = cnt[i] + 1;
				acum[d] = acum[i] * p;
				break;
			}else{ //gcd(i, p) = 1
				f[d] = f[i] * g(p, 1);
				cnt[d] = 1;
				acum[d] = p;
			}
		}
	}
	return f;
}

vector<int> segmented_sieve(int limit){
	const int L1D_CACHE_SIZE = 32768;
	int raiz = sqrt(limit);
	int segment_size = max(raiz, L1D_CACHE_SIZE);
	int s = 3, n = 3;
	vector<int> primes(1, 2), tmp, next;
	vector<char> sieve(segment_size);
	vector<bool> is_prime(raiz + 1, 1);
	for(int i = 2; i * i <= raiz; i++)
		if(is_prime[i])
			for(int j = i * i; j <= raiz; j += i)
				is_prime[j] = 0;
	for(int low = 0; low <= limit; low += segment_size){
		fill(sieve.begin(), sieve.end(), 1);
		int high = min(low + segment_size - 1, limit);
		for(; s * s <= high; s += 2){
			if(is_prime[s]){
				tmp.push_back(s);
				next.push_back(s * s - low);
			}
		}
		for(size_t i = 0; i < tmp.size(); i++){
			int j = next[i];
			for(int k = tmp[i] * 2; j < segment_size; j += k)
				sieve[j] = 0;
			next[i] = j - segment_size;
		}
		for(; n <= high; n += 2)
			if(sieve[n - low])
				primes.push_back(n);
	}
	return primes;
}

//number of ways to partition a set of n elements
//the nth bell number is at Bell[n][0]
vector<vector<int>> Bell;
void bellNumbers(int n){
	Bell.resize(n + 1);
	Bell[0] = {1};
	for(int i = 1; i <= n; ++i){
		Bell[i].resize(i + 1);
		Bell[i][0] = Bell[i - 1][i - 1];
		for(int j = 1; j <= i; ++j)
			Bell[i][j] = Bell[i][j - 1] + Bell[i - 1][j - 1];
	}
}

//finds x such that x^2 = a mod p
lli sqrtMod(lli a, lli p){
	a %= p;
	if(a < 0) a += p;
	if(a == 0) return 0;
	assert(powerMod(a, (p - 1) / 2, p) == 1);
	if(p % 4 == 3) return powerMod(a, (p + 1) / 4, p);
	lli s = p - 1;
	int r = 0;
	while((s & 1) == 0) ++r, s >>= 1;
	lli n = 2;
	while(powerMod(n, (p - 1) / 2, p) != p - 1) ++n;
	lli x = powerMod(a, (s + 1) / 2, p);
	lli b = powerMod(a, s, p);
	lli g = powerMod(n, s, p);
	while(true){
		lli t = b;
		int m = 0;
		for(; m < r; ++m){
			if(t == 1) break;
			t = t * t % p;
		}
		if(m == 0) return x;
		lli gs = powerMod(g, 1 << (r - m - 1), p);
		g = gs * gs % p;
		x = x * gs % p;
		b = b * g % p;
		r = m;
	}
}

vector<int> greatestPrime;
void greatestPrimeSieve(int n){
	greatestPrime.resize(n + 1, 1);
	greatestPrime[0] = greatestPrime[1] = 0;
	for(int i = 2; i <= n; ++i) greatestPrime[i] = i;
	for(int i = 2; i <= n; i++)
		if(greatestPrime[i] == i)
			for(int j = i; j <= n; j += i)
				greatestPrime[j] = i;
}

//generalized chinese remainder theorem
//the modulos doesn't need to be pairwise coprime
pair<lli, lli> crt(const vector<lli> & a, const vector<lli> & m){
	lli a0 = a[0] % m[0], m0 = m[0], a1, m1, s, t, d, M;
	for(int i = 1; i < a.size(); ++i){
		a1 = a[i] % m[i], m1 = m[i];
		d = extendedGcd(m0, m1, s, t);
		if((a0 - a1) % d != 0) return {0, 0}; //error, no solution
		M = m0 * (m1 / d);
		a0 = a0 * t % M * (m1 / d) % M + a1 * s % M * (m0 / d) % M;
		while(a0 >= M) a0 -= M; while(a0 < 0) a0 += M;
		m0 = M;
	}
	while(a0 >= m0) a0 -= m0; while(a0 < 0) a0 += m0;
	return {a0, m0};
}

//s(n, k) represents the number of permutations
//of n elements with k disjoint cycles
vector<vector<lli>> stirling1;
void stirlingNumber1stKind(lli n){
	stirling1.resize(n+1, vector<lli>(n+1));
	stirling1[0][0] = 1;
	for(int i = 1; i <= n; ++i)
		for(int j = 1; j <= i; ++j)
			stirling1[i][j] = (i-1) * stirling1[i-1][j] + stirling1[i-1][j-1];
}

//S(n, k) represents the number of ways to
//partition a set of n object into k non-empty
//distinct subsets
vector<vector<lli>> stirling2;
void stirlingNumber2ndKind(lli n){
	stirling2.resize(n+1, vector<lli>(n+1));
	stirling2[0][0] = 1;
	for(int i = 1; i <= n; ++i)
		for(int j = 1; j <= i; ++j)
			stirling2[i][j] = j * stirling2[i-1][j] + stirling2[i-1][j-1];
}

//euler(n, k) represents the number of permutations
//of 1,...,n with exactly k numbers greater than
//the previous number
vector<vector<lli>> euler;
void eulerianNumbers(lli n){
	euler.resize(n+1, vector<lli>(n+1));
	for(int i = 1; i <= n; ++i){
		euler[i][0] = 1;
		for(int j = 1; j < i; ++j)
			euler[i][j] = (i-j) * euler[i-1][j-1] + (j+1) * euler[i-1][j];
	}
}

//finds sum(floor(p*i/q), 1<=i<=n)
lli floorsSum(lli p, lli q, lli n){
    lli t = gcd(p, q);
    p /= t, q /= t;
    lli s = 0, z = 1;
    while(q && n){
        t = p/q;
        s += z*t*n*(n+1)/2;
        p -= q*t;
        t = n/q;
        s += z*p*t*(n+1) - z*t*(p*q*t + p + q - 1)/2;
        n -= q*t;
        t = n*p/q;
        s += z*t*n;
        n = t;
        swap(p, q);
        z = -z;
    }
    return s;
}

//Finds (n!/p^m) mod p^s, where m is the largest power of p
//that divides n!, p must be prime
lli factmod(lli n, lli p, int s){
	lli ans = 1;
	lli ps = power(p, s);
	while(n > 1){
		lli q = n / ps, r = n % ps;
		ans = ans * (q % 2 == 1 && !(p == 2 && s >= 3) ? ps-1 : 1) % ps;
		for(lli i = 2; i <= r; ++i){
			if(i % p == 0) continue;
			ans = ans * i % ps;
		}
		n /= p;
	}
	return ans;
}

uint64_t mul_mod(uint64_t a, uint64_t b, uint64_t m){
	if(a >= m) a %= m;
	if(b >= m) b %= m;
	uint64_t c = (long double)a * b / m;
	int64_t r = (int64_t)(a * b - c * m) % (int64_t)m;
	return r < 0 ? r + m : r;
}

lli pisano_prime(lli p){
	if(p == 2) return 3;
	if(p == 5) return 20;
	lli order = 0;
	if(p%10 == 1 || p%10 == 9) order = p - 1;
	else order = 2*p + 2;
	auto fact = factorize(order);
	for(auto par : fact){
		lli q; int a;
		tie(q, a) = par;
		order /= power(q, a);
		while(!(fibo(order, p) == 0 && fibo(order+1, p) == 1)){
			order *= q;
		}
	}
	return order;
}

lli pisano(lli mod){
	lli ans = 1;
	auto fact = factorize(mod);
	for(auto par : fact){
		lli p; int a;
		tie(p, a) = par;
		ans = lcm(ans, power(p, a-1) * pisano_prime(p));
	}
	return ans;
}

ostream &operator<<(ostream &os, const __int128 & value){
	char buffer[64];
	char *pos = end(buffer) - 1;
	*pos = '\0';
	__int128 tmp = value < 0 ? -value : value;
	do{
		--pos;
		*pos = tmp % 10 + '0';
		tmp /= 10;
	}while(tmp != 0);
	if(value < 0){
		--pos;
		*pos = '-';
	}
	return os << pos;
}

istream &operator>>(istream &is, __int128 & value){
	char buffer[64];
	is >> buffer;
	char *pos = begin(buffer);
	int sgn = 1;
	value = 0;
	if(*pos == '-'){
		sgn = -1;
		++pos;
	}else if(*pos == '+'){
		++pos;
	}
	while(*pos != '\0'){
		value = (value << 3) + (value << 1) + (*pos - '0');
		++pos;
	}
	value *= sgn;
	return is;
}

void info_ntt(lli e, lli l, lli r){
	lli n = (lli)1 << e;
	for(lli k = l; k <= r; ++k){
		lli p = k * n + 1;
		if(isPrimeMillerRabin(p)){
			cout << "g = " << findFirstGenerator(p) << "\n";
			lli w = findFirstPrimitiveKthRootUnity(n, p);
			cout << "w = " << w << "\nw^-1 = " << modularInverse(w, p) << "\nk = " << k << "\nn = " << n << "\np = " << p << "\n\n";
		}
	}
}

int main(){
	srand(time(NULL));
	primesSieve(1e5);
	ncrSieve(300);
	muSieve(300);
	bellNumbers(50);
	stirlingNumber1stKind(50);
	stirlingNumber2ndKind(50);
	eulerianNumbers(50);

	/*lli sum = 0, p, q, n;
	cin >> p >> q >> n;
	for(lli i = 1; i <= n; ++i) sum += p*i/q;
	assert(sum == floorsSum(p, q, n));
	cout << sum << "\n";*/

	/*for(int i = 0; i <= 10; ++i){
		for(int j = 0; j <= i; ++j){
			cout << euler[i][j] << " ";
		}
		cout << "\n";
	}*/

	/*int l = 1e7;
	clock_t start = clock();
	for(int i = 1; i <= l; ++i) fibo(i, mod);
	cout << fixed << setprecision(4) << (clock() - start) / (double)CLOCKS_PER_SEC << "\n";*/

	/*clock_t start = clock();
	vector<int> ans = segmented_sieve(1e9);
	cout << ans.size() << "\n" << ans.back() << "\n";
	cout << fixed << setprecision(4) << (clock() - start)/(double)CLOCKS_PER_SEC << "s\n";*/

	/*auto g = [&](int p, int a){
		return power(p, a) - power(p, a-1);
	};
	vector<int> pre = generalSieve(105, g);
	for(int i = 1; i <= 105; ++i){
		cout << "f(" << i << ") = " << pre[i] << "\n";
	}*/

	/*for(int i = 1; i <= 300; ++i){
		cout << "mu(" << i << ") = " << Mu[i] << "\n";
	}*/

	/*lli x;
	cin >> x;
	cout << isPrimeMillerRabin(x) << "\n";
	factorizePollardRho(x);
	for(auto & it : fact){
		cout << it.first << " " << it.second << "\n";
	}*/

	for(lli i = 1; i <= 10000; ++i){
		cout << i << " " << pisano(i) << "\n";
	}

	/*lli p, n, q;
	cin >> p >> n >> q;
	auto cf = ContinuedFraction(p, n, q);
	auto eq = PellEquation(n);
	cout << "Period: " << cf.second << "\n";
	for(lli & coef : cf.first) cout << coef << " ";
	cout << "\n" << eq.first << ", " << eq.second;*/
	/*int N = 54;
	for(int i = 1; i < N; i++){
		cout << i << " " << multiplicativeOrder(i, N) << " " << testPrimitiveRoot(i, N) << "\n";
	}*/

	/*int N = 7340033, k = 1 << 20;
	for(int i = 1; i < 100; i++){
		cout << i << " " << testPrimitiveKthRootUnity(i, k, N) << "\n";
	}*/

	/*long long int k, m;
	cin >> k >> m;
	cout << (long long int)findFirstPrimitiveKthRootUnity(k, m);*/

	/*for(lli n = 1; n <= 100; ++n){
		cout << n << " " << factmod(n, 7, 2) << "\n";
	}*/

	/*lli a, b, m;
	cin >> a >> b >> m;
	pair<lli, lli> ans = discreteLogarithm(a, b, m);
	cout << ans.first << " + " << ans.second << "t";*/

	/*lli k, b, m;
	cin >> k >> b >> m;
	vector<lli> roots = discreteRoot(k, b, m);
	for(lli & root : roots){
		cout << root << " ";
	}*/

	/*calculateFunctionP(1e5);
	for(int i = 99900; i <= 100000; i++){
		cout << "P(" << i << ") = " << P[i] << "\n";
	}*/

	/*calculateFunctionQ(1e5);
	for(int i = 99900; i <= 100000; i++){
		cout << "Q(" << i << ") = " << Q[i] << "\n";
	}*/

	/*lli e;
	cin >> e;
	info_ntt(e, 1, 1000);*/

	/*lli p;
	cin >> p;
	vector<lli> inverses = allInverses(p);
	for(lli i = 1; i < p; ++i){
		cout << i << "^-1 mod " << p << " = " << inverses[i] << "\n";
	}*/

	/*lli n;
	cin >> n;
	cout << unorderedFactorizations(n, n);*/

	/*lli a, p;
	cin >> a >> p;
	lli ans = sqrtMod(a, p);
	cout << ans << "^2 = " << ans*ans%p << " mod " << p << "\n";*/

	
	/*int n;
	cin >> n;
	vector<lli> a(n), m(n);
	for(int i = 0; i < n; ++i){
		cin >> a[i] >> m[i];
	}
	auto answer = crt(a, m);
	cout << "x = " << answer.first << " mod " << answer.second << "\n";*/
	
	return 0;
}