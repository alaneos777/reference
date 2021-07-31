#include <bits/stdc++.h>
using namespace std;
using lli = long long int;
using comp = complex<double>;
using poly = vector<int>;
const double PI = acos(-1.0);
int nearestPowerOfTwo(int n){
	int ans = 1;
	while(ans < n) ans <<= 1;
	return ans;
}

void fft(vector<comp> & X, int inv){
	int n = X.size();
	for(int i = 1, j = 0; i < n - 1; ++i){
		for(int k = n >> 1; (j ^= k) < k; k >>= 1);
		if(i < j) swap(X[i], X[j]);
	}
	vector<comp> wp(n>>1);
	for(int k = 1; k < n; k <<= 1){
		for(int j = 0; j < k; ++j)
			wp[j] = polar(1.0, PI * j / k * inv); // best precision but slower
		for(int i = 0; i < n; i += k << 1){
			for(int j = 0; j < k; ++j){
				comp t = X[i + j + k] * wp[j];
				X[i + j + k] = X[i + j] - t;
				X[i + j] += t;
			}
		}
	}
	if(inv == -1)
		for(int i = 0; i < n; ++i)
			X[i] /= n;
}

lli powerMod(lli b, lli e, lli m){
	lli ans = 1;
	e %= m-1;
	if(e < 0) e += m-1;
	while(e){
		if(e & 1) ans = ans * b % m;
		e >>= 1;
		b = b * b % m;
	}
	return ans;
}

template<int p, int g>
void ntt(poly & X, int inv){
	int n = X.size();
	for(int i = 1, j = 0; i < n - 1; ++i){
		for(int k = n >> 1; (j ^= k) < k; k >>= 1);
		if(i < j) swap(X[i], X[j]);
	}
	vector<lli> wp(n>>1, 1);
	for(int k = 1; k < n; k <<= 1){
		lli wk = powerMod(g, inv * (p - 1) / (k<<1), p);
		for(int j = 1; j < k; ++j)
			wp[j] = wp[j - 1] * wk % p;
		for(int i = 0; i < n; i += k << 1){
			for(int j = 0; j < k; ++j){
				int u = X[i + j], v = X[i + j + k] * wp[j] % p;
				X[i + j] = u + v < p ? u + v : u + v - p;
				X[i + j + k] = u - v < 0 ? u - v + p : u - v;
			}
		}
	}
	if(inv == -1){
		lli nrev = powerMod(n, p - 2, p);
		for(int i = 0; i < n; ++i)
			X[i] = X[i] * nrev % p;
	}
}

vector<comp> convolution(vector<comp> A, vector<comp> B){
	int sz = A.size() + B.size() - 1;
	int size = nearestPowerOfTwo(sz);
	A.resize(size), B.resize(size);
	fft(A, 1), fft(B, 1);
	for(int i = 0; i < size; i++)
		A[i] *= B[i];
	fft(A, -1);
	A.resize(sz);
	return A;
}

template<int p, int g>
poly convolution(poly A, poly B){
	int sz = A.size() + B.size() - 1;
	int size = nearestPowerOfTwo(sz);
	A.resize(size), B.resize(size);
	ntt<p, g>(A, 1), ntt<p, g>(B, 1);
	for(int i = 0; i < size; i++)
		A[i] = (lli)A[i] * B[i] % p;
	ntt<p, g>(A, -1);
	A.resize(sz);
	return A;
}

const int p = 7340033, g = 3; //default values for NTT

string multiplyNumbers(const string & a, const string & b){
	int sgn = 1;
	int pos1 = 0, pos2 = 0;
	while(pos1 < a.size() && (a[pos1] < '1' || a[pos1] > '9')){
		if(a[pos1] == '-') sgn *= -1;
		++pos1;
	}
	while(pos2 < b.size() && (b[pos2] < '1' || b[pos2] > '9')){
		if(b[pos2] == '-') sgn *= -1;
		++pos2;
	}
	poly X(a.size() - pos1), Y(b.size() - pos2);
	if(X.empty() || Y.empty()) return "0";
	for(int i = pos1, j = X.size() - 1; i < a.size(); ++i)
		X[j--] = a[i] - '0';
	for(int i = pos2, j = Y.size() - 1; i < b.size(); ++i)
		Y[j--] = b[i] - '0';
	X = convolution<p, g>(X, Y);
	stringstream ss;
	if(sgn == -1) ss << "-";
	int carry = 0;
	for(int i = 0; i < X.size(); ++i){
		X[i] += carry;
		carry = X[i] / 10;
		X[i] %= 10;
	}
	while(carry){
		X.push_back(carry % 10);
		carry /= 10;
	}
	for(int i = X.size() - 1; i >= 0; --i)
		ss << X[i];
	return ss.str();
}

poly inversePolynomial(const poly & A){
	poly R(1, powerMod(A[0], p - 2, p));
	//R(x) = 2R(x)-A(x)R(x)^2
	while(R.size() < A.size()){
		size_t c = 2 * R.size();
		R.resize(c);
		poly R2 = R;
		poly a(min(c, A.size()));
		for(int i = 0; i < a.size(); ++i)
			a[i] = A[i];
		R2 = convolution<p, g>(R2, R2);
		R2.resize(c);
		R2 = convolution<p, g>(R2, a);
		for(int i = 0; i < c; ++i){
			R[i] = R[i] + R[i] - R2[i];
			if(R[i] < 0) R[i] += p;
			if(R[i] >= p) R[i] -= p;
		}
	}
	R.resize(A.size());
	return R;
}

const int inv2 = powerMod(2, p - 2, p);

poly sqrtPolynomial(const poly & A){
	int r0 = 1; //verify that r0^2 = A[0] mod p
	poly R(1, r0);
	//R(x) = R(x)/2 + A(x)/(2R(x))
	while(R.size() < A.size()){
		size_t c = 2 * R.size();
		R.resize(c);
		poly a(min(c, A.size()));
		for(int i = 0; i < a.size(); ++i)
			a[i] = A[i];
		a = convolution<p, g>(a, inversePolynomial(R));
		for(int i = 0; i < c; ++i){
			R[i] = R[i] + a[i];
			if(R[i] >= p) R[i] -= p;
			R[i] = (lli)R[i] * inv2 % p;
		}
	}
	R.resize(A.size());
	return R;
}

poly derivative(poly A){
	for(int i = 0; i < A.size(); ++i)
		A[i] = (lli)A[i] * i % p;
	if(!A.empty()) A.erase(A.begin());
	return A;
}

poly integral(poly A){
	for(int i = 0; i < A.size(); ++i)
		A[i] = (lli)A[i] * (powerMod(i+1, p-2, p)) % p;
	A.insert(A.begin(), 0);
	return A;
}

poly logarithm(poly A){
	assert(A[0] == 1);
	int n = A.size();
	A = convolution<p, g>(derivative(A), inversePolynomial(A));
	A.resize(n);
	A = integral(A);
	A.resize(n);
	return A;
}

poly exponential(const poly & A){
	assert(A[0] == 0);
	//E(x) = E(x)(1-ln(E(x))+A(x))
	poly E(1, 1);
	while(E.size() < A.size()){
		size_t c = 2*E.size();
		E.resize(c);
		poly S = logarithm(E);
		for(int i = 0; i < c && i < A.size(); ++i){
			S[i] = A[i] - S[i];
			if(S[i] < 0) S[i] += p;
		}
		S[0] = 1;
		E = convolution<p, g>(E, S);
		E.resize(c);
	}
	E.resize(A.size());
	return E;
}

//returns Q(x), where A(x)=B(x)Q(x)+R(x)
poly quotient(poly A, poly B){
	int n = A.size(), m = B.size();
	if(n < m) return poly{};
	reverse(A.begin(), A.end());
	reverse(B.begin(), B.end());
	A.resize(n-m+1), B.resize(n-m+1);
	A = convolution<p, g>(A, inversePolynomial(B));
	A.resize(n-m+1);
	reverse(A.begin(), A.end());
	return A;
}

//returns R(x), where A(x)=B(x)Q(x)+R(x)
poly remainder(poly A, const poly & B){
	int n = A.size(), m = B.size();
	if(n >= m){
		poly C = convolution<p, g>(quotient(A, B), B);
		A.resize(m-1);
		for(int i = 0; i < m-1; ++i){
			A[i] -= C[i];
			if(A[i] < 0) A[i] += p;
		}
	}
	return A;
}

//evaluates all the points in P(x)
vector<int> multiEvaluate(const poly & P, const vector<int> & points){
	int n = points.size();
	vector<poly> t(n<<1), r(n<<1); vector<vector<int>> e(n<<1);
	vector<bool> calc(n<<1);
	vector<int> ans(n);
	for(int i = 0; i < n; ++i){
		t[n+i] = {(p - points[i]) % p, 1};
		e[n+i].push_back(i);
	}
	for(int i = n-1; i > 0; --i){
		t[i] = convolution<p, g>(t[i<<1], t[i<<1|1]);
		e[i] = e[i<<1];
		e[i].insert(e[i].end(), e[i<<1|1].begin(), e[i<<1|1].end());
	}
	auto naive = [&](const poly& P, int x){
		int y = 0;
		for(int i = (int)P.size()-1; i >= 0; --i){
			y = ((lli)y*x + P[i]) % p;
		}
		return y;
	};
	r[1] = remainder(P, t[1]);
	for(int i = 1; i < n; ++i){
		if(calc[i]){
			calc[i<<1] = calc[i<<1|1] = true;
		}else if(e[i].size() < 400){
			for(int pos : e[i]){
				r[n+pos] = {naive(r[i], points[pos])};
			}
			calc[i<<1] = calc[i<<1|1] = true;
		}else{
			r[i<<1] = remainder(r[i], t[i<<1]);
			r[i<<1|1] = remainder(r[i], t[i<<1|1]);
		}
	}
	for(int i = 0; i < n; ++i){
		ans[i] = r[n+i][0];
	}
	return ans;
}

//finds a polynomial P(x) such that P(x[i]) = y[i]
poly interpolate(const vector<int>& x, const vector<int>& y){
	int n = x.size();
	vector<poly> t(n<<1), r(n<<1);
	for(int i = 0; i < n; ++i){
		t[n+i] = {(p - x[i]) % p, 1};
	}
	for(int i = n-1; i > 0; --i){
		t[i] = convolution<p, g>(t[i<<1], t[i<<1|1]);
	}
	vector<int> Q = multiEvaluate(derivative(t[1]), x);
	for(int i = 0; i < n; ++i){
		r[n+i] = {y[i] * powerMod(Q[i], p-2, p) % p};
	}
	for(int i = n-1; i > 0; --i){
		r[i] = convolution<p, g>(r[i<<1], t[i<<1|1]);
		poly rhs = convolution<p, g>(r[i<<1|1], t[i<<1]);
		r[i].resize(max(r[i].size(), rhs.size()));
		for(int j = 0; j < rhs.size(); ++j){
			r[i][j] += rhs[j];
			if(r[i][j] >= p) r[i][j] -= p;
		}
	}
	return r[1];
}

void clean(poly& A){
	while(!A.empty() && A.back() == 0) A.pop_back();
}

poly operator+(const poly& a, const poly& b){
	poly c(max(a.size(), b.size()));
	for(int i = 0; i < c.size(); ++i){
		if(i < a.size()) c[i] = a[i];
		if(i < b.size()) c[i] += b[i];
		if(c[i] >= p) c[i] -= p;
	}
	clean(c);
	return c;
}

poly operator-(const poly& a, const poly& b){
	poly c(max(a.size(), b.size()));
	for(int i = 0; i < c.size(); ++i){
		if(i < a.size()) c[i] = a[i];
		if(i < b.size()) c[i] -= b[i];
		if(c[i] < 0) c[i] += p;
	}
	clean(c);
	return c;
}

const poly zero, one = {1};
poly operator*(const poly& a, const poly& b){
	if(a.empty() || b.empty()) return {};
	poly ans = convolution<p,g>(a, b);
	clean(ans);
	return ans;
}

using mat = array<poly, 4>;
using arr = array<poly, 2>;
mat operator*(const mat& A, const mat& B){
	return {A[0]*B[0] + A[1]*B[2], A[0]*B[1] + A[1]*B[3], A[2]*B[0] + A[3]*B[2], A[2]*B[1] + A[3]*B[3]};
}

arr operator*(const mat& A, const arr& b){
	return {A[0]*b[0] + A[1]*b[1], A[2]*b[0] + A[3]*b[1]};
}

mat pgcd(arr a){
	assert(a[0].size() > a[1].size() && !a[1].empty());
	int m = a[0].size()/2;
	if(a[1].size() <= m) return {one, zero, zero, one};
	auto R = pgcd({poly(a[0].begin() + m, a[0].end()), poly(a[1].begin() + m, a[1].end())});
	a = R*a;
	if(a[1].size() <= m) return R;
	mat Q = {zero, one, one, zero - quotient(a[0], a[1])};
	R = Q*R, a = Q*a;
	if(a[1].size() <= m) return R;
	int k = 2*m + 1 - a[0].size();
	return pgcd({poly(a[0].begin() + k, a[0].end()), poly(a[1].begin() + k, a[1].end())}) * R;
}

mat egcd(arr a){
	assert(a[0].size() > a[1].size() && !a[1].empty());
	auto m0 = pgcd(a);
	a = m0*a;
	if(a[1].empty()) return m0;
	mat Q = {zero, one, one, zero - quotient(a[0], a[1])};
	m0 = Q*m0, a = Q*a;
	if(a[1].empty()) return m0;
	return egcd(a) * m0;
}

array<poly, 3> extgcd(const poly& a, const poly& b){
	mat Q = {zero, one, one, zero - quotient(a, b)};
	auto m = Q;
	auto ap = Q*arr{a, b};
	if(!ap[1].empty()) m = egcd(ap) * m;
	return {a*m[0] + b*m[1], m[0], m[1]};
}

//it evaluates 1, w, w^2, ..., w^(n-1) on the polynomial a(x)
//in this example we do a DFT with arbitrary size
vector<comp> bluestein(vector<comp> A){
	int n = A.size(), m = nearestPowerOfTwo(2*n-1);
	comp w = polar(1.0, 2*PI/n), w1 = 1, w2 = 1;
	vector<comp> p(m), q(m), b(n);
	for(int k = 0; k < n; ++k, w2 *= w1, w1 *= w){
		b[k] = w2;
		p[n-1-k] = A[k] / b[k];
		q[k] = b[k];
		if((n&1) == 1 && k < n-1) q[k+n] = q[k];
		else if((n&1) == 0 && k < n-1) q[k+n] = -q[k]; // q[k]*w^(n/2)
	}
	fft(p, 1), fft(q, 1);
	for(int i = 0; i < m; i++)
		p[i] *= q[i];
	fft(p, -1);
	for(int k = 0; k < n; ++k)
		A[k] = p[k+n-1] / b[k];
	return A;
}

//A and B are real-valued vectors, just 2 fft's instead of 3
vector<double> convolutionTrick(const vector<double> & A, const vector<double> & B){
	int sz = A.size() + B.size() - 1;
	int size = nearestPowerOfTwo(sz);
	vector<comp> C(size);
	comp I(0, 1);
	for(int i = 0; i < A.size() || i < B.size(); ++i){
		if(i < A.size()) C[i] += A[i];
		if(i < B.size()) C[i] += I*B[i];
	}
	fft(C, 1);
	vector<comp> D(size);
	for(int i = 0, j = 0; i < size; ++i){
		j = (size-1) & (size-i);
		D[i] = (conj(C[j]*C[j]) - C[i]*C[i]) * 0.25 * I;
	}
	fft(D, -1);
	vector<double> E;
	for_each(D.begin(), D.begin() + sz, [&](comp x){E.push_back(x.real());});
	return E;
}

//convolution with arbitrary modulo using only 4 fft's
poly convolutionMod(const poly & A, const poly & B, int mod){
	int s = sqrt(mod);
	int sz = A.size() + B.size() - 1;
	int size = nearestPowerOfTwo(sz);
	vector<comp> a(size), b(size);
	for(int i = 0; i < A.size(); ++i)
		a[i] = comp(A[i] % s, A[i] / s);
	for(int i = 0; i < B.size(); ++i)
		b[i] = comp(B[i] % s, B[i] / s);
	fft(a, 1), fft(b, 1);
	comp I(0, 1);
	vector<comp> c(size), d(size);
	for(int i = 0, j = 0; i < size; ++i){
		j = (size-1) & (size-i);
		comp e = (a[i] + conj(a[j])) * 0.5;
		comp f = (conj(a[j]) - a[i]) * 0.5 * I;
		comp g = (b[i] + conj(b[j])) * 0.5;
		comp h = (conj(b[j]) - b[i]) * 0.5 * I;
		c[i] = e * g + I * (e * h + f * g);
		d[i] = f * h;
	}
	fft(c, -1), fft(d, -1);
	poly D(sz);
	for(int i = 0, j = 0; i < sz; ++i){
		j = (size-1) & (size-i);
		int p0 = (lli)round(real(c[i])) % mod;
		int p1 = (lli)round(imag(c[i])) % mod;
		int p2 = (lli)round(real(d[i])) % mod;
		D[i] = p0 + s*(p1 + (lli)p2*s % mod) % mod;
		if(D[i] >= mod) D[i] -= mod;
		if(D[i] < 0) D[i] += mod;
	}
	return D;
}

//convolution with arbitrary modulo using CRT
//slower but with no precision errors
const int a = 998244353, b = 985661441, c = 754974721;
const lli a_b = powerMod(a, b-2, b), a_c = powerMod(a, c-2, c), b_c = powerMod(b, c-2, c);
poly convolutionModCRT(const poly & A, const poly & B, int mod){
	poly P = convolution<a, 3>(A, B);
	poly Q = convolution<b, 3>(A, B);
	poly R = convolution<c, 11>(A, B);
	poly D(P.size());
	for(int i = 0; i < D.size(); ++i){
		int x1 = P[i] % a;
		if(x1 < 0) x1 += a;
		int x2 = a_b * (Q[i] - x1) % b;
		if(x2 < 0) x2 += b;
		int x3 = (a_c * (R[i] - x1) % c - x2) * b_c % c;
		if(x3 < 0) x3 += c;
		D[i] = x1 % mod + a*(x2 + (lli)x3*b % mod) % mod;
		if(D[i] >= mod) D[i] -= mod;
		if(D[i] < 0) D[i] += mod;
	}
	return D;
}

//Fast Walsh-Hadamard transform, works with any modulo p
//op: 0(OR), 1(AND), 2(XOR), A.size() must be power of 2
void fwt(vector<int> & A, int op, int inv){
	int n = A.size();
	for(int k = 1; k < n; k <<= 1)
		for(int i = 0; i < n; i += k << 1)
			for(int j = 0; j < k; ++j){
				int u = A[i + j], v = A[i + j + k];
				int sum = u + v < p ? u + v : u + v - p;
				int rest = u - v < 0 ? u - v + p : u - v;
				if(inv == -1){
					if(op == 0)
						A[i + j + k] = rest ? p - rest : 0;
					else if(op == 1)
						A[i + j] = rest;
					else if(op == 2)
						A[i + j] = sum, A[i + j + k] = rest;
				}else{
					if(op == 0)
						A[i + j + k] = sum;
					else if(op == 1)
						A[i + j] = sum;
					else if(op == 2)
						A[i + j] = sum, A[i + j + k] = rest;
				}
			}
	if(inv == -1 && op == 2){
		lli nrev = powerMod(n, p-2, p);
		for(int i = 0; i < n; ++i)
			A[i] = A[i] * nrev % p;
	}
}

mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());

int aleatorio_int(int a, int b){
	uniform_int_distribution<int> dist(a, b);
	return dist(rng);
}

double aleatorio_double(double a, double b){
	uniform_real_distribution<double> dist(a, b);
	return dist(rng);
}

void test_fft(){
	int sz = 1<<20;
	vector<comp> A(sz);
	for(int i = 0; i < sz; ++i){
		A[i] = comp(aleatorio_double(0, 1), aleatorio_double(0, 1));
	}
	clock_t start = clock();
	fft(A, 1);
	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << duration << "\n";
}

void test_ntt(){
	int sz = 1<<20;
	poly A(sz);
	for(int i = 0; i < sz; ++i){
		A[i] = aleatorio_int(0, 9);
	}
	clock_t start = clock();
	ntt<p,g>(A, 1);
	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << duration << "\n";
}

void test_random_conv_fft(){
	int deg = 1e6;
	vector<comp> A(deg + 1), B(deg + 1);
	for(int i = 0; i <= deg; i++){
		A[i] = rand() % 10;
		B[i] = rand() % 10;
	}
	clock_t start = clock();
	A = convolution(A, B);
	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << duration << "\n";
}

void test_random_conv_ntt(){
	int deg = 1e6;
	poly A(deg + 1), B(deg + 1);
	for(int i = 0; i <= deg; i++){
		A[i] = rand() % 10;
		B[i] = rand() % 10;
	}
	clock_t start = clock();
	A = convolution<p, g>(A, B);
	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << duration << "\n";
}

void test_random_mult(){
	int digits = 1e6;
	stringstream ss1, ss2;
	for(int i = 1; i <= digits; i++){
		ss1 << rand() % 10;
		ss2 << rand() % 10;
	}
	clock_t start = clock();
	string res = multiplyNumbers(ss1.str(), ss2.str());
	double duration = (clock() - start) / (double) CLOCKS_PER_SEC;
	cout << duration << "\n";
}

int main(){
	srand(time(NULL));
	/*string a, b;
	cin >> a >> b;
	cout << multiplyNumbers(a, b) << "\n";*/
	test_fft();
	test_ntt();
	//test_random_mult();
	//test_random_conv_fft();
	//test_random_conv_ntt();
	//test_fft();
	//test_ntt();

	/*int m; lli n;
	cin >> m >> n;
	vector<int> fact(m + 2, 1), invfact(m + 2, 1), den(m + 1, 1);
	for(int i = 1; i <= m + 1; ++i){
		fact[i] = (lli)fact[i - 1] * i % p;
		invfact[i] = powerMod(fact[i], p-2, p);
	}
	for(int i = 0; i <= m; ++i){
		den[i] = invfact[i + 1];
		if(i & 1) den[i] = (p - den[i]) % p;
	}
	vector<int> bernoulli = inversePolynomial(den);
	for(int i = 0; i <= m; ++i){
		bernoulli[i] = (lli)bernoulli[i] * fact[i] % p;
	}
	int sum = 0;
	for(int k = 0; k <= m; ++k){
		sum += (lli)fact[m + 1] * invfact[k] % p * invfact[m + 1 - k] % p * bernoulli[k] % p * powerMod(n % p, m + 1 - k, p) % p;
		if(sum >= p) sum -= p;
	}
	sum = (lli)sum * powerMod(m + 1, p-2, p) % p;
	cout << sum << "\n";*/

	/*int n;
	cin >> n;
	vector<int> den(n + 1);
	den[0] = 1; den[1] = p - 4;
	den = sqrtPolynomial(den);
	den[0]++;
	vector<int> catalan = inversePolynomial(den);
	for(int i = 0; i <= n; ++i){
		catalan[i] = (lli)2 * catalan[i] % p;
		cout << catalan[i] << " ";
	}*/

	/*vector<int> A = {569675680, 478964123, 346798452, 146739485, 649785142}, B = {126741258, 700174685, 115649658};
	A = convolutionModCRT(A, B, 1e9+7);
	for(auto c : A) cout << c << " ";*/

	/*vector<comp> test = {comp(5,-3), comp(2,1), comp(0,7), comp(-4,9), comp(8,0)};
	test = bluestein(bluestein(test));
	for(auto t : test) cout << t/comp(test.size()) << " "; cout << "\n";*/

	/*vector<int> A = {9, 7, 2, 11, 3, 4, 5, 1}, B = {3, 7, 2, 5, 1};
	auto Q = quotient(A, B);
	for(auto q : Q) cout << q << " "; cout << "\n";
	auto R = remainder(A, B);
	for(auto r : R) cout << r << " ";*/

	/*int n = 60000;
	set<int> pts;
	while(pts.size() < n){
		pts.insert(aleatorio_int(0, p-1));
	}
	vector<int> P(n), points(pts.begin(), pts.end());
	for(int i = 0; i < n; ++i){
		P[i] = aleatorio_int(0, p-1);
	}
	vector<int> evals = multiEvaluate(P, points);
	vector<int> naive(n);
	for(int i = 0; i < n; ++i){
		for(int j = n-1; j >= 0; --j){
			naive[i] = ((lli)naive[i] * points[i] % p + P[j]) % p;
		}
	}
	assert(naive == evals);
	vector<int> Q = interpolate(points, evals);
	assert(P == Q);*/

	/*int M = 100;
	vector<int> fact(M+1, 1), invfact(M+1, 1);
	for(int i = 1; i <= M; ++i){
		fact[i] = (lli)i * fact[i-1] % p;
		invfact[i] = powerMod(fact[i], p-2, p);
	}*/
	
	/*vector<int> C(M+1, 1);
	for(int i = 2; i <= M; ++i){
		C[i] = powerMod(2, (lli)i * (i-1)/2, p) * invfact[i] % p;
	}
	vector<int> D = logarithm(C);
	D[0]++;
	for(int i = 1; i <= M; ++i){
		D[i] = (lli)D[i] * fact[i] % p;
	}
	for(int i = 1; i <= M; ++i){
		cout << i << " " << D[i] << "\n";
	}*/

	/*vector<int> B(M+1);
	B[1] = 1;
	B = exponential(B);
	B[0]--;
	B = exponential(B);
	for(int i = 1; i <= M; ++i){
		B[i] = (lli)B[i] * fact[i] % p;
	}
	for(int i = 1; i <= M; ++i){
		cout << i << " " << B[i] << "\n";
	}*/
	return 0;
}