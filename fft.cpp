#include <bits/stdc++.h>
using namespace std;
using lli = long long int;
using comp = complex<double>;
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
			wp[j] = polar(1.0, PI * j / k * inv);
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

int inverse(int a, int n){
	int r0 = a, r1 = n, ri, s0 = 1, s1 = 0, si;
	while(r1){
		si = s0 - s1 * (r0 / r1), s0 = s1, s1 = si;
		ri = r0 % r1, r0 = r1, r1 = ri;
	}
	if(s0 < 0) s0 += n;
	return s0;
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

template<int prime, int gen>
void ntt(vector<int> & X, int inv){
	int n = X.size();
	for(int i = 1, j = 0; i < n - 1; ++i){
		for(int k = n >> 1; (j ^= k) < k; k >>= 1);
		if(i < j) swap(X[i], X[j]);
	}
	vector<lli> wp(n>>1, 1);
	for(int k = 1; k < n; k <<= 1){
		lli wk = powerMod(gen, inv * (prime - 1) / (k<<1), prime);
		for(int j = 1; j < k; ++j)
			wp[j] = wp[j - 1] * wk % prime;
		for(int i = 0; i < n; i += k << 1){
			for(int j = 0; j < k; ++j){
				int u = X[i + j], v = X[i + j + k] * wp[j] % prime;
				X[i + j] = u + v < prime ? u + v : u + v - prime;
				X[i + j + k] = u - v < 0 ? u - v + prime : u - v;
			}
		}
	}
	if(inv == -1){
		lli nrev = inverse(n, prime);
		for(int i = 0; i < n; ++i)
			X[i] = X[i] * nrev % prime;
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

template<int prime, int gen>
vector<int> convolution(vector<int> A, vector<int> B){
	int sz = A.size() + B.size() - 1;
	int size = nearestPowerOfTwo(sz);
	A.resize(size), B.resize(size);
	ntt<prime, gen>(A, 1), ntt<prime, gen>(B, 1);
	for(int i = 0; i < size; i++)
		A[i] = (lli)A[i] * B[i] % prime;
	ntt<prime, gen>(A, -1);
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
	vector<int> X(a.size() - pos1), Y(b.size() - pos2);
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

vector<int> inversePolynomial(const vector<int> & A){
	vector<int> R(1, inverse(A[0], p));
	//R(x) = 2R(x)-A(x)R(x)^2
	while(R.size() < A.size()){
		int c = 2 * R.size();
		R.resize(c);
		vector<int> TR = R;
		TR.resize(2 * c);
		vector<int> TF(TR.size());
		for(int i = 0; i < c && i < A.size(); ++i)
			TF[i] = A[i];
		ntt<p, g>(TR, 1);
		ntt<p, g>(TF, 1);
		for(int i = 0; i < TR.size(); ++i)
			TR[i] = (lli)TR[i] * TR[i] % p * TF[i] % p;
		ntt<p, g>(TR, -1);
		for(int i = 0; i < c; ++i){
			R[i] = R[i] + R[i] - TR[i];
			if(R[i] < 0) R[i] += p;
			if(R[i] >= p) R[i] -= p;
		}
	}
	R.resize(A.size());
	return R;
}

const int inv2 = inverse(2, p);

vector<int> sqrtPolynomial(const vector<int> & A){
	int r0 = 1; //verify that r0^2 = A[0] mod p
	vector<int> R(1, r0);
	//R(x) = R(x)/2 + A(x)/(2R(x))
	while(R.size() < A.size()){
		int c = 2 * R.size();
		R.resize(c);
		vector<int> TF(c);
		for(int i = 0; i < c && i < A.size(); ++i)
			TF[i] = A[i];
		vector<int> IR = inversePolynomial(R);
		TF = convolution<p, g>(TF, IR);
		for(int i = 0; i < c; ++i){
			R[i] = R[i] + TF[i];
			if(R[i] >= p) R[i] -= p;
			R[i] = (lli)R[i] * inv2 % p;
		}
	}
	R.resize(A.size());
	return R;
}

vector<int> derivative(vector<int> A){
	for(int i = 0; i < A.size(); ++i)
		A[i] = (lli)A[i] * i % p;
	if(!A.empty()) A.erase(A.begin());
	return A;
}

vector<int> integral(vector<int> A){
	for(int i = 0; i < A.size(); ++i)
		A[i] = (lli)A[i] * (inverse(i+1, p)) % p;
	A.insert(A.begin(), 0);
	return A;
}

vector<int> logarithm(vector<int> A){
	assert(A[0] == 1);
	int n = A.size();
	A = convolution<p, g>(derivative(A), inversePolynomial(A));
	A.resize(n);
	A = integral(A);
	A.resize(n);
	return A;
}

vector<int> exponential(const vector<int> & A){
	assert(A[0] == 0);
	//E(x) = E(x)(1-ln(E(x))+A(x))
	vector<int> E(1, 1);
	while(E.size() < A.size()){
		int c = 2*E.size();
		E.resize(c);
		vector<int> S = logarithm(E);
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
vector<int> quotient(vector<int> A, vector<int> B){
	int n = A.size(), m = B.size();
	if(n < m) return vector<int>{0};
	reverse(A.begin(), A.end());
	reverse(B.begin(), B.end());
	A.resize(n-m+1), B.resize(n-m+1);
	A = convolution<p, g>(A, inversePolynomial(B));
	A.resize(n-m+1);
	reverse(A.begin(), A.end());
	return A;
}

//returns R(x), where A(x)=B(x)Q(x)+R(x)
vector<int> remainder(vector<int> A, const vector<int> & B){
	int n = A.size(), m = B.size();
	if(n >= m){
		vector<int> C = convolution<p, g>(quotient(A, B), B);
		A.resize(m-1);
		for(int i = 0; i < m-1; ++i){
			A[i] -= C[i];
			if(A[i] < 0) A[i] += p;
		}
	}
	return A;
}

//evaluates all the points in P(x), both the size of P and points must be the same
vector<int> multiEvaluate(const vector<int> & P, const vector<int> & points){
	int n = points.size();
	vector<vector<int>> prod(2*n - 1);
	function<void(int, int, int)> pre = [&](int v, int l, int r){
		if(l == r) prod[v] = vector<int>{(p - points[l]) % p, 1};
		else{
			int y = (l + r) / 2;
			int z = v + (y - l + 1) * 2;
			pre(v + 1, l, y);
			pre(z, y + 1, r);
			prod[v] = convolution<p, g>(prod[v + 1], prod[z]);
		}
	};
	pre(0, 0, n - 1);

	function<int(const vector<int>&, int)> eval = [&](const vector<int> & poly, int x0){
		int ans = 0;
		for(int i = (int)poly.size()-1; i >= 0; --i){
			ans = (lli)ans * x0 % p + poly[i];
			if(ans >= p) ans -= p;
		}
		return ans;
	};

	vector<int> res(n);
	function<void(int, int, int, vector<int>)> evaluate = [&](int v, int l, int r, vector<int> poly){
		poly = remainder(poly, prod[v]);
		if(poly.size() < 400){
			for(int i = l; i <= r; ++i)
				res[i] = eval(poly, points[i]);
		}else{
			if(l == r)
				res[l] = poly[0];
			else{
				int y = (l + r) / 2;
				int z = v + (y - l + 1) * 2;
				evaluate(v + 1, l, y, poly);
				evaluate(z, y + 1, r, poly);
			}
		}
	};
	evaluate(0, 0, n - 1, P);
	return res;
}

//it evaluates 1, w^2, w^4, ..., w^(2n-2) on the polynomial a(x)
//in this example we do a DFT with arbitrary size
vector<comp> bluestein(vector<comp> A){
	int n = A.size();
	int m = nearestPowerOfTwo(2*n-1);
	comp w = polar(1.0, PI / n), w1 = w, w2 = 1;
	vector<comp> p(m), q(m), b(n);
	for(int k = 0; k < n; ++k, w2 *= w1, w1 *= w*w){
		b[k] = w2;
		p[k] = A[k] * b[k];
		q[k] = (comp)1 / b[k];
		if(k) q[m-k] = q[k];
	}
	fft(p, 1), fft(q, 1);
	for(int i = 0; i < m; i++)
		p[i] *= q[i];
	fft(p, -1);
	for(int k = 0; k < n; ++k)
		A[k] = b[k] * p[k];
	return A;
}

//A and B are real-valued vectors
//just do 2 fft's instead of 3
vector<comp> convolutionTrick(const vector<comp> & A, const vector<comp> & B){
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
	D.resize(sz);
	return D;
}

//convolution with arbitrary modulo using only 4 fft's
vector<int> convolutionMod(const vector<int> & A, const vector<int> & B, int mod){
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
	vector<int> D(sz);
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
const lli a_b = inverse(a, b), a_c = inverse(a, c), b_c = inverse(b, c);
vector<int> convolutionModCRT(const vector<int> & A, const vector<int> & B, int mod){
	vector<int> P = convolution<a, 3>(A, B);
	vector<int> Q = convolution<b, 3>(A, B);
	vector<int> R = convolution<c, 11>(A, B);
	vector<int> D(P.size());
	for(int i = 0; i < D.size(); ++i){
		int x1 = P[i] % a;
		if(x1 < 0) x1 += a;
		int x2 = a_b * (Q[i] - x1) % b;
		if(x2 < 0) x2 += b;
		int x3 = (a_c * (R[i] - x1) % c - x2) * b_c % c;
		if(x3 < 0) x3 += c;
		D[i] = x1 + a*(x2 + (lli)x3*b % mod) % mod;
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
		lli nrev = inverse(n, p);
		for(int i = 0; i < n; ++i)
			A[i] = A[i] * nrev % p;
	}
}

void test_fft(){
	int degX, degY;	
	cin >> degX >> degY;
	vector<comp> X(degX + 1), Y(degY + 1);

	for(int i = 0; i <= degX; i++) cin >> X[i];
	for(int i = 0; i <= degY; i++) cin >> Y[i];

	std::clock_t start;
	double duration;
	start = std::clock();

	X = convolution(X, Y);

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	for(int i = 0; i < X.size(); i++) cout << (int)round(X[i].real()) << " ";

	cout << "\n" << duration << "\n";
}

void test_ntt(){
	int degX, degY;
	cin >> degX >> degY;
	vector<int> X(degX + 1), Y(degY + 1);

	for(int i = 0; i <= degX; i++) cin >> X[i];
	for(int i = 0; i <= degY; i++) cin >> Y[i];

	std::clock_t start;
	double duration;
	start = std::clock();

	X = convolution<p, g>(X, Y);

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

	for(int i = 0; i < X.size(); i++) cout << X[i] << " ";

	cout << "\n" << duration << "\n";
}

void test_random_fft(){
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

void test_random_ntt(){
	int deg = 1e6;
	vector<int> A(deg + 1), B(deg + 1);
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
	//test_random_mult();
	//test_random_fft();
	//test_random_ntt();
	//test_fft();
	//test_ntt();

	/*int m; lli n;
	cin >> m >> n;
	vector<int> fact(m + 2, 1), invfact(m + 2, 1), den(m + 1, 1);
	for(int i = 1; i <= m + 1; ++i){
		fact[i] = (lli)fact[i - 1] * i % p;
		invfact[i] = inverse(fact[i], p);
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
		sum += (lli)fact[m + 1] * invfact[k] % p * invfact[m + 1 - k] % p * bernoulli[k] % p * powerMod(n, m + 1 - k, p) % p;
		if(sum >= p) sum -= p;
	}
	sum = (lli)sum * inverse(m + 1, p) % p;
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
	for(auto t : test) cout << t << " "; cout << "\n";*/

	/*vector<int> A = {9, 7, 2, 11, 3, 4, 5, 1}, B = {3, 7, 2, 5, 1};
	auto Q = quotient(A, B);
	for(auto q : Q) cout << q << " "; cout << "\n";
	auto R = remainder(A, B);
	for(auto r : R) cout << r << " ";*/

	/*int n = 60000;
	vector<int> P(n), points(n);
	for(int i = 0; i < n; ++i){
		P[i] = rand() % p;
		points[i] = rand() % p;
	}
	vector<int> evals = multiEvaluate(P, points);
	vector<int> naive(n);
	for(int i = 0; i < n; ++i){
		for(int j = n-1; j >= 0; --j){
			naive[i] = ((lli)naive[i] * points[i] % p + P[j]) % p;
		}
	}
	for(int i = 0; i < n; ++i){
		assert(naive[i] == evals[i]);
	}*/

	/*int M = 100;
	vector<int> fact(M+1, 1), invfact(M+1, 1);
	for(int i = 1; i <= M; ++i){
		fact[i] = (lli)i * fact[i-1] % p;
		invfact[i] = inverse(fact[i], p);
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