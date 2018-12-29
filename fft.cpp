#include <bits/stdc++.h>
using namespace std;
typedef complex<double> comp;
typedef long long int lli;
double PI = acos(-1.0);

int nearestPowerOfTwo(int n){
	int ans = 1;
	while(ans < n) ans <<= 1;
	return ans;
}

void fft(vector<comp> & X, int inv){
	int n = X.size();
	int len, len2, i, j, k;
	for(i = 1, j = 0; i < n - 1; ++i){
		for(k = n >> 1; (j ^= k) < k; k >>= 1);
		if(i < j) swap(X[i], X[j]);
	}
	double ang;
	comp t, u, v;
	vector<comp> wlen_pw(n >> 1);
	wlen_pw[0] = 1;
	for(len = 2; len <= n; len <<= 1){
		ang = inv == -1 ? -2 * PI / len : 2 * PI / len;
		len2 = len >> 1;
		comp wlen(cos(ang), sin(ang));
		for(i = 1; i < len2; ++i){
			wlen_pw[i] = wlen_pw[i - 1] * wlen;
		}
		for(i = 0; i < n; i += len){
			for(j = 0; j < len2; ++j){
				t = X[i + j + len2] * wlen_pw[j];
				X[i + j + len2] = X[i + j] - t;
				X[i + j] += t;
			}
		}
	}
	if(inv == -1){
		for(i = 0; i < n; ++i){
			X[i] /= n;
		}
	}
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

const int p = 7340033;
const int root = 5;
const int root_1 = inverse(root, p);
const int root_pw = 1 << 20;

void ntt(vector<int> & X, int inv){
	int n = X.size();
	int len, len2, wlen, i, j, k, u, v, w;
	for(i = 1, j = 0; i < n - 1; ++i){
		for(k = n >> 1; (j ^= k) < k; k >>= 1);
		if(i < j) swap(X[i], X[j]);
	}
	for(len = 2; len <= n; len <<= 1){
		len2 = len >> 1;
		wlen = (inv == -1) ? root_1 : root;
		for(i = len; i < root_pw; i <<= 1){
			wlen = (lli)wlen * wlen % p;
		}
		for(i = 0; i < n; i += len){
			w = 1;
			for(j = 0; j < len2; ++j){
				u = X[i + j], v = (lli)X[i + j + len2] * w % p;
				X[i + j] = u + v < p ? u + v : u + v - p;
				X[i + j + len2] = u - v < 0 ? u - v + p : u - v;
				w = (lli)w * wlen % p;
			}
		}
	}
	if(inv == -1){
		int nrev = inverse(n, p);
		for(i = 0; i < n; ++i){
			X[i] = (lli)X[i] * nrev % p;
		}
	}
}

void convolution(vector<comp> & A, vector<comp> & B){
	int degree = A.size() + B.size() - 2;
	int size = nearestPowerOfTwo(degree + 1);
	A.resize(size);
	B.resize(size);
	fft(A, 1);
	fft(B, 1);
	for(int i = 0; i < size; i++){
		A[i] *= B[i];
	}
	fft(A, -1);
	A.resize(degree + 1);
}

void convolution(vector<int> & A, vector<int> & B){
	int degree = A.size() + B.size() - 2;
	int size = nearestPowerOfTwo(degree + 1);
	A.resize(size);
	B.resize(size);
	ntt(A, 1);
	ntt(B, 1);
	for(int i = 0; i < size; i++){
		A[i] = (lli)A[i] * B[i] % p;
	}
	ntt(A, -1);
	A.resize(degree + 1);
}

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
	for(int i = pos1, j = X.size() - 1; i < a.size(); ++i){
		X[j--] = a[i] - '0';
	}
	for(int i = pos2, j = Y.size() - 1; i < b.size(); ++i){
		Y[j--] = b[i] - '0';
	}
	convolution(X, Y);
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
	for(int i = X.size() - 1; i >= 0; --i){
		ss << X[i];
	}
	return ss.str();
}

vector<int> inversePolynomial(vector<int> & A){
	vector<int> R(1, inverse(A[0], p));
	while(R.size() < A.size()){
		int c = 2 * R.size();
		R.resize(c);
		vector<int> TR = R;
		TR.resize(nearestPowerOfTwo(2 * c));
		vector<int> TF(TR.size());
		for(int i = 0; i < c; ++i){
			TF[i] = A[i];
		}
		ntt(TR, 1);
		ntt(TF, 1);
		for(int i = 0; i < TR.size(); ++i){
			TR[i] = (lli)TR[i] * TR[i] % p * TF[i] % p;
		}
		ntt(TR, -1);
		TR.resize(2 * c);
		for(int i = 0; i < c; ++i){
			R[i] = R[i] + R[i] - TR[i];
			while(R[i] < 0) R[i] += p;
			while(R[i] >= p) R[i] -= p;
		}
	}
	R.resize(A.size());
	return R;
}

const int inv2 = inverse(2, p);

vector<int> sqrtPolynomial(vector<int> & A){
	int r0 = 1; //r0^2 = A[0] mod p
	vector<int> R(1, r0);
	while(R.size() < A.size()){
		int c = 2 * R.size();
		R.resize(c);
		vector<int> TF(c);
		for(int i = 0; i < c; ++i){
			TF[i] = A[i];
		}
		vector<int> IR = inversePolynomial(R);
		convolution(TF, IR);
		for(int i = 0; i < c; ++i){
			R[i] = R[i] + TF[i];
			if(R[i] >= p) R[i] -= p;
			R[i] = (lli)R[i] * inv2 % p;
		}
	}
	R.resize(A.size());
	return R;
}

void bluestein(vector<comp> & x, int inv){
	int n = x.size();
	comp w = polar(1.0, PI * inv / n), w1 = w, w2 = 1;
	vector<comp> p(n), q(2*n-1), b(n);
	for(int k = 0; k < n; ++k, w2 *= w1, w1 *= w*w){
		b[k] = w2;
		p[k] = x[k] * b[k];
		q[n-1-k] = q[n-1+k] = (comp)1 / b[k];
	}
	convolution(p, q);
	for(int k = 0; k < n; ++k){
		if(inv == -1) x[k] = b[k] * p[n-1+k] / (comp)n;
		else x[k] = b[k] * p[n-1+k];
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

	convolution(X, Y);

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

	convolution(X, Y);

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
	convolution(A, B);
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
	convolution(A, B);
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
	cout << multiplyNumbers(a, b) << "\n";
	test_random_mult();
	test_random_fft();
	test_random_ntt();
	test_fft();
	test_ntt();*/
	int n;
	cin >> n;
	vector<int> test(n);
	for(int i = 0; i < n; ++i) cin >> test[i];
	vector<int> R = sqrtPolynomial(test);
	for(int i = 0; i < R.size(); ++i) cout << R[i] << " ";
	cout << "\n";
	vector<int> R_ = R;
	convolution(R, R_);
	for(int i = 0; i < R.size(); ++i) cout << R[i] << " ";
	return 0;
}