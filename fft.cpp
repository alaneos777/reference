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

bool isZero(comp z){
	return abs(z.real()) < 1e-3;
}

const int p = 7340033;
const int root = 5;
const int root_1 = 4404020;
const int root_pw = 1 << 20;

int inverse(int a, int n){
    int r0 = a, r1 = n, ri, s0 = 1, s1 = 0, si;
    while(r1){
        si = s0 - s1 * (r0 / r1), s0 = s1, s1 = si;
        ri = r0 % r1, r0 = r1, r1 = ri;
    }
    if(s0 < 0) s0 += n;
    return s0;
}

template<typename T>
void swapPositions(vector<T> & X){
	int n = X.size();
	int bit;
	for (int i = 1, j = 0; i < n; ++i) {
		bit = n >> 1;
		while(j & bit){
			j ^= bit;
			bit >>= 1;
		}
		j ^= bit;
		if (i < j){
			swap (X[i], X[j]);
		}
	}
}

void fft(vector<comp> & X, int inv){
	int n = X.size();
    swapPositions<comp>(X);
    int len, len2, i, j;
    double ang;
    comp t, u, v;
    vector<comp> wlen_pw(n >> 1);
    wlen_pw[0] = 1;
    for(len = 2; len <= n; len <<= 1) {
        ang = inv == -1 ? -2 * PI / len : 2 * PI / len;
        len2 = len >> 1;
        comp wlen(cos(ang), sin(ang));
        for(i = 1; i < len2; ++i){
            wlen_pw[i] = wlen_pw[i - 1] * wlen;
        }
        for(i = 0; i < n; i += len) {
            for(j = 0; j < len2; ++j) {
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

void ntt(vector<int> & X, int inv) {
	int n = X.size();
	swapPositions<int>(X);
	int len, len2, wlen, i, j, u, v, w;
	for (len = 2; len <= n; len <<= 1) {
		len2 = len >> 1;
		wlen = (inv == -1) ? root_1 : root;
		for (i = len; i < root_pw; i <<= 1){
			wlen = wlen * 1ll * wlen % p;
		}
		for (i = 0; i < n; i += len) {
			w = 1;
			for (j = 0; j < len2; ++j) {
				u = X[i + j], v = X[i + j + len2] * 1ll * w % p;
				X[i + j] = u + v < p ? u + v : u + v - p;
				X[i + j + len2] = u - v < 0 ? u - v + p : u - v;
				w = w * 1ll * wlen % p;
			}
		}
	}
	if (inv == -1) {
		int nrev = inverse(n, p);
		for (i = 0; i < n; ++i){
			X[i] = X[i] * 1ll * nrev % p;
		}
	}
}

void multiplyPolynomials(vector<comp> & A, vector<comp> & B){
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

void multiplyPolynomials(vector<int> & A, vector<int> & B){
	int degree = A.size() + B.size() - 2;
	int size = nearestPowerOfTwo(degree + 1);
	A.resize(size);
	B.resize(size);
	ntt(A, 1);
	ntt(B, 1);
	for(int i = 0; i < size; i++){
		A[i] = A[i] * 1ll * B[i] % p;
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
	multiplyPolynomials(X, Y);
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

void test_fft(){
	int degX, degY;
	cin >> degX >> degY;
	vector<comp> X(degX + 1), Y(degY + 1);

	for(int i = 0; i <= degX; i++) cin >> X[i];
	for(int i = 0; i <= degY; i++) cin >> Y[i];

	std::clock_t start;
    double duration;
    start = std::clock();

	multiplyPolynomials(X, Y);

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

	multiplyPolynomials(X, Y);

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
	multiplyPolynomials(A, B);
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
	multiplyPolynomials(A, B);
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
	string a, b;
	cin >> a >> b;
	cout << multiplyNumbers(a, b) << "\n";
	test_random_mult();
	test_random_fft();
	test_random_ntt();
	test_fft();
	test_ntt();
	return 0;
}