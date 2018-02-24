#include <bits/stdc++.h>
using namespace std;
typedef complex<double> comp;
double PI = acos(-1.0);

int nearestPowerOfTwo(int n){
	int ans = 1;
	while(ans < n) ans <<= 1;
	return ans;
}

bool isZero(comp z){
	return abs(z.real()) < 1e-3;
}

void swapPositions(vector<comp> & X){
	int n = X.size();
	int j = n >> 1;
	for(int i = 1; i < n - 1; i++){
		if(i < j) swap(X[i], X[j]);
		int k = n >> 1;
		while(j >= k){
			j -= k;
			k >>= 1;
		}
		if(j < k){
			j += k;
		}
	}
}

void fft(vector<comp> & X, int inv){
	int n = X.size();
	swapPositions(X);
	comp w1 = polar(1.0, 2.0 * PI * inv / n);
	vector<comp> w(n);
	w[0] = 1;
	for(int i = 1; i < n; i++){
		w[i] = w[i - 1] * w1;
	}
	int pot = n >> 1;
	for(int i = 1; i < n; i <<= 1){
		for(int j = 0; j < i; j++){
			for(int k = 0; k < pot; k++){
				int first = j + 2 * i * k, second = first + i;
				comp r = w[pot * j] * X[second];
				X[second] = X[first] - r;
				X[first] += r;
			}
		}
		pot >>= 1;
	}
	if(inv == -1){
		for(int i = 0; i < n; i++){
			X[i] /= n;
		}
	}
}

void quitar(vector<comp> & X){
	while(isZero(X.back())) X.pop_back();
	if(X.size() == 0) X.push_back(0);
}

void multiplyPolynomials(vector<comp> & A, vector<comp> & B){
	int degree = nearestPowerOfTwo(A.size() + B.size() - 1);
	A.resize(degree);
	B.resize(degree);
	fft(A, 1);
	fft(B, 1);
	for(int i = 0; i < degree; i++){
		A[i] *= B[i];
	}
	fft(A, -1);
	quitar(A);
}

int main(){
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

	cout << "\n" << duration;
	return 0;
}