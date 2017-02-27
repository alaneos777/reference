#include <iostream>
#include <vector>
#include <sstream>

using namespace std;

typedef long long int ull;

ull gcd(ull a, ull b){
    ull r;
    while(b!=0) r = a%b, a = b, b = r;
    return a;
}

ull mod(ull a, ull b){
    ull ans = a % b;
    if(ans < 0) ans += b;
    return ans;
}

struct enteroModular{
    ull a, n;
    enteroModular(ull x, ull y){
        a = mod(x, y), n = y;
    }
    enteroModular(ull x){
        a = mod(x, 2), n = 2;
    }
    enteroModular(){
        a = 0, n = 2;
    }
    enteroModular operator+(const enteroModular & e) const{
        return enteroModular(a + e.a, n);
    }
    enteroModular operator-() const{
        return enteroModular(-a, n);
    }
    enteroModular operator-(const enteroModular & e) const{
        return *this + (-e);
    }
    enteroModular operator*(const enteroModular & e) const{
        return enteroModular(a * e.a, n);
    }
    enteroModular inverso() const{
        ull q, r0 = a, r1 = n, ri, s0 = 1, s1 = 0, si;
        while(r1 != 0){
            q = r0 / r1, ri = r0 % r1;
            si = s0 - s1 * q;
            r0 = r1, r1 = ri;
            s0 = s1, s1 = si;
        }
        return enteroModular(s0, n);
    }
    enteroModular operator/(const enteroModular & e) const{
        return enteroModular(a * e.inverso().a, n);
    }
    enteroModular operator+=(const enteroModular & e){
        *this = *this + e;
        return *this;
    }
    enteroModular operator-=(const enteroModular & e){
        *this = *this - e;
        return *this;
    }
    enteroModular operator*=(const enteroModular & e){
        *this = *this * e;
        return *this;
    }
    enteroModular operator/=(const enteroModular & e){
        *this = *this / e;
        return *this;
    }
    enteroModular operator++(int xd){
        *this = *this + enteroModular(1, n);
        return * this;
    }
    enteroModular operator--(int xd){
        *this = *this - enteroModular(1, n);
        return * this;
    }
    bool operator==(const enteroModular & e) const{
        return mod(a, n) == mod(e.a, n);
    }
    bool operator!=(const enteroModular & e) const{
        return mod(a, n) != mod(e.a, n);
    }
};

struct fraccion{
    ull num, den;
    fraccion(){
        num = 0, den = 1;
    }
    fraccion(ull x, ull y){
        if(y < 0){
            x *= -1, y *=-1;
        }
        ull d = gcd(abs(x), abs(y));
        num = x/d, den = y/d;
    }
    fraccion(ull v){
        num = v;
        den = 1;
    }
    fraccion operator+(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return fraccion(num*(f.den/d) + f.num*(den/d), den*(f.den/d));
    }
    fraccion operator-() const{
        return fraccion(-num, den);
    }
    fraccion operator-(const fraccion& f) const{
        return *this + (-f);
    }
    fraccion operator*(const fraccion& f) const{
        return fraccion(num*f.num, den*f.den);
    }
    fraccion operator/(const fraccion& f) const{
        return fraccion(num*f.den, den*f.num);
    }
    fraccion operator+=(const fraccion& f){
        *this = *this + f;
        return *this;
    }
    fraccion operator-=(const fraccion& f){
        *this = *this - f;
        return *this;
    }
    fraccion operator++(int xd){
        *this = *this + 1;
        return *this;
    }
    fraccion operator--(int xd){
        *this = *this - 1;
        return *this;
    }
    fraccion operator*=(const fraccion& f){
        *this = *this * f;
        return *this;
    }
    fraccion operator/=(const fraccion& f){
        *this = *this / f;
        return *this;
    }
    bool operator==(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num*(f.den/d) == (den/d)*f.num);
    }
    bool operator!=(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num*(f.den/d) != (den/d)*f.num);
    }
    bool operator >(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num*(f.den/d) > (den/d)*f.num);
    }
    bool operator <(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num*(f.den/d) < (den/d)*f.num);
    }
    bool operator >=(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num*(f.den/d) >= (den/d)*f.num);
    }
    bool operator <=(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num*(f.den/d) <= (den/d)*f.num);
    }
    fraccion inverso() const{
        return fraccion(den, num);
    }
    fraccion fabs() const{
        fraccion nueva;
        nueva.num = abs(num);
        nueva.den = den;
        return nueva;
    }
    string str() const{
        stringstream ss;
        ss << num;
        if(den != 1) ss << "/" << den;
        return ss.str();
    }
};

ostream &operator<<(ostream &os, const enteroModular & e) { 
    return os << e.a;
}

ostream &operator<<(ostream &os, const fraccion & f) { 
    return os << f.str();
}

istream &operator>>(istream &is, fraccion & f){
    ull num = 0, den = 1;
    string str;
    is >> str;
    size_t pos = str.find("/");
    if(pos == string::npos){
        istringstream(str) >> num;
    }else{
        istringstream(str.substr(0, pos)) >> num;
        istringstream(str.substr(pos + 1)) >> den;
    }
    fraccion nueva(num, den);
    f = nueva;
    return is;
}

template <typename entrada>
using fila = vector<entrada>;

template <typename entrada>
using matrix = vector< fila<entrada> >;

template <typename entrada>
struct sistema{
	matrix<entrada> A;
	fila<entrada> b;
	int m, n;

	sistema(int _m, int _n){
		m = _m, n = _n;
		A.resize(m, fila<entrada>(n));
		b.resize(m);
	}

	void multiplicarFilaPorEscalar(int k, entrada c){
		for(int j = 0; j < n; j++){
			A[k][j] *= c;
		}
		b[k] *= c;
	}

	void intercambiarFilas(int k, int l){
		swap(A[k], A[l]);
		swap(b[k], b[l]);
	}

	void sumaMultiploFilaAOtra(int k, int l, entrada c){
		for(int j = 0; j < n; j++){
			A[k][j] += c * A[l][j];
		}
		b[k] += c * b[l];
	}

    void imprime_matriz(){
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                cout << A[i][j] << " ";
            }
            cout << b[i] << "\n";
        }
        cout << "\n";
    }

	void gauss_jordan(){
        cout << "Matriz original:\n";
        imprime_matriz();
		int i = 0, j = 0;
		while(i < m && j < n){
			if(A[i][j] == 0){
				for(int f = i + 1; f < m; f++){
					if(A[f][j] != 0){
						intercambiarFilas(f, i);
                        cout << "F_" << (f + 1) << " <-> F_" << (i + 1) << ":\n";
                        imprime_matriz();
						break;
					}
				}
			}
			if(A[i][j] != 0){
                entrada inv_mult = A[i][j].inverso();
				multiplicarFilaPorEscalar(i, inv_mult);
                cout << "(" << inv_mult << ")F_" << (i + 1) << " -> F_" << (i + 1) << ":\n";
                imprime_matriz();
				for(int f = 0; f < m; f++){
                    if(f != i){
                        entrada inv_adit = -A[f][j];
                        sumaMultiploFilaAOtra(f, i, inv_adit);
                        cout << "F_" << (f + 1) << " + (" << inv_adit << ")F_" << (i + 1) << " -> F_" << (f + 1) << ":\n";
                        imprime_matriz();
                    }
                }
				i++;
			}
			j++;
		}
	}
};

void pedirValores(sistema<fraccion> & S){
    for(int i = 0; i < S.m; i++){
        cout << "Introduce los coeficientes y t\202rmino independiente de la ecuaci\242n " << (i + 1) << ": ";
        for(int j = 0; j < S.n; j++){
            cin >> S.A[i][j];
        }
        cin >> S.b[i];
    }
    S.gauss_jordan();
}

void pedirValores(sistema<enteroModular> & S, ull p){
    ull valor;
    for(int i = 0; i < S.m; i++){
        cout << "Introduce los coeficientes y t\202rmino independiente de la ecuaci\242n " << (i + 1) << ": ";
        for(int j = 0; j < S.n; j++){
            cin >> valor;
            S.A[i][j] = enteroModular(valor, p);
        }
        cin >> valor;
        S.b[i] = enteroModular(valor, p);
    }
    cout << "\n";
    S.gauss_jordan();
}

int main()
{
    int m, n;
    ull p;
    string campo;
    cout << "Introduce el n\243mero de ecuaciones: ";
    cin >> m;
    cout << "Introduce el n\243mero de inc\242gnitas: ";
    cin >> n;
    cout << "Introduce Q para trabajar en los racionales, o un numero primo p para trabajar en F_p: ";
    cin >> campo;
    if(campo == "Q"){
        sistema<fraccion> S(m, n);
        pedirValores(S);
    }else{
        istringstream(campo) >> p;
        sistema<enteroModular> S(m, n);
        pedirValores(S, p);
    }
    return 0;
}