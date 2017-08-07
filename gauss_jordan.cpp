#include <iostream>
#include <cmath>
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
        a = x, n = 0;
    }
    enteroModular(){
        a = 0, n = 0;
    }
    enteroModular operator+(const enteroModular & e) const{
        return enteroModular(a + e.a, max(n, e.n));
    }
    enteroModular operator-() const{
        return enteroModular(-a, n);
    }
    enteroModular operator-(const enteroModular & e) const{
        return *this + (-e);
    }
    enteroModular operator*(const enteroModular & e) const{
        return enteroModular(a * e.a, max(n, e.n));
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
        return enteroModular(a * e.inverso().a, max(n, e.n));
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
        return mod(a, max(n, e.n)) == mod(e.a, max(n, e.n));
    }
    bool operator!=(const enteroModular & e) const{
        return mod(a, max(n, e.n)) != mod(e.a, max(n, e.n));
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
struct matrix{
    vector< vector<entrada> > A;
    int m, n;

    matrix(int _m, int _n){
        m = _m, n = _n;
        A.resize(m, vector<entrada>(n, 0));
    }

    vector<entrada> & operator[] (int i){
        return A[i];
    }

    void multiplicarFilaPorEscalar(int k, entrada c){
        for(int j = 0; j < n; j++) A[k][j] *= c;
    }

    void intercambiarFilas(int k, int l){
        swap(A[k], A[l]);
    }

    void sumaMultiploFilaAOtra(int k, int l, entrada c){
        for(int j = 0; j < n; j++) A[k][j] += c * A[l][j];
    }

    void gauss_jordan(matrix & dest){
        cout << "Matriz original:\n" << str() << dest;
        int i = 0, j = 0;
        while(i < m && j < n){
            if(A[i][j] == 0){
                for(int f = i + 1; f < m; f++){
                    if(A[f][j] != 0){
                        intercambiarFilas(f, i);
                        dest.intercambiarFilas(f, i);
                        cout << "F_" << (f + 1) << " <-> F_" << (i + 1) << ":\n" << str() << dest;
                        break;
                    }
                }
            }
            if(A[i][j] != 0){
                if(A[i][j] != 1){
                    entrada inv_mult = A[i][j].inverso();
                    multiplicarFilaPorEscalar(i, inv_mult);
                    dest.multiplicarFilaPorEscalar(i, inv_mult);
                    cout << "(" << inv_mult << ")F_" << (i + 1) << " -> F_" << (i + 1) << ":\n" << str() << dest;
                }
                for(int f = 0; f < m; f++){
                    if(f != i && A[f][j] != 0){
                        entrada inv_adit = -A[f][j];
                        sumaMultiploFilaAOtra(f, i, inv_adit);
                        dest.sumaMultiploFilaAOtra(f, i, inv_adit);
                        cout << "F_" << (f + 1) << " + (" << inv_adit << ")F_" << (i + 1) << " -> F_" << (f + 1) << ":\n" << str() << dest;
                    }
                }
                i++;
            }
            j++;
        }
    }

    void gauss_jordan(){
        matrix xd(m, n);
        gauss_jordan(xd);
    }

    entrada eliminacion_gaussiana(matrix & dest){
        //cout << "Matriz original:\n" << str() << dest;
        int i = 0, j = 0;
        entrada determinante = 1;
        while(i < m && j < n){
            if(A[i][j] == 0){
                for(int f = i + 1; f < m; f++){
                    if(A[f][j] != 0){
                        intercambiarFilas(f, i);
                        dest.intercambiarFilas(f, i);
                        determinante *= -1;
                        //cout << "F_" << (f + 1) << " <-> F_" << (i + 1) << ":\n" << str() << dest;
                        break;
                    }
                }
            }
            determinante *= A[i][j];
            if(A[i][j] != 0){
                if(A[i][j] != 1){
                    entrada inv_mult = A[i][j].inverso();
                    multiplicarFilaPorEscalar(i, inv_mult);
                    dest.multiplicarFilaPorEscalar(i, inv_mult);
                    //cout << "(" << inv_mult << ")F_" << (i + 1) << " -> F_" << (i + 1) << ":\n" << str() << dest;
                }
                for(int f = i + 1; f < m; f++){
                    if(A[f][j] != 0){
                        entrada inv_adit = -A[f][j];
                        sumaMultiploFilaAOtra(f, i, inv_adit);
                        dest.sumaMultiploFilaAOtra(f, i, inv_adit);
                        //cout << "F_" << (f + 1) << " + (" << inv_adit << ")F_" << (i + 1) << " -> F_" << (f + 1) << ":\n" << str() << dest;
                    }
                }
                i++;
            }
            j++;
        }
        return determinante;
    }

    entrada eliminacion_gaussiana(){
        matrix xd(m, n);
        return eliminacion_gaussiana(xd);
    }

    static entrada delta(int i, int j){
        if(i == j) return 1;
        else return 0;
    }

    static matrix identidad(int n){
        matrix<entrada> id(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                id[i][j] = delta(i, j);
            }
        }
        return id;
    }

    static matrix elemental_1(int n, int k, entrada c){
        matrix<entrada> E(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == k) E[i][j] = c * delta(k, j);
                else E[i][j] = delta(i, j);
            }
        }
        return E;
    }

    static matrix elemental_2(int n, int k, int l){
        matrix<entrada> E(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == k) E[i][j] = delta(l, j);
                else if(i == l) E[i][j] = delta(k, j);
                else E[i][j] = delta(i, j);
            }
        }
        return E;
    }

    static matrix elemental_3(int n, int k, int l, entrada c){
        matrix<entrada> E(n, n);
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == k) E[i][j] = delta(k, j) + c * delta(l, j);
                else E[i][j] = delta(i, j);
            }
        }
        return E;
    }

    matrix operator+(const matrix & B) const{
        if(m == B.m && n == B.n){
            matrix<entrada> C(m, n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++){
                    C[i][j] = A[i][j] + B.A[i][j];
                }
            }
            return C;
        }else{
            return *this;
        }
    }

    matrix operator+=(const matrix & M){
        *this = *this + M;
        return *this;
    }

    matrix operator-() const{
        matrix<entrada> C(m, n);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                C[i][j] = -A[i][j];
            }
        }
        return C;
    }

    matrix operator-(const matrix & B) const{
        return *this + (-B);
    }

    matrix operator-=(const matrix & M){
        *this = *this + (-M);
        return *this;
    }

    matrix operator*(const matrix & B) const{
        if(n == B.m){
            matrix<entrada> C(m, B.n);
            for(int i = 0; i < m; i++){
                for(int j = 0; j < B.n; j++){
                    for(int k = 0; k < n; k++){
                        C[i][j] += A[i][k] * B.A[k][j];
                    }
                }
            }
            return C;
        }else{
            return *this;
        }
    }

    matrix operator*(const entrada & c) const{
        matrix<entrada> C(m, n);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                C[i][j] = A[i][j] * c;
            }
        }
        return C;
    }

    matrix operator*=(const matrix & M){
        *this = *this * M;
        return *this;
    }

    matrix operator*=(const entrada & c){
        *this = *this * c;
        return *this;
    }

    bool operator==(const matrix & B) const{
        if(m == B.m && n == B.n){
            for(int i = 0; i < m; i++){
                for(int j = 0; j < n; j++){
                    if(A[i][j] != B.A[i][j]) return false;
                }
            }
            return true;
        }else{
            return false;
        }
    }

    bool operator!=(const matrix & B) const{
        return !(*this == B);
    }

    matrix<entrada> escalonada_reducida_por_filas(){
        matrix<entrada> asoc = *this;
        asoc.gauss_jordan();
        return asoc;
    }

    matrix<entrada> escalonada_por_filas(){
        matrix<entrada> asoc = *this;
        asoc.eliminacion_gaussiana();
        return asoc;
    }

    bool invertible(){
        if(m == n){
            return escalonada_reducida_por_filas() == matrix<entrada>::identidad(n);
        }else{
            return false;
        }
    }

    matrix<entrada> transpuesta(){
        matrix<entrada> T(n, m);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++){
                T[j][i] = A[i][j];
            }
        }
        return T;
    }

    matrix<entrada> inversa(){
        if(m == n){
            matrix<entrada> tmp = *this;
            matrix<entrada> inv = matrix<entrada>::identidad(n);
            tmp.gauss_jordan(inv);
            if(tmp == matrix<entrada>::identidad(n)){
                return inv;
            }else{
                return *this;
            }
        }else{
            return *this;
        }
    }

    entrada determinante(){
        if(m == n){
            matrix<entrada> tmp = *this;
            return tmp.eliminacion_gaussiana();
        }else{
            return 0;
        }
    }

    matrix<entrada> menor(int x, int y){
        matrix<entrada> M(0, 0);
        for(int i = 0; i < m; i++){
            if(i != x){
                M.A.push_back(vector<entrada>());
                for(int j = 0; j < n; j++){
                    if(j != y){
                        M.A[M.A.size() - 1].push_back(A[i][j]);
                    }
                }
            }
        }
        M.m = m - 1;
        M.n = n - 1;
        return M;
    }

    entrada cofactor(int x, int y){
        entrada ans = menor(x, y).determinante();
        if((x + y) % 2 == 1) ans *= -1;
        return ans;
    }

    matrix<entrada> adjunta(){
        return inversa() * determinante();
    }

    matrix<entrada> cofactores(){
        return adjunta().transpuesta();
    }

    string str() const{
        stringstream ss;
        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++) ss << A[i][j] << " ";
            ss << "\n";
        }
        ss << "\n";
        return ss.str();
    }

};

template <typename entrada>
ostream &operator<<(ostream &os, const matrix<entrada> & M) { 
    return os << M.str();
}

void pedirValores(matrix<fraccion> & S){
    for(int i = 0; i < S.m; i++){
        cout << "Introduce la fila " << (i + 1) << ": ";
        for(int j = 0; j < S.n; j++){
            cin >> S.A[i][j];
        }
    }
}

void pedirValores(matrix<enteroModular> & S, ull p){
    ull valor;
    for(int i = 0; i < S.m; i++){
        cout << "Introduce la fila " << (i + 1) << ": ";
        for(int j = 0; j < S.n; j++){
            cin >> valor;
            S.A[i][j] = enteroModular(valor, p);
        }
        cin >> valor;
    }
}

int main()
{
    /*int m, n;
    ull p;
    string campo;
    cout << "Introduce el n\243mero de filas: ";
    cin >> m;
    cout << "Introduce el n\243mero de columnas: ";
    cin >> n;
    cout << "Introduce Q para trabajar en los racionales, o un numero primo p para trabajar en F_p: ";
    cin >> campo;
    if(campo == "Q"){
        matrix<fraccion> M(m, n);
        pedirValores(M);
        matrix<fraccion> I_n = matrix<fraccion>::identidad(n);
        cout << "Determinante: " << M.determinante();
    }else{
        istringstream(campo) >> p;
        matrix<enteroModular> M(m, n);
        pedirValores(M, p);
    }*/
    int n = 202;
    matrix<fraccion> A(n, n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A.A[i][j] = (i == j ? 0 : 1);
        }
    }
    cout << A.determinante();
    return 0;
}