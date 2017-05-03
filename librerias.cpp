#include <bits/stdc++.h>

using namespace std;

# define M_PI 3.14159265358979323846
typedef long long int ull;
typedef complex<double> comp;

vector<ull> suma_divisores, primos, phi;
vector< vector<ull> > factores_primos, divisores, pascal;
vector<bool> es_primo;

bool igual(double a, double b){
    return abs(a - b) < 1e-7;
}

comp precision(comp n){
    if(igual(n.real(), 0.0)) n.real(0);
    if(igual(n.imag(), 0.0)) n.imag(0);
    return n;
}

ull piso(ull a, ull b){
    if((a >= 0 && b > 0) || (a < 0 && b < 0)){
        return a / b;
    }else{
        if(a % b == 0) return a / b;
        else return a / b - 1;
    }
}

ull techo(ull a, ull b){
    if((a >= 0 && b > 0) || (a < 0 && b < 0)){
        if(a % b == 0) return a / b;
        else return a / b + 1;
    }else{
        return a / b;
    }
}

ull fast_pow(ull b, ull e){
    ull ans = 1;
    while(e){
        if(e & 1){
            ans *= b;
        }
        e >>= 1;
        b *= b;
    }
    return ans;
}

ull fast_pow_mod(ull b, ull e, ull m){
    ull ans = 1;
    while(e){
        if(e & 1){
            ans = (ans * b) % m;
        }
        e >>= 1;
        b = (b * b) % m;
    }
    return ans;
}

ull gcd(ull a, ull b){
    ull r;
    while(b != 0) r = a % b, a = b, b = r;
    return a;
}

ull lcm(ull a, ull b){
    ull d = gcd(a, b);
    if(a > b) return b * (a / d);
    else return a * (b / d);
}

ull gcd_multiple(list<ull> nums){
    ull a, b;
    while(nums.size() > 1){
        a = nums.back();
        nums.pop_back();
        b = nums.back();
        nums.pop_back();
        nums.push_back(gcd(a, b));
    }
    return nums.back();
}

ull lcm_multiple(list<ull> nums){
    ull a, b;
    while(nums.size() > 1){
        a = nums.back();
        nums.pop_back();
        b = nums.back();
        nums.pop_back();
        nums.push_back(lcm(a, b));
    }
    return nums.back();
}

ull phi_single(ull n){
	ull resultado = n;
    ull bound = sqrt(n);
	for(ull i = 2; i <= bound; i++){
		if(n % i == 0){
			while(n % i == 0) n /= i;
			resultado -= resultado/i;
		}
	}
    if(n > 1) resultado -= resultado/n;
	return resultado;
}

ull carmichael_lambda(ull n){
    if(n == 1) return 1;
    ull ans, a, p;
    list<ull> f;
    for(ull i = 2; i <= sqrt(n); i++){
        if(n%i == 0){
            a = 0;
            while(n%i == 0){
                n /= i;
                a++;
            }
            p = fast_pow(i, a);
            p -= p/i;
            if(a <= 2 || i >= 3) f.push_back(p);
            else f.push_back(p / 2);
        }
    }
    if(n > 1) f.push_back(n - 1);
    return lcm_multiple(f);
}

vector<ull> euclides(ull a, ull b){
    ull x0 = 1, y0 = 0, x1 = 0, y1 = 1, q, r, xn, yn;
    while(b != 0){
        q = a / b, r = a % b;
        xn = x0 - x1 * q, yn = y0 - y1 * q;
        x0 = x1, x1 = xn;
        y0 = y1, y1 = yn;
        a = b, b = r;
    }
    return {a, x0, y0};
}

ull inverso(ull a, ull m){
    ull inv = euclides(a, m)[1];
    if(inv < 0) inv += abs(m);
    return inv;
}

ull inverso2(ull a, ull m){
    return fast_pow_mod(a, phi_single(m) - 1, m);
}

vector<ull> chinese(vector<ull> a, vector<ull> n){
     ull prod = 1, p, ans = 0;
     for(ull ni : n) prod *= ni;
     for(size_t i = 0; i < a.size(); i++){
        p = prod / n[i];
        ans = (ans + (a[i] % n[i]) * inverso(p, n[i]) * p) % prod;
     }
     return {prod, ans};
}

void criba_primos(ull n){
    es_primo.resize(n + 1, true);
    es_primo[0] = es_primo[1] = false;
    for(ull i = 4; i <= n; i += 2){
        es_primo[i] = false;
    }
    ull bound = sqrt(n);
    for(ull i = 3; i <= bound; i += 2){
        if(es_primo[i]){
            for(ull j = i*i; j <= n; j += 2*i){
                es_primo[j] = false;
            }
        }
    }
    for(ull i = 0; i <= n; i++) if(es_primo[i]) primos.push_back(i);
}

void criba_phi(ull n){
    for(ull i = 0; i <= n; i++) phi.push_back(i);
    for(ull i = 1; i <= n; i++){
        if(es_primo[i]){
            for(ull j = i; j <= n; j += i){
                phi[j] -= phi[j] / i;
            }
        }
    }
}

void criba_divisores(ull n){
    suma_divisores.resize(n + 1, 0);
    divisores.resize(n + 1, vector<ull>());
    for(ull i = 1; i <= n; i++){
        for(ull j = i; j <= n; j += i){
            suma_divisores[j] += i;
            divisores[j].push_back(i);
        }
    }
}

void criba_factores_primos(ull n){
    factores_primos.resize(n + 1, vector<ull>());
    for(ull i = 1; i <= n; i++){
        if(es_primo[i]){
            for(ull j = i; j <= n; j += i){
                factores_primos[j].push_back(i);
            }
        }
    }
}

void factorizar_map(ull n, ull m, map<ull, ull> & f){
    if(n == 1) return;
    ull d = 2;
    while(d <= sqrt(n)){
        if(n % d == 0){
            f[d] += m;
            n /= d;
        }else{
            d++;
        }
    }
    f[n] += m;
}

size_t factorizar_criba(ull n, ull m, vector<ull> & f){
    size_t i = 0;
    for(i = 0; i < primos.size(); i++){
        ull d = primos[i];
        if(d > n) break;
        ull pot = 0;
        while(n % d == 0){
            pot++;
            n /= d;
        }
        f[i] += pot * m;
    }
    return i - 1;
}

ull mu_map(ull n){
    if(n == 0) return 0;
    ull ans = 1;
    map<ull, ull> f;
    factorizar_map(n, 1, f);
    for(pair<const ull, ull> & p : f){
        if(p.second > 1) return 0;
        ans *= -1;
    }
    return ans;
}

ull mu_criba(ull n){
    if(n == 0) return 0;
    ull ans = 1;
    vector<ull> f(primos.size() + 1, 0);
    ull max_p = factorizar_criba(n, 1, f);
    for(ull i = 0; i <= max_p; i++){
        if(f[i] > 0){
            if(f[i] > 1) return 0;
            ans *= -1;
        }
    }
    return ans;
}

size_t factorizar_factorial_criba(ull n, ull m, vector<ull> & f){
    size_t i = 0;
    for(i = 0; i < primos.size(); i++){
        ull d = primos[i];
        if(d > n) break;
        ull pot = 0, contador = 1, tmp = 0;
        while(true){
            tmp = n / fast_pow(d, contador);
            if(tmp == 0) break;
            pot += tmp;
            contador++;
        }
        f[i] += pot*m;
    }
    return i - 1;
}

vector<ull> coprimos_map(ull n){
    map<ull, ull> f;
    vector<ull> ans;
    factorizar_map(n, 1, f);
    for(ull i = 1 ; i <= n; i++){
        bool test = true;
        for(pair<const ull, ull> & p : f){
            if(i % p.first == 0){
                test = false;
                break;
            }
        }
        if(test) ans.push_back(i);
    }
    return ans;
}

vector<ull> coprimos_criba(ull n){
    vector<ull> f(primos.size() + 1, 0);
    vector<ull> ans;
    ull max_p = factorizar_criba(n, 1, f);
    for(ull i = 1 ; i <= n; i++){
        bool test = true;
        for(ull j = 0; j <= max_p; j++){
            if(f[j] > 0){
                if(i % primos[j] == 0){
                    test = false;
                    break;
                }
            }
        }
        if(test) ans.push_back(i);
    }
    return ans;
}

vector<ull> coprimos_gcd(ull n){
    vector<ull> ans;
    for(ull i = 1; i <= n; i++){
        if(gcd(n, i) == 1) ans.push_back(i);
    }
    return ans;
}

string decimal_a_base(ull n, ull b){
    string ans = "";
    ull digito;
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

ull base_a_decimal(string n, ull b){
    ull ans = 0;
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

ull ncr(ull n, ull r){
    ull a = max(r, n - r), b = min(r, n - r);
    ull ans = 1;
    for(ull i = n; i >= a + 1; i--){
        ans *= i;
    }
    for(ull i = 1; i <= b; i++){
        ans /= i;
    }
    return ans;
}

void criba_pascal(ull n){
    pascal.resize(n + 1, vector<ull>());
    pascal[0] = {1};
    for(ull i = 1; i <= n; i++){
        pascal[i].resize(i + 1);
        pascal[i][0] = 1;
        for(ull j = 1; j <= i / 2; j++){
            pascal[i][i - j] = pascal[i][j] = pascal[i - 1][j - 1] + pascal[i - 1][j];
        }
        pascal[i][i] = 1;
    }
}

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
        num = x / d, den = y / d;
    }
    fraccion(ull v){
        num = v;
        den = 1;
    }
    fraccion operator+(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return fraccion(num * (f.den / d) + f.num * (den / d), den * (f.den / d));
    }
    fraccion operator-() const{
        return fraccion(-num, den);
    }
    fraccion operator-(const fraccion& f) const{
        return *this + (-f);
    }
    fraccion operator*(const fraccion& f) const{
        return fraccion(num * f.num, den * f.den);
    }
    fraccion operator/(const fraccion& f) const{
        return fraccion(num * f.den, den * f.num);
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
        return (num * (f.den / d) == (den / d) * f.num);
    }
    bool operator!=(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num * (f.den / d) != (den / d) * f.num);
    }
    bool operator >(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num * (f.den / d) > (den / d) * f.num);
    }
    bool operator <(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num * (f.den / d) < (den / d) * f.num);
    }
    bool operator >=(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num * (f.den / d) >= (den / d) * f.num);
    }
    bool operator <=(const fraccion& f) const{
        ull d = gcd(den, f.den);
        return (num * (f.den / d) <= (den / d) * f.num);
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

vector< vector<fraccion> > gauss;

void criba_gauss(ull n){
    gauss.resize(n + 1, vector<fraccion>());
    gauss[0] = {fraccion(1, 1), fraccion()};
    for(ull i = 1; i <=n; i++){
        gauss[i].resize(i + 2);
        for(ull j = 0; j <= i + 1; j++){
            gauss[i][j] = fraccion(pascal[i + 1][j], i + 1);
        }
        gauss[i][i + 1] = gauss[i][i + 1] - fraccion(1, i + 1);
        for(ull j = 0; j <= i - 1; j++){
            fraccion coef = fraccion(pascal[i + 1][j], i + 1);
            vector<fraccion> pj = gauss[j];
            for(size_t k = 0; k < pj.size(); k++){
                gauss[i][i + 1 - (pj.size() - k - 1)] = gauss[i][i + 1 - (pj.size() - k - 1)] - coef * pj[k];
            }
        }
    }
}

string numero_a_romano(ull n){
    ull digito, base = 0;
    string ans = "";
    vector< vector<char> > datos = {{'I', 'V'}, {'X', 'L'}, {'C', 'D'}, {'M', '\0'}};
    do{
        string tmp = "";
        digito = n % 10;
        n /= 10;
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
        ans = tmp + ans;
        base++;
    }while(n != 0);
    return ans;
}

ull romano_a_numero(string n){
    ull ans = 0;
    char actual, anterior;
    bool f = false;
    map<char, int> datos = {{'I', 1}, {'V', 5}, {'X', 10}, {'L', 50}, {'C', 100}, {'M', 1000}};
    for(int i = n.size() - 1; i >= 0; i--){
        actual = n[i];
        if(i > 0) anterior = n[i - 1];
        if(actual == 'V' && anterior=='I') ans+=4, f=true;
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

struct polinomio{
    vector<fraccion> coeficientes;
    polinomio(){
        coeficientes.push_back(0);
    }
    polinomio(vector<fraccion> coef){
        for(fraccion x : coef) coeficientes.push_back(x);
        quitar_ceros();
    }
    polinomio(ull x0){
        coeficientes.push_back(x0);
    }
    polinomio(fraccion x0){
        coeficientes.push_back(x0);
    }
    size_t grado() const{
        return coeficientes.size() - 1;
    }
    void quitar_ceros(){
        if(coeficientes.empty()){
            coeficientes.push_back(0);
        }else{
            if(coeficientes[coeficientes.size() - 1] == 0){
                coeficientes.pop_back();
                quitar_ceros();
            }
        }
    }
    polinomio operator+(const polinomio& p) const{
        size_t g = max(grado(), p.grado());
        vector<fraccion> nuevos(g + 1, 0);
        for(size_t i = 0; i <= g; i++){
            if(i <= coeficientes.size() - 1) nuevos[i] += coeficientes[i];
            if(i <= p.coeficientes.size() - 1) nuevos[i] += p.coeficientes[i];
        }
        return polinomio(nuevos);
    }
    polinomio operator-() const{
        size_t g = grado();
        vector<fraccion> nuevos(g + 1);
        for(size_t i = 0; i <= g; i++) nuevos[i] = -coeficientes[i];
        return polinomio(nuevos);
    }
    polinomio operator-(const polinomio& p) const{
        return *this + (-p);
    }
    polinomio operator*(const polinomio& p) const{
        size_t g = grado() + p.grado();
        vector<fraccion> nuevos(g + 1, 0);
        for(size_t i = 0; i <= grado(); i++){
            for(size_t j = 0; j <= p.grado(); j++){
                nuevos[i + j] += coeficientes[i] * p.coeficientes[j];
            }
        }
        return polinomio(nuevos);
    }
    pair<polinomio, polinomio> operator/(const polinomio& B) const{
        polinomio Q, R;
        Q = 0;
        R = polinomio(coeficientes);
        while(R.grado() >= B.grado() && !(R.grado() == 0 && R.coeficientes[0] == 0)){
            fraccion q = R.coeficientes[R.grado()] / B.coeficientes[B.grado()];
            size_t g = R.grado() - B.grado();
            vector<fraccion> tmp(g + 1, 0);
            tmp[g] = q;
            polinomio nuevo(tmp);
            Q += nuevo;
            R -= B*nuevo;
        }
        return make_pair(Q, R);
    }
    polinomio operator+=(const polinomio& p){
        *this = *this + p;
        return *this;
    }
    polinomio operator-=(const polinomio& p){
        *this = *this - p;
        return *this;
    }
    polinomio operator*=(const polinomio& p){
        *this = *this * p;
        return *this;
    }
    fraccion evaluar(fraccion x0){
        fraccion ans = 0;
        for(auto it = coeficientes.rbegin(); it != coeficientes.rend(); ++it){
            ans = ans * x0 + (*it);
        }
        return ans;
    }
    string str(){
        bool primero = true;
        stringstream exp;
        for(int pot = coeficientes.size() - 1; pot >= 0; pot--){
            fraccion coef = coeficientes[pot];
            if(coef > 0){
                if(primero){
                    primero = false;
                }else{
                    exp << " + ";
                }
            }else if(coef < 0){
                if(primero){
                    exp << "-";
                    primero = false;
                }else{
                    exp << " - ";
                }
            }
            if(coef != 1 && coef != -1 && coef != 0){
                if(coef > 0) exp << coef.str();
                else exp << (-coef).str();
            }else{
                if(pot == 0){
                    if(coef > 0) exp << coef.str();
                    else if(coef < 0) exp << (-coef).str();
                }
            }
            if(coef != 0){
                if(pot == 1){
                    exp << "x";
                }
                else if(pot > 1){
                    exp << "x^" << pot;
                }
            }
        }
        string linea = exp.str();
        if(linea.size() == 0) return "0";
        else return linea;
    }
};

vector<polinomio> bezout_polinomio(polinomio A, polinomio B){
    polinomio Q, R, S0 = 1, T0 = 0, S1 = 0, T1 = 1, Si, Ti;
    while(!(B.grado() == 0 && B.coeficientes[0] == 0)){
        pair<polinomio, polinomio> div = A / B;
        Q = div.first, R = div.second;
        Si = S0 - S1 * Q, Ti = T0 - T1 * Q;
        S0 = S1, S1 = Si;
        T0 = T1, T1 = Ti;
        A = B, B = R;
    }
    return {A, S0, T0};
}

polinomio generar_ecuacion(vector<ull> raices){
    polinomio ans = 1;
    for(ull & raiz : raices) ans *= polinomio(vector<fraccion>{-raiz, 1});
    return ans;
}

vector<comp> ec_1(comp a, comp b){ //ax+b=0
    return {precision(-b / a)};
}

vector<comp> ec_2(comp a, comp b, comp c){ //ax^2+bx+c=0
    if(igual(abs(c), 0)){
        vector<comp> tmp = ec_1(a, b);
        tmp.push_back(0);
        return tmp;
    }
    comp D = b*b - 4.0*a*c;
    D = polar(sqrt(abs(D)), arg(D)/2);
    vector<comp> raices = {precision((-b + D)/(2.0 * a)), precision((-b - D)/(2.0 * a))};
    return raices;
}

vector<comp> ec_3(comp a, comp b, comp c, comp d){ //ax^3+bx^2+cx+d=0
    if(igual(abs(d), 0)){
        vector<comp> tmp = ec_2(a, b, c);
        tmp.push_back(0);
        return tmp;
    }
    //ecuación auxiliar: hacemos x=y-b/3a, y^3+Py+Q=0
    comp P = (3.0 * a * c - b * b)/(3.0 * a * a);
    comp Q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d)/(27.0 * a * a * a);
    //cuadrática auxiliar: z^2+Qz-(P/3)^3=0
    vector<comp> aux = ec_2(1, Q, -pow(P / 3.0, 3));
    //u=cbrt(z1), v=cbrt(z2)
    for(size_t i = 0; i < 2; i++) aux[i] = polar(cbrt(abs(aux[i])), arg(aux[i]) / 3);
    comp u, v, w = comp(-0.5, sqrt(0.75));
    //y=u+v
    vector<comp> raices;
    if(igual(abs(P), 0)){
        raices = {precision(aux[0] + aux[1] - b / (3.0 * a)), precision(aux[0] * w + aux[1] * w * w - b / (3.0 * a)), precision(aux[0] * w * w + aux[1] * w - b / (3.0 * a))};
    }else{
        for(size_t i = 0; i < 3; i++){
            for(size_t j = 0; j < 3; j++){
                u = aux[0] * pow(w, i), v = aux[1] * pow(w, j);
                if(igual(abs(3.0 * u * v + P), 0)) raices.push_back(precision(u + v - b / (3.0 * a)));
            }
        }
    }
    return raices;
}

vector<comp> ec_4(comp a, comp b, comp c, comp d, comp e){ //ax^4+bx^3+cx^2+dx+e=0
    if(igual(abs(e), 0)){
        vector<comp> tmp = ec_3(a, b, c, d);
        tmp.push_back(0);
        return tmp;
    }
    //ecuacion auxiliar: hacemos x=y-b/4a, y^4+Ax^2+Bx+C=0
    comp A = (8.0 * a * c - 3.0 * b * b)/(8.0 * a * a);
    comp B = (b * b * b - 4.0 * a * b * c + 8.0 * a * a * d)/(8.0 * a * a * a);
    comp C = (16.0 * a * b * b * c - 3.0 * b * b * b * b - 64.0 * a * a * b * d + 256.0 * a * a * a * e)/(256.0 * a * a * a * a);
    //cúbica auxiliar: p^3+2Ap^2+(A^2-4C)p-B^2=0
    vector<comp> aux = ec_3(1, 2.0 * A, A * A - 4.0 * C, -B * B);
    for(int i = 0 ; i < 3; i++) aux[i] = polar(sqrt(abs(aux[i])), arg(aux[i]) / 2);
    comp p1, p2, p3, w = -1;
    //y=(p1+p2+p3)/2
    vector<comp> raices;
    if(igual(abs(B), 0)){
        raices = {precision((aux[0] + aux[1] - aux[2]) / 2.0 - b / (4.0 * a)), precision((aux[0] - aux[1] + aux[2]) / 2.0 - b / (4.0 * a)), precision((-aux[0] + aux[1] + aux[2]) / 2.0 - b / (4.0 * a)), precision((-aux[0] - aux[1] - aux[2]) / 2.0 - b / (4.0 * a))};
    }else{
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    p1 = aux[0] * pow(w, i), p2 = aux[1] * pow(w, j), p3 = aux[2] * pow(w, k);
                    if(igual(abs(p1 * p2 * p3 + B), 0)) raices.push_back(precision((p1 + p2 + p3) / 2.0 - b / (4.0 * a)));
                }
            }
        }
    }
    return raices;
}

typedef fraccion entrada;
typedef vector<entrada> fila;
typedef vector<fila> matrix;

void imprime_matriz(matrix matriz){
    for(fila & renglon : matriz){
        for(entrada & valor : renglon){
            cout << valor.str() << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

matrix transpuesta(matrix matriz){
    int m = matriz.size(), n = matriz[0].size();
    matrix ans(n, fila(m));
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            ans[j][i] = matriz[i][j];
        }
    }
    return ans;
}

matrix suma_matrices(matrix A, matrix B){
    int m = A.size(), n = A[0].size();
    matrix ans(m, fila(n));
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
            ans[i][j] = A[i][j] + B[i][j];
        }
    }
    return ans;
}

matrix mult_matrices(matrix A, matrix B){
    int m = A.size(), n = B.size(), p = B[0].size();
    matrix C(m, fila(p, 0));
    for(int i = 0; i < m; i++){
        for(int j = 0; j < p; j++){
            for(int k = 0; k < n; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

matrix eliminacion_gaussiana(matrix matriz, entrada & p, bool salir){
    int m = matriz.size();
    int n = matriz[0].size();
    int offset = 0;
    for(int j=0;j-offset<m && j<n;j++){
        int i = j - offset;
        if(matriz[i][j] == 0){
            for(int f=i+1;f<m;f++){
                if(matriz[f][j] != 0){
                    swap(matriz[f], matriz[i]);
                    p *= -1;
                    break;
                }
            }
        }
        if(matriz[i][j] == 0){
            if(salir){
                p = 0;
                return matriz;
            }else{
                offset++;
                continue;
            }
        }
        for(int f=i+1;f<m;f++){
            if(matriz[f][j] != 0){
                fraccion c = matriz[f][j]/matriz[i][j];
                for(int l=n-1;l>=j;l--){
                    matriz[f][l] -= c*matriz[i][l];
                }
            }
        }
    }
    return matriz;
}

matrix eliminacion_gaussiana(matrix matriz){
    fraccion xd = 0;
    return eliminacion_gaussiana(matriz, xd, false);
}

matrix gauss_jordan(matrix matriz){
    int m = matriz.size();
    int n = matriz[0].size();
    int i = 0, j = 0;
    while(i < m && j < n){
        if(matriz[i][j] == 0){
            for(int f = i + 1; f < m; f++){
                if(matriz[f][j] != 0){
                    swap(matriz[f], matriz[i]);
                    break;
                }
            }
        }
        if(matriz[i][j] != 0){
            for(int l = n - 1; l >= j; l--) matriz[i][l] /= matriz[i][j];
            for(int f = 0; f < m; f++) if(f != i && matriz[f][j] != 0) for(int l = n - 1; l >= j; l--) matriz[f][l] -= matriz[f][j] * matriz[i][l];
            i++;
        }
        j++;
    }
    return matriz;
}

entrada determinante(matrix matriz){
    entrada ans = 1;
    matriz = eliminacion_gaussiana(matriz, ans, true);
    if(ans != 0) for(size_t i = 0; i < matriz.size(); i++) ans *= matriz[i][i];
    return ans;
}

matrix matriz_inversa(matrix matriz){
    int m = matriz.size();
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            if(i == j) matriz[i].push_back(1);
            else matriz[i].push_back(0);
        }
    }
    matriz = gauss_jordan(matriz);
    matrix ans(m, fila(m));
    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            ans[i][j] = matriz[i][j + m];
        }
    }
    return ans;
}

matrix matriz_adjunta(matrix matriz){
    entrada det = determinante(matriz);
    if(det == 0) return matriz;
    matriz = matriz_inversa(matriz);
    for(size_t i = 0; i < matriz.size(); i++){
        for(size_t j = 0; j < matriz.size(); j++){
            matriz[i][j] *= det;
        }
    }
    return matriz;
}

matrix matriz_cofactores(matrix matriz){
    return transpuesta(matriz_adjunta(matriz));
}

matrix menor(int i, int j, matrix matriz){
    matrix ans;
    int s = matriz.size();
    for(int m = 0; m<s; m++){
        if(m != i){
            ans.push_back(fila());
            for(int n = 0; n < s; n++){
                if(n != j){
                    ans[ans.size() - 1].push_back(matriz[m][n]);
                }
            }
        }
    }
    return ans;
}

entrada cofactor(int i, int j, matrix matriz){
    entrada ans = 1;
    if((i + j) % 2 == 1) ans = -1;
    ans *= determinante(menor(i, j, matriz));
    return ans;
}

polinomio cyclotomic(ull n){
    vector<comp> raices;
    for(int i = 1; i <= n;i++){
        if(gcd(i, n) == 1){
            raices.push_back(polar(1.0, 2.0 * M_PI * i / n));
        }
    }
    int deg = raices.size();
    vector<comp> ans(deg + 1, 0);
    ans[0] = 1;
    for(int i = 0; i < deg; i++){
        vector<comp> tmp(i + 1);
        for(int j = 0; j < i + 1; j++) tmp[j] = -raices[i] * ans[j];
        for(size_t j = 0; j < tmp.size(); j++) ans[j + 1] += tmp[j];
    }
    vector<fraccion> coef(deg + 1);
    for(int i = 0; i <= deg; i++){
        coef[deg - i] = round(real(ans[i]));
    }
    return polinomio(coef);
}

int main()
{
    criba_divisores(101);
    criba_primos(101);
    criba_phi(101);
    criba_factores_primos(101);
    criba_pascal(50);
    criba_gauss(20);
    /*map<ull, ull> f;
    factorizar_map(63, 1, &f);
    for(pair<ull, ull> p:f) cout << p.first << " " << p.second << endl;
    vector<ull> f2(primos.size()+1, 0);
    cout << endl;
    ull max_p = factorizar_criba(63, 1, &f2);
    for(ull i=0;i<=max_p;i++) if(f2[i]>0) cout << primos[i] << " " << f2[i] << endl;*/

    //for(ull i=0;i<=101;i++) cout << "mu(" << i << ") = " << mu_map(i) << " = " << mu_criba(i) << endl;

    /*vector<ull> f3(primos.size()+1, 0);
    ull max_p = factorizar_factorial_criba(59, 1, &f3);
    for(ull i=0;i<=max_p;i++) if(f3[i]>0) cout << primos[i] << " " << f3[i] << endl;*/

    /*vector<ull> c = coprimos_criba(60);
    for(ull ci:c) cout << ci << endl;*/

    /*for(ull i=1;i<=100;i++){
        vector<ull> pol = cyclotomic(i);
        cout << "Phi_(" << i << ")(x) = " << str_polinomio(pol) << endl;
    }*/

    /*for(ull i=1;i<factores_primos.size();i++){
        cout << i << ": ";
        for(ull fi:factores_primos[i]) cout << fi << " ";
        cout << endl;
    }*/

    /*for(ull i=1;i<divisores.size();i++){        int g = max(grado(), p.grado());
        vector<fraccion> nuevo(g+1, fraccion());
        for(int i=0;i<=g;i++){
            if(i<=polinomio.coeficientes.size()-1) nuevo[i] += polinomio.coeficientes[i];
            if(i<=p.coeficientes.size()-1) nuevo[i] += p.coeficientes[i];
        }
        return polinomio(nuevo);
        cout << i << ": ";
        for(ull fi:divisores[i]) cout << fi << " ";
        cout << endl;
    }*/

    //cout << decimal_a_base(494159, 36) << endl << base_a_decimal("alan", 36);

    /*for(vector<ull> fila:pascal){
        for(ull valor:fila){
            cout << valor << " ";
        }
        cout << endl;
    }*/

    /*for(size_t i=0;i<gauss.size();i++){
        cout << "P_" << i << " = ";
        for(size_t j=0;j<gauss[i].size();j++){
            cout << "+(" << gauss[i][j].num << "/" << gauss[i][j].den << ")n^" << (gauss[i].size()-j-1);
        }
        cout << endl;
    }*/

    //cout << numero_a_romano(456) << endl << romano_a_numero("CDLVI");

    /*comp a, b, c, d;
    cin >> a >> b >> c >> d;
    vector<comp> ec = ec_3(a,b,c,d);
    for(comp raiz:ec) cout << raiz << " ";*/

    /*comp a, b, c, d, e;
    cin >> a >> b >> c >> d >> e;
    vector<comp> ec = ec_4(a,b,c,d,e);
    for(comp raiz:ec) cout << raiz << " ";*/

    //cout << phi_single(100);

    //cout << inverso2(73, 1508);

    /*polinomio A({5, 9, 7, -4, 1});
    polinomio B({-2, 1, 7, 3});
    pair<polinomio, polinomio> info = A/B;
    cout << A.str() << " = (" << B.str() << ")(" << info.first.str() << ") + (" << info.second.str() << ")" << endl;
    polinomio ec = generar_ecuacion({1, 2, 3, -6});
    cout << "P(x)=" << ec.str() << endl;
    cout << cyclotomic(7).str() << endl;*/

    /*polinomio A(vector<fraccion>{-6, 11, -6, 1});
    polinomio B(vector<fraccion>{-24, 26, -9, 1});
    vector<polinomio> info = bezout_polinomio(A, B);
    cout << "(" << info[1].str() << ")(" << A.str() << ") + (" << info[2].str() << ")(" << B.str() << ") = " << info[0].str();*/

    /*vector< vector<fraccion> > C = mult_matrices({{5,8,4,9},{6,1,7,2}}, {{3,4},{6,8},{1,5},{9,2}});
    imprime_matriz(C);*/

    vector< vector<fraccion> > o {{2,-3,-1,4,-27},{-3,5,6,2,44},{1,3,-8,9,-13},{7,-4,5,3,-34}};
    o = eliminacion_gaussiana(o);
    imprime_matriz(o);

    vector< vector<fraccion> > p {{2,-3,-1,4,-27},{-3,5,6,2,44},{1,3,-8,9,-13},{7,-4,5,3,-34}};
    p = gauss_jordan(p);
    imprime_matriz(p);

    vector< vector<fraccion> > q {{2,-3,-1,4},{-3,5,6,2},{1,3,-8,9},{7,-4,5,3}};
    cout << "D=" << determinante(q).str() << endl << endl;

    vector< vector<fraccion> > r {{2,-3,-1,4},{-3,5,6,2},{1,3,-8,9},{7,-4,5,3}};
    r = matriz_inversa(r);
    imprime_matriz(r);

    vector< vector<fraccion> > s {{2,-3,-1,4},{-3,5,6,2},{1,3,-8,9},{7,-4,5,3}};
    s = matriz_adjunta(s);
    imprime_matriz(s);

    vector< vector<fraccion> > t {{2,-3,-1,4},{-3,5,6,2},{1,3,-8,9},{7,-4,5,3}};
    t = matriz_cofactores(t);
    imprime_matriz(t);

    vector< vector<fraccion> > u {{2,-3,-1,4},{-3,5,6,2},{1,3,-8,9},{7,-4,5,3}};
    for(int i=0;i<=3;i++) for(int j=0;j<=3;j++) cout << "C(" << i << "," << j << ")=" << cofactor(i, j, u).str() << "\n";

    /*imprime_matriz(eliminacion_gaussiana({{2,-5,3,4},{1,-2,1,3},{5,1,7,11}}));

    imprime_matriz(gauss_jordan({{1,5,1,1},{-2,4,7,9}}));

    imprime_matriz(gauss_jordan({{0,1,1},{1,2,4},{-7,3,-11},{2,2,6}}));*/

    /*cout << "\n\n\n";
    for(ull i = 1; i <= 100; i++) cout << carmichael_lambda(i) << ", "; */

    return 0;
}
