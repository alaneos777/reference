#include <bits/stdc++.h>
#include "numberTheory.cpp"
#include "fraccion.cpp"
using namespace std;
typedef complex<long double> comp;

bool igual(long double a, long double b){
    return abs(a - b) < 1e-9;
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
    polinomio(ll x0){
        coeficientes.push_back(x0);
    }
    polinomio(fraccion x0){
        coeficientes.push_back(x0);
    }
    size_t grado() const{
        return coeficientes.size() - 1;
    }
    void quitar_ceros(){
        while(coeficientes.size() && coeficientes.back() == 0){
            coeficientes.pop_back();
        }
        if(coeficientes.size() == 0) coeficientes.push_back(0);
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
    polinomio operator*(const fraccion& f) const{
        size_t g = grado();
        vector<fraccion> nuevos(g + 1);
        for(size_t i = 0; i <= g; i++) nuevos[i] = f * coeficientes[i];
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
    polinomio operator*=(const fraccion& f){
        *this = *this * f;
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

vector<polinomio> bezout_polinomio(polinomio & A, polinomio & B){
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

polinomio generar_ecuacion(vector<fraccion> raices){
    polinomio ans = 1;
    for(fraccion & raiz : raices) ans *= polinomio(vector<fraccion>{-raiz, 1});
    return ans;
}

polinomio interpolar(vector< pair<fraccion, fraccion> > puntos){
    polinomio ans = 0;
    for(size_t i = 0; i < puntos.size(); i++){
        fraccion k = puntos[i].second;
        polinomio p = 1;
        for(size_t j = 0; j < puntos.size(); j++){
            if(i != j){
                p *= polinomio(vector<fraccion>{-puntos[j].first, 1});
                k /= puntos[i].first - puntos[j].first;
            }
        }
        ans += p * k;
    }
    return ans;
}

polinomio cyclotomic(ll n){
    polinomio num = 1;
    polinomio den = 1;
    for(ll d : divisors[n]){
        ll pot = mu(n / d);
        vector<fraccion> coef(d + 1);
        coef[d] = 1;
        coef[0] = -1;
        if(pot == 1){
            num *= polinomio(coef);
        }else if(pot == -1){
            den *= polinomio(coef);
        }
    }
    return (num / den).first;
}

vector<comp> ec_1(comp a, comp b){ //ax+b=0
    return {-b / a};
}

vector<comp> ec_2(comp a, comp b, comp c){ //ax^2+bx+c=0
    comp D = b * b - 4.0L * a * c;
    D = polar(sqrt(abs(D)), arg(D)/2);
    vector<comp> raices = {(-b + D)/(2.0L * a), (-b - D)/(2.0L * a)};
    return raices;
}

vector<comp> ec_3(comp a, comp b, comp c, comp d){ //ax^3+bx^2+cx+d=0
    //ecuación auxiliar: hacemos x=y-b/3a, y^3+Py+Q=0
    comp P = (3.0L * a * c - b * b)/(3.0L * a * a);
    comp Q = (2.0L * b * b * b - 9.0L * a * b * c + 27.0L * a * a * d)/(27.0L * a * a * a);
    //cuadrática auxiliar: z^2+Qz-(P/3)^3=0
    vector<comp> aux = ec_2(1.0L, Q, -pow(P / 3.0L, 3));
    //u=cbrt(z1), v=cbrt(z2)
    for(size_t i = 0; i < 2; i++) aux[i] = polar(cbrt(abs(aux[i])), arg(aux[i]) / 3);
    comp u, v, w = comp(-0.5, sqrt(0.75));
    //y=u+v
    vector<comp> raices;
    if(igual(abs(P), 0)){
        raices = {aux[0] + aux[1] - b / (3.0L * a), aux[0] * w + aux[1] * w * w - b / (3.0L * a), aux[0] * w * w + aux[1] * w - b / (3.0L * a)};
    }else{
        for(size_t i = 0; i < 3; i++){
            for(size_t j = 0; j < 3; j++){
                u = aux[0] * pow(w, i), v = aux[1] * pow(w, j);
                if(igual(abs(3.0L * u * v + P), 0)) raices.push_back(u + v - b / (3.0L * a));
            }
        }
    }
    return raices;
}

vector<comp> ec_4(comp a, comp b, comp c, comp d, comp e){ //ax^4+bx^3+cx^2+dx+e=0
    //ecuacion auxiliar: hacemos x=y-b/4a, y^4+Ax^2+Bx+C=0
    comp A = (8.0L * a * c - 3.0L * b * b)/(8.0L * a * a);
    comp B = (b * b * b - 4.0L * a * b * c + 8.0L * a * a * d)/(8.0L * a * a * a);
    comp C = (16.0L * a * b * b * c - 3.0L * b * b * b * b - 64.0L * a * a * b * d + 256.0L * a * a * a * e)/(256.0L * a * a * a * a);
    //cúbica auxiliar: p^3+2Ap^2+(A^2-4C)p-B^2=0
    vector<comp> aux = ec_3(1.0L, 2.0L * A, A * A - 4.0L * C, -B * B);
    for(int i = 0 ; i < 3; i++) aux[i] = polar(sqrt(abs(aux[i])), arg(aux[i]) / 2);
    comp p1, p2, p3, w = -1;
    //y=(p1+p2+p3)/2
    vector<comp> raices;
    if(igual(abs(B), 0)){
        raices = {(aux[0] + aux[1] - aux[2]) / 2.0L - b / (4.0L * a), (aux[0] - aux[1] + aux[2]) / 2.0L - b / (4.0L * a), (-aux[0] + aux[1] + aux[2]) / 2.0L - b / (4.0L * a), (-aux[0] - aux[1] - aux[2]) / 2.0L - b / (4.0L * a)};
    }else{
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k < 2; k++){
                    p1 = aux[0] * pow(w, i), p2 = aux[1] * pow(w, j), p3 = aux[2] * pow(w, k);
                    if(igual(abs(p1 * p2 * p3 + B), 0)) raices.push_back((p1 + p2 + p3) / 2.0L - b / (4.0L * a));
                }
            }
        }
    }
    return raices;
}

/*int main(){
    comp a, b, c, d;
    cin >> a >> b >> c >> d;
    vector<comp> ans = ec_3(a, b, c, d);
    function<long double(comp, comp)> cross = [&](comp a, comp b){
        return imag(conj(a) * b);
    };
    for(comp root : ans){
        cout << root << "\n";
    }
    long double area = abs((cross(ans[0], ans[1]) + cross(ans[1], ans[2]) + cross(ans[2], ans[0])) * 0.5L);
    cout << fixed << setprecision(10) << area << "\n";
    return 0;
}*/