#include <bits/stdc++.h>

using namespace std;

template <typename tipo>
using matriz = vector< vector<tipo> >;

matriz<bool> composicion(matriz<bool> B, matriz<bool> A){
    size_t a = A.size(), b = B.size(), c = B[0].size();
    matriz<bool> C(a, vector<bool>(c, 0));
    for(size_t i=0;i<a;i++){
        for(size_t j=0;j<c;j++){
            for(size_t k=0;k<b;k++){
                C[i][j] = C[i][j] || (A[i][k] && B[k][j]);
            }
        }
    }
    return C;
}

matriz<bool> inversa(matriz<bool> A){
    size_t a = A.size(), b = A[0].size();
    matriz<bool> B(b, vector<bool>(a));
    for(size_t i=0;i<a;i++){
        for(size_t j=0;j<b;j++){
            B[j][i] = A[i][j];
        }
    }
    return B;
}

matriz<bool> unir(matriz<bool> A, matriz<bool> B){
    size_t a = A.size(), b = A[0].size();
    matriz<bool> C(a, vector<bool>(b, 0));
    for(size_t i=0;i<a;i++){
        for(size_t j=0;j<b;j++){
            C[i][j] = A[i][j] || B[i][j];
        }
    }
    return C;
}

matriz<bool> clr(matriz<bool> A){
    for(size_t i=0;i<A.size();i++) A[i][i] = 1;
    return A;
}

matriz<bool> cls(matriz<bool> A){
    for(size_t i=0;i<A.size();i++){
        for(size_t j=0;j<A.size();j++){
            if(A[i][j] == 1) A[j][i] = 1;
        }
    }
    return A;
}

matriz<bool> clt(matriz<bool> A){
    for(size_t i=0;i<A.size();i++){
        for(size_t j=0;j<A.size();j++){
            for(size_t k=0;k<A.size();k++){
                A[j][k] = A[j][k] || (A[j][i] && A[i][k]);
            }
        }
    }
    return A;
}

vector<int> clase_equivalencia(matriz<bool> A, int a){
    vector<int> ans;
    for(size_t i=0;i<A.size();i++){
        if(A[a][i] == 1) ans.push_back(i);
    }
    return ans;
}

vector< vector<int> > conjunto_cociente(matriz<bool> A){
    vector< vector<int> > ans;
    vector<bool> pendientes(A.size(), true);
    for(size_t i=0;i<A.size();i++){
        if(pendientes[i]){
            vector<int> tmp = clase_equivalencia(A, i);
            ans.push_back(tmp);
            for(int x:tmp) pendientes[x] = false;
        }
    }
    return ans;
}

template<typename T>
void imprime_matriz(matriz<T> M, vector<string> A, vector<string> B){
    cout << "\\{";
    bool primero = true;
    for(size_t i=0;i<M.size();i++){
        for(size_t j=0;j<M[i].size();j++){
            if(M[i][j] == 1){
                if(primero){
                    primero = false;
                }else{
                    cout << ",";
                }
                cout << "(" << A[i] << "," << B[j] << ")";
            }
        }
    }
    cout << "\\}$" << endl;
}

int main()
{
    vector<string> A = {"1", "2", "3", "4", "5", "6", "7"};

    matriz<bool> R(7, vector<bool>(7, 0));
    matriz<bool> S(7, vector<bool>(7, 0));
    matriz<bool> T(7, vector<bool>(7, 0));

    R[0][1] = R[0][2] = R[1][2] = R[2][2] = R[2][3] = R[3][1] = R[4][5] = R[6][5] = 1;
    S[0][0] = S[1][4] = S[3][2] = S[5][1] = 1;
    T[0][4] = T[0][6] = T[1][0] = T[1][2] = T[3][5] = T[3][6] = T[4][2] = T[6][6] = 1;
    cout << "R:" << endl;
    imprime_matriz<bool>(R, A, A);
    cout << endl << "S:" << endl;
    imprime_matriz<bool>(S, A, A);
    cout << endl << "T:" << endl;
    imprime_matriz<bool>(T, A, A);

    cout << endl << "\\item $cl_r(\\mathcal{R})=";
    imprime_matriz<bool>(clr(R), A, A);
    cout << endl << "\\item $cl_s(\\mathcal{R})=" << endl;
    imprime_matriz<bool>(cls(R), A, A);
    cout << endl << "\\item $cl_t(\\mathcal{R})=" << endl;
    imprime_matriz<bool>(clt(R), A, A);

    cout << endl << "\\item $cl_r(\\mathcal{S})=" << endl;
    imprime_matriz<bool>(clr(S), A, A);
    cout << endl << "\\item $cl_s(\\mathcal{S})=" << endl;
    imprime_matriz<bool>(cls(S), A, A);
    cout << endl << "\\item $cl_t(\\mathcal{S})=" << endl;
    imprime_matriz<bool>(clt(S), A, A);

    cout << endl << "\\item $cl_r(\\mathcal{T})=" << endl;
    imprime_matriz<bool>(clr(T), A, A);
    cout << endl << "\\item $cl_s(\\mathcal{T})=" << endl;
    imprime_matriz<bool>(cls(T), A, A);
    cout << endl << "\\item $cl_t(\\mathcal{T})=" << endl;
    imprime_matriz<bool>(clt(T), A, A);

    cout << endl << "\\item $cl_r(cl_s(\\mathcal{S}))=" << endl;
    imprime_matriz<bool>(clr(cls(S)), A, A);

    cout << endl << "\\item $cl_r(cl_t(\\mathcal{S}))=" << endl;
    imprime_matriz<bool>(clr(clt(S)), A, A);

    cout << endl << "\\item $cl_s(cl_t(cl_r(\\mathcal{R})))=" << endl;
    imprime_matriz<bool>(cls(clt(clr(R))), A, A);

    cout << endl << "\\item $cl_t(cl_s(cl_r(\\mathcal{T})))=" << endl;
    imprime_matriz<bool>(clt(cls(clr(T))), A, A);

    cout << endl << "\\item $cl_t(cl_s(cl_r(\\mathcal{R})))=" << endl;
    imprime_matriz<bool>(clt(cls(clr(R))), A, A);

    cout << endl;

    vector< vector<int> > quotient = conjunto_cociente(clt(cls(clr(R))));
    for(vector<int> cto:quotient){
        for(int x:cto) cout << x << " ";
        cout << "\n";
    }

    /*matriz<bool> R(4, vector<bool>(5, 0));
    matriz<bool> S(5, vector<bool>(3, 0));
    matriz<bool> T(4, vector<bool>(4, 0));
    matriz<bool> U(5, vector<bool>(5, 0));
    matriz<bool> V(3, vector<bool>(3, 0));

    vector<string> A = {"1", "2", "3", "4"};
    vector<string> B = {"a", "b", "c", "d", "e"};
    vector<string> C = {"\\alpha", "\\beta", "\\gamma"};

    R[0][1] = R[0][2] = R[1][2] = R[2][0] = R[2][4] = R[3][0] = R[3][1] = R[3][3] = R[3][4] = 1; //{(1, b),(1, c),(2, c),(3, a),(3, e),(4, a),(4, b),(4, d),(4, e)}
    S[0][1] = S[0][2] = S[1][2] = S[2][1] = S[4][0] = S[4][1] = 1; //{(a, β),(a, γ),(b, γ),(c, β),(e, α),(e, β)}
    T[0][2] = T[1][1] = T[1][2] = T[1][3] = T[2][0] = T[3][0] = T[3][3] = 1; //{(1, 3),(2, 2),(2, 3),(2, 4),(3, 1),(4, 1),(4, 4)}
    U[0][0] = U[0][4] = U[1][2] = U[2][1] = U[3][0] = U[3][1] = U[3][4] = U[4][0] = U[4][4] = 1; //{(a, a),(a, e),(b, c),(c, b),(d, a),(d, b),(d, e),(e, a),(e, e)}
    V[0][2] = 1; //{(α, γ)}

    cout << "$R="; imprime_matriz<bool>(R, A, B);
    cout << "$S="; imprime_matriz<bool>(S, B, C);
    cout << "$T="; imprime_matriz<bool>(T, A, A);
    cout << "$U="; imprime_matriz<bool>(U, B, B);
    cout << "$V="; imprime_matriz<bool>(V, C, C);
    cout << "\n";
    cout << "$R^{-1}="; imprime_matriz(inversa(R), B, A);
    cout << "$S o R="; imprime_matriz(composicion(S, R), A, C);
    cout << "$R^{-1} o S^{-1}="; imprime_matriz(composicion(inversa(R), inversa(S)), C, A);
    cout << "$R^{-1} o R="; imprime_matriz(composicion(inversa(R), R), A, A);
    cout << "$(R o T)^{-1}="; imprime_matriz(inversa(composicion(R, T)), B, A);
    cout << "$V o (S o R)="; imprime_matriz(composicion(V, composicion(S, R)), A, C);
    cout << "$(U o R) o T="; imprime_matriz(composicion(composicion(U, R), T), A, B);*/
    return 0;
}
