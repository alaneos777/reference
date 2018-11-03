#include <bits/stdc++.h>
#include "enteroModular.cpp"
#include "fraccion.cpp"
using namespace std;
typedef long long int lli;

template <typename T>
struct matrix{
	vector<vector<T>> A;
	int m, n;

	matrix(int m, int n): m(m), n(n){
		A.resize(m, vector<T>(n, 0));
	}

	vector<T> & operator[] (int i){
		return A[i];
	}

	const vector<T> & operator[] (int i) const{
		return A[i];
	}

	static matrix identity(int n){
		matrix<T> id(n, n);
		for(int i = 0; i < n; i++)
			id[i][i] = 1;
		return id;
	}

	matrix operator+(const matrix & B) const{
		assert(m == B.m && n == B.n); //same dimensions
		matrix<T> C(m, n);
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++)
				C[i][j] = A[i][j] + B[i][j];
		return C;
	}

	matrix operator+=(const matrix & M){
		*this = *this + M;
		return *this;
	}

	matrix operator-() const{
		matrix<T> C(m, n);
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++)
				C[i][j] = -A[i][j];
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
		assert(n == B.m); //#columns of 1st matrix = #rows of 2nd matrix
		matrix<T> C(m, B.n);
		for(int i = 0; i < m; i++)
			for(int j = 0; j < B.n; j++)
				for(int k = 0; k < n; k++)
					C[i][j] += A[i][k] * B[k][j];
		return C;
	}

	matrix operator*(const T & c) const{
		matrix<T> C(m, n);
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++)
				C[i][j] = A[i][j] * c;
		return C;
	}

	matrix operator*=(const matrix & M){
		*this = *this * M;
		return *this;
	}

	matrix operator*=(const T & c){
		*this = *this * c;
		return *this;
	}

	matrix operator^(lli b) const{
		matrix<T> ans = matrix<T>::identity(n);
		matrix<T> A = *this;
		while(b){
			if(b & 1) ans *= A;
			b >>= 1;
			if(b) A *= A;
		}
		return ans;
	}

	matrix operator^=(lli n){
		*this = *this ^ n;
		return *this;
	}

	bool operator==(const matrix & B) const{
		if(m != B.m || n != B.n) return false;
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++)
				if(A[i][j] != B[i][j]) return false;
		return true;
	}

	bool operator!=(const matrix & B) const{
		return !(*this == B);
	}

	void scaleRow(int k, T c){
		for(int j = 0; j < n; j++)
			A[k][j] *= c;
	}

	void swapRows(int k, int l){
		swap(A[k], A[l]);
	}

	void addRow(int k, int l, T c){
		for(int j = 0; j < n; j++)
			A[k][j] += c * A[l][j];
	}

	matrix<T> transpose(){
		matrix<T> tr(n, m);
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++)
				tr[j][i] = A[i][j];
		return tr;
	}

	T trace(){
		T sum = 0;
		for(int i = 0; i < min(m, n); i++)
			sum += A[i][i];
		return sum;
	}

	//full: true: reduce above and below the diagonal, false: reduce only below
	//makeOnes: true: make the elements in the diagonal ones, false: leave the diagonal unchanged
	//For every elemental operation that we apply to the matrix,
	//we will call to callback(operation, k, l, value).
	//operation 1: multiply row "k" by "value"
	//operation 2: swap rows "k" and "l"
	//operation 3: add "value" times the row "l" to the row "k"
	//It returns the rank of the matrix, and modifies it
	int gauss_jordan(bool full = true, bool makeOnes = true, function<void(int, int, int, T)>callback = NULL){
		int i = 0, j = 0;
		while(i < m && j < n){
			if(A[i][j] == 0){
				for(int f = i + 1; f < m; f++){
					if(A[f][j] != 0){
						swapRows(i, f);
						if(callback) callback(2, i, f, 0);
						break;
					}
				}
			}
			if(A[i][j] != 0){
				T inv_mult = A[i][j].inverso();
				if(makeOnes && A[i][j] != 1){
					scaleRow(i, inv_mult);
					if(callback) callback(1, i, 0, inv_mult);
				}
				for(int f = (full ? 0 : (i + 1)); f < m; f++){
					if(f != i && A[f][j] != 0){
						T inv_adit = -A[f][j];
						if(!makeOnes) inv_adit *= inv_mult;
						addRow(f, i, inv_adit);
						if(callback) callback(3, f, i, inv_adit);
					}
				}
				i++;
			}
			j++;
		}
		return i;
	}

	void gaussian_elimination(){
		gauss_jordan(false);
	}

	matrix<T> reducedRowEchelonForm(){
		matrix<T> asoc = *this;
		asoc.gauss_jordan();
		return asoc;
	}

	matrix<T> rowEchelonForm(){
		matrix<T> asoc = *this;
		asoc.gaussian_elimination();
		return asoc;
	}

	bool invertible(){
		assert(m == n); //this is defined only for square matrices 
		matrix<T> tmp = *this;
		return tmp.gauss_jordan(false) == n;
	}

	matrix<T> inverse(){
		assert(m == n); //this is defined only for square matrices 
		matrix<T> tmp = *this;
		matrix<T> inv = matrix<T>::identity(n);
		auto callback = [&](int op, int a, int b, T e){
			if(op == 1){
				inv.scaleRow(a, e);
			}else if(op == 2){
				inv.swapRows(a, b);
			}else if(op == 3){
				inv.addRow(a, b, e);
			}
		};
		assert(tmp.gauss_jordan(true, true, callback) == n); //check non-invertible
		return inv;
	}

	T determinant(){
		assert(m == n); //only square matrices have determinant
		matrix<T> tmp = *this;
		T det = 1;
		auto callback = [&](int op, int a, int b, T e){
			if(op == 1){
				det /= e;
			}else if(op == 2){
				det *= -1;
			}
		};
		if(tmp.gauss_jordan(false, true, callback) != n) det = 0;
		return det;
	}

	matrix<T> minor(int x, int y){
		matrix<T> M(m-1, n-1);
		for(int i = 0; i < m-1; ++i)
			for(int j = 0; j < n-1; ++j)
				M[i][j] = A[i < x ? i : i+1][j < y ? j : j+1];
		return M;
	}

	T cofactor(int x, int y){
		T ans = minor(x, y).determinant();
		if((x + y) % 2 == 1) ans *= -1;
		return ans;
	}

	matrix<T> cofactorMatrix(){
		matrix<T> C(m, n);
		for(int i = 0; i < m; i++)
			for(int j = 0; j < n; j++)
				C[i][j] = cofactor(i, j);
		return C;
	}

	matrix<T> adjugate(){
		if(invertible()) return inverse() * determinant();
		return cofactorMatrix().transpose();
	}

	tuple<matrix<T>, matrix<T>, matrix<T>> PA_LU(){
		matrix<T> U = *this;
		matrix<T> L = matrix<T>::identity(n);
		matrix<T> P = matrix<T>::identity(n);
		auto callback = [&](int op, int a, int b, T e){
			if(op == 2){
				L.swapRows(a, b);
				P.swapRows(a, b);
				L[a][a] = L[b][b] = 1;
				L[a][a + 1] = L[b][b - 1] = 0;
			}else if(op == 3){
				L[a][b] = -e;
			}
		};
		U.gauss_jordan(false, false, callback);
		return {P, L, U};
	}

	vector<T> characteristicPolynomial(){
		matrix<T> M(n, n);
		vector<T> coef(n + 1);
		matrix<T> I = matrix<T>::identity(n);
		coef[n] = 1;
		for(int i = 1; i <= n; i++){
			M = (*this) * M + I * coef[n - i + 1];
			coef[n - i] = -((*this) * M).trace() / i;
		}
		return coef;
	}

	matrix<T> gram_schmidt(){
		//vectors are rows of the matrix (also in the answer)
		//the answer doesn't have the vectors normalized
		matrix<T> B = (*this) * (*this).transpose();
		matrix<T> ans = *this;
		auto callback = [&](int op, int a, int b, T e){
			if(op == 1){
				ans.scaleRow(a, e);
			}else if(op == 2){
				ans.swapRows(a, b);
			}else if(op == 3){
				ans.addRow(a, b, e);
			}
		};
		B.gauss_jordan(false, false, callback);
		return ans;
	}
};

template <typename T>
ostream &operator<<(ostream & os, const matrix<T> & A){
	for(int i = 0; i < A.m; i++){
		for(int j = 0; j < A.n; j++)
			os << A[i][j] << " ";
		os << "\n";
	}
	os << "\n";
	return os;
}

void pedirValores(matrix<fraccion> & S){
	for(int i = 0; i < S.m; i++){
		cout << "Introduce la fila " << (i + 1) << ": ";
		for(int j = 0; j < S.n; j++){
			cin >> S[i][j];
		}
	}
}

void pedirValores(matrix<enteroModular> & S, lli p){
	lli valor;
	for(int i = 0; i < S.m; i++){
		cout << "Introduce la fila " << (i + 1) << ": ";
		for(int j = 0; j < S.n; j++){
			cin >> valor;
			S[i][j] = enteroModular(valor, p);
		}
	}
}

int main()
{
	int m, n;
	lli p;
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
		cout << "\nDeterminante: " << M.determinant() << "\n\n";
		cout << "Inversa:\n" << M.inverse() << "\n";
		cout << "Adjunta:\n" << M.adjugate() << "\n";
		auto LU = M.PA_LU();
		cout << "P:\n" << get<0>(LU) << "L:\n" << get<1>(LU) << "U:\n" << get<2>(LU);
		cout << "Polinomio caracteristico: ";
		vector<fraccion> polinomio = M.characteristicPolynomial();
		for(int i = 0; i < polinomio.size(); i++){
			cout << polinomio[i] << "x^" << i << ", ";
		}
		cout << "\nGram-Schmidt:\n" << M.gram_schmidt();
	}else{
		istringstream(campo) >> p;
		matrix<enteroModular> M(m, n);
		pedirValores(M, p);
		cout << "\nDeterminante: " << M.determinant() << "\n\n";
		cout << "Inversa:\n" << M.inverse() << "\n";
		cout << "Adjunta:\n" << M.adjugate() << "\n";
		auto LU = M.PA_LU();
		cout << "P:\n" << get<0>(LU) << "L:\n" << get<1>(LU) << "U:\n" << get<2>(LU);
		cout << "Polinomio caracteristico: ";
		vector<enteroModular> polinomio = M.characteristicPolynomial();
		for(int i = 0; i < polinomio.size(); i++){
			cout << polinomio[i] << "x^" << i << ", ";
		}
	}
	return 0;
}