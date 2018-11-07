#include <bits/stdc++.h>
using namespace std;
typedef long long int lli;

int lis(vector<int> & arr){
	if(arr.size() == 0) return 0;
	vector<int> aux(arr.size());
	int ans = 1;
	aux[0] = arr[0];
	for(int i = 1; i < arr.size(); ++i){
		if(arr[i] < aux[0])
			aux[0] = arr[i];
		else if(arr[i] > aux[ans - 1])
			aux[ans++] = arr[i];
		else
			aux[lower_bound(aux.begin(), aux.begin() + ans, arr[i]) - aux.begin()] = arr[i];
	}
	return ans;
}

int lcs(string & a, string & b){
	int m = a.size(), n = b.size();
	vector<vector<int>> aux(m + 1, vector<int>(n + 1));
	for(int i = 1; i <= m; ++i){
		for(int j = 1; j <= n; ++j){
			if(a[i - 1] == b[j - 1])
				aux[i][j] = 1 + aux[i - 1][j - 1];
			else
				aux[i][j] = max(aux[i - 1][j], aux[i][j - 1]);
		}
	}
	return aux[m][n];
}

//0:saturday, 1:sunday, ..., 6:friday
int dayOfWeek(int d, int m, lli y){
	if(m == 1 || m == 2){
		m += 12;
		--y;
	}
	int k = y % 100;
	lli j = y / 100;
	return (d + 13*(m+1)/5 + k + k/4 + j/4 + 5*j) % 7;
}

//cout for __int128
ostream &operator<<(ostream &os, const __int128 & value){
	char buffer[64];
	char *pos = end(buffer) - 1;
	*pos = '\0';
	__int128 tmp = value < 0 ? -value : value;
	do{
		--pos;
		*pos = tmp % 10 + '0';
		tmp /= 10;
	}while(tmp != 0);
	if(value < 0){
		--pos;
		*pos = '-';
	}
	return os << pos;
}

//cin for __int128
istream &operator>>(istream &is, __int128 & value){
	char buffer[64];
	is >> buffer;
	char *pos = begin(buffer);
	int sgn = 1;
	value = 0;
	if(*pos == '-'){
		sgn = -1;
		++pos;
	}else if(*pos == '+'){
		++pos;
	}
	while(*pos != '\0'){
		value = (value << 3) + (value << 1) + (*pos - '0');
		++pos;
	}
	value *= sgn;
	return is;
}

struct satisfiability_twosat{
	int n;
	vector<vector<int>> imp;

	satisfiability_twosat(int n) : n(n), imp(2 * n) {}

	void add_edge(int u, int v){imp[u].push_back(v);}

	int neg(int u){return (n << 1) - u - 1;}

	void implication(int u, int v){
		add_edge(u, v);
		add_edge(neg(v), neg(u));
	}

	vector<bool> solve(){
		int size = 2 * n;
		vector<int> S, B, I(size);

		function<void(int)> dfs = [&](int u){
			B.push_back(I[u] = S.size());
			S.push_back(u);

			for(int v : imp[u])
				if(!I[v]) dfs(v);
				else while (I[v] < B.back()) B.pop_back();

			if(I[u] == B.back())
				for(B.pop_back(), ++size; I[u] < S.size(); S.pop_back())
					I[S.back()] = size;
		};

		for(int u = 0; u < 2 * n; ++u)
			if(!I[u]) dfs(u);

		vector<bool> values(n);

		for(int u = 0; u < n; ++u)
			if(I[u] == I[neg(u)]) return {};
			else values[u] = I[u] < I[neg(u)];

		return values;
	}
};

//gray code
int gray(int n){
	return n ^ (n >> 1);
}

//inverse gray code
int inv_gray(int g){
	int n = 0;
	while(g){
		n ^= g;
		g >>= 1;
	}
	return n;
}

int LevenshteinDistance(string & a, string & b){
	int m = a.size(), n = b.size();
	vector<vector<int>> aux(m + 1, vector<int>(n + 1));
	for(int i = 1; i <= m; ++i)
		aux[i][0] = i;
	for(int j = 1; j <= n; ++j)
		aux[0][j] = j;
	for(int j = 1; j <= n; ++j)
		for(int i = 1; i <= m; ++i)
			aux[i][j] = min({aux[i-1][j] + 1, aux[i][j-1] + 1, aux[i-1][j-1] + (a[i-1] != b[j-1])});
	return aux[m][n];
}

//count the number of 1's in the i-th bit of all
//representations in binary of numbers in [1,n]
lli count(lli n, int i){
	if(n <= 0) return 0ll;
	lli ans = ((n + 1) >> (i + 1)) << i;
	ans += max(((n + 1) & ((1ll << (i + 1)) - 1)) - (1ll << i), 0ll);
	return ans;
}

/*int main(){
	int a, b, c;
	cin >> a >> b >> c;
	vector<string> days = {"Sabado", "Domingo", "Lunes", "Martes", "Miercoles", "Jueves", "Viernes"};
	cout << days[dayOfWeek(a, b, c)] << "\n";
	return 0;
}*/

int main(){
	string a, b;
	cin >> a >> b;
	cout << LevenshteinDistance(a, b) << "\n";
}