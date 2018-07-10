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

int main(){
	int a, b, c;
	cin >> a >> b >> c;
	vector<string> days = {"Sabado", "Domingo", "Lunes", "Martes", "Miercoles", "Jueves", "Viernes"};
	cout << days[dayOfWeek(a, b, c)] << "\n";
	return 0;
}