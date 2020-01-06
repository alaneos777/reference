#include <bits/stdc++.h>
using namespace std;

struct kmp{
	vector<int> aux;
	string pattern;

	kmp(string pattern){
		this->pattern = pattern;
		aux.resize(pattern.size());
		int i = 1, j = 0;
		while(i < pattern.size()){
			if(pattern[i] == pattern[j])
				aux[i++] = ++j;
			else{
				if(j == 0) aux[i++] = 0;
				else j = aux[j - 1];
			}
		}
	}

	vector<int> search(string & text){
		vector<int> ans;
		int i = 0, j = 0;
		while(i < text.size() && j < pattern.size()){
			if(text[i] == pattern[j]){
				++i, ++j;
				if(j == pattern.size()){
					ans.push_back(i - j);
					j = aux[j - 1];
				}
			}else{
				if(j == 0) ++i;
				else j = aux[j - 1];
			}
		}
		return ans;
	}
};

const int M = 26;
struct node{
	vector<int> child;
	int p = -1;
	char c = 0;
	int suffixLink = -1, endLink = -1;
	int id = -1;

	node(int p = -1, char c = 0) : p(p), c(c){
		child.resize(M, -1);
	}
};

struct AhoCorasick{
	vector<node> t;
	vector<int> lenghts;
	int wordCount = 0;

	AhoCorasick(){
		t.emplace_back();
	}

	void add(const string & s){
		int u = 0;
		for(char c : s){
			if(t[u].child[c-'a'] == -1){
				t[u].child[c-'a'] = t.size();
				t.emplace_back(u, c);
			}
			u = t[u].child[c-'a'];
		}
		t[u].id = wordCount++;
		lenghts.push_back(s.size());
	}

	void link(int u){
		if(u == 0){
			t[u].suffixLink = 0;
			t[u].endLink = 0;
			return;
		}
		if(t[u].p == 0){
			t[u].suffixLink = 0;
			if(t[u].id != -1) t[u].endLink = u;
			else t[u].endLink = t[t[u].suffixLink].endLink;
			return;
		}
		int v = t[t[u].p].suffixLink;
		char c = t[u].c;
		while(true){
			if(t[v].child[c-'a'] != -1){
				t[u].suffixLink = t[v].child[c-'a'];
				break;
			}
			if(v == 0){
				t[u].suffixLink = 0;
				break;
			}
			v = t[v].suffixLink;
		}
		if(t[u].id != -1) t[u].endLink = u;
		else t[u].endLink = t[t[u].suffixLink].endLink;
	}

	void build(){
		queue<int> Q;
		Q.push(0);
		while(!Q.empty()){
			int u = Q.front(); Q.pop();
			link(u);
			for(int v = 0; v < M; ++v)
				if(t[u].child[v] != -1)
					Q.push(t[u].child[v]);
		}
	}

	int match(const string & text){
		int u = 0;
		int ans = 0;
		for(int j = 0; j < text.size(); ++j){
			int i = text[j] - 'a';
			while(true){
				if(t[u].child[i] != -1){
					u = t[u].child[i];
					break;
				}
				if(u == 0) break;
				u = t[u].suffixLink;
			}
			int v = u;
			while(true){
				v = t[v].endLink;
				if(v == 0) break;
				++ans;
				int idx = j + 1 - lenghts[t[v].id];
				cout << "Found word #" << t[v].id << " at position " << idx << "\n";
				v = t[v].suffixLink;
			}
		}
		return ans;
	}
};

struct Node{
  	bool isWord = false;
	map<char, Node*> letters;
};

struct Trie{
	Node* root;

	Trie(){
		root = new Node();
	}

	inline bool exists(Node * actual, const char & c){
		return actual->letters.find(c) != actual->letters.end();
	}

	void InsertWord(const string& word){
		Node* current = root;
		for(auto & c : word){
			if(!exists(current, c))
				current->letters[c] = new Node();
			current = current->letters[c];
		}
		current->isWord = true;
	}

	bool FindWord(const string& word){
		Node* current = root;
		for(auto & c : word){
			if(!exists(current, c))
				return false;
			current = current->letters[c];
		}
		return current->isWord;
	}

	void printRec(Node * actual, string acum){
		if(actual->isWord){
			cout << acum << "\n";
		}
		for(auto & next : actual->letters)
			printRec(next.second, acum + next.first);
	}

	void printWords(const string & prefix){
		Node * actual = root;
		for(auto & c : prefix){
			if(!exists(actual, c)) return;
			actual = actual->letters[c];
		}
		printRec(actual, prefix);
	}
};

struct state{
	int len, link;
	vector<int> child;
	state(int len = 0, int link = -1): len(len), link(link), child(M, -1){}
	state(int len, int link, const vector<int> & child): len(len), link(link), child(child){}
};

struct SuffixAutomaton{
	vector<state> st;
	int last = 0;

	SuffixAutomaton(){
		st.emplace_back();
	}

	void extend(char c){
		int curr = st.size();
		st.emplace_back(st[last].len + 1);
		int p = last;
		while(p != -1 && st[p].child[c-'A'] == -1){
			st[p].child[c-'A'] = curr;
			p = st[p].link;
		}
		if(p == -1){
			st[curr].link = 0;
		}else{
			int q = st[p].child[c-'A'];
			if(st[p].len + 1 == st[q].len){
				st[curr].link = q;
			}else{
				int clone = st.size();
				st.emplace_back(st[p].len + 1, st[q].link, st[q].child);
				while(p != -1 && st[p].child[c-'A'] == q){
					st[p].child[c-'A'] = clone;
					p = st[p].link;
				}
				st[q].link = st[curr].link = clone;
			}
		}
		last = curr;
	}
};

vector<int> z_function(const string & s){
    int n = s.size();
    vector<int> z(n);
    for(int i = 1, l = 0, r = 0; i < n; ++i){
        if(i <= r)
            z[i] = min(r - i + 1, z[i - l]);
        while(i + z[i] < n && s[z[i]] == s[i + z[i]])
            ++z[i];
        if(i + z[i] - 1 > r)
            l = i, r = i + z[i] - 1;
    }
    return z;
}

/*int main(){
	AhoCorasick ac;
	int patterns;
	cin >> patterns;
	for(int i = 0; i < patterns; ++i){
		string pattern;
		cin >> pattern;
		ac.add(pattern);
	}
	ac.build();
	string text;
	cin >> text;
	ac.match(text);
	return 0;
}*/

/*int main(){
	int t;
	string w;
	Trie T;
	while(cin >> t && t != -1){
		if(t == 0){
			cin >> w;
			T.InsertWord(w);
		}else if(t == 1){
			cin >> w;
			if(T.FindWord(w)) cout << "Found\n";
			else cout << "Not found\n";
		}else if(t == 2){
			T.printWords("");
		}
	}
}*/