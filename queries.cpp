#include <bits/stdc++.h>
using namespace std;

template<typename T>
struct SegmentTree{
	int N;
	vector<T> ST;

	//build from an array in O(n)
	SegmentTree(int N, vector<T> & arr): N(N){
		ST.resize(N << 1);
		for(int i = 0; i < N; ++i)
			ST[N + i] = arr[i];
		for(int i = N - 1; i > 0; --i)
			ST[i] = ST[i << 1] + ST[i << 1 | 1];
	}

	//single element update in i
	void update(int i, T value){
		ST[i += N] = value; //update the element accordingly
		while(i >>= 1)
			ST[i] = ST[i << 1] + ST[i << 1 | 1];
	}

	//single element update in [l, r]
	void update(int l, int r, T value){
		l += N, r += N;
		for(int i = l; i <= r; ++i)
			ST[i] = value;
		l >>= 1, r >>= 1;
		while(l >= 1){
			for(int i = r; i >= l; --i)
				ST[i] = ST[i << 1] + ST[i << 1 | 1];
			l >>= 1, r >>= 1;
		}
	}

	//range query, [l, r]
	T query(int l, int r){
		T res = 0;
		for(l += N, r += N; l <= r; l >>= 1, r >>= 1){
			if(l & 1) res += ST[l++];
			if(!(r & 1)) res += ST[r--];
		}
		return res;
	}
};

template<typename T>
struct SegmentTreeDin{
	SegmentTreeDin *left, *right;
	int l, r;
	T sum, lazy;
 
	SegmentTreeDin(int start, int end, vector<T> & arr): left(NULL), right(NULL), l(start), r(end), sum(0), lazy(0){
		if(l == r) sum = arr[l];
		else{
			int half = l + ((r - l) >> 1);
			left = new SegmentTreeDin(l, half, arr);
			right = new SegmentTreeDin(half+1, r, arr);
			sum = left->sum + right->sum;
		}
	}
 
	void propagate(T dif){
		sum += (r - l + 1) * dif;
		if(l != r){
			left->lazy += dif;
			right->lazy += dif;
		}
	}
 
	T sum_query(int start, int end){
		if(lazy != 0){
			propagate(lazy);
			lazy = 0;
		}
		if(end < l || r < start) return 0;
		if(start <= l && r <= end) return sum;
		else return left->sum_query(start, end) + right->sum_query(start, end);
	}
 
	void add_range(int start, int end, T dif){
		if(lazy != 0){
			propagate(lazy);
			lazy = 0;
		}
		if(end < l || r < start) return;
		if(start <= l && r <= end) propagate(dif);
		else{
			left->add_range(start, end, dif);
			right->add_range(start, end, dif);
			sum = left->sum + right->sum;
		}
	}

	void add_pos(int i, T sum){
		add_range(i, i, sum);
	}
};

template<typename T>
struct SegmentTreeEst{
	int size;
	vector<T> sum, lazy;
 
	void rec(int pos, int l, int r, vector<T> & arr){
		if(l == r) sum[pos] = arr[l];
		else{
			int half = l + ((r - l) >> 1);
			rec(2*pos+1, l, half, arr);
			rec(2*pos+2, half+1, r, arr);
			sum[pos] = sum[2*pos+1] + sum[2*pos+2];
		}
	}

	SegmentTreeEst(int n, vector<T> & arr): size(n){
		int h = ceil(log2(n));
		sum.resize((1 << (h + 1)) - 1);
		lazy.resize((1 << (h + 1)) - 1);
		rec(0, 0, n - 1, arr);
	}
 
	void propagate(int pos, int l, int r, T dif){
		sum[pos] += (r - l + 1) * dif;
		if(l != r){
			lazy[2*pos+1] += dif;
			lazy[2*pos+2] += dif;
		}
	}
 
	T sum_query_rec(int start, int end, int pos, int l, int r){
		if(lazy[pos] != 0){
			propagate(pos, l, r, lazy[pos]);
			lazy[pos] = 0;
		}
		if(end < l || r < start) return 0;
		if(start <= l && r <= end) return sum[pos];
		else{
			int half = l + ((r - l) >> 1);
			return sum_query_rec(start, end, 2*pos+1, l, half) + sum_query_rec(start, end, 2*pos+2, half+1, r);
		}
	}

	T sum_query(int start, int end){
		return sum_query_rec(start, end, 0, 0, size - 1);
	}
 
	void add_range_rec(int start, int end, int pos, int l, int r, T dif){
		if(lazy[pos] != 0){
			propagate(pos, l, r, lazy[pos]);
			lazy[pos] = 0;
		}
		if(end < l || r < start) return;
		if(start <= l && r <= end) propagate(pos, l, r, dif);
		else{
			int half = l + ((r - l) >> 1);
			add_range_rec(start, end, 2*pos+1, l, half, dif);
			add_range_rec(start, end, 2*pos+2, half+1, r, dif);
			sum[pos] = sum[2*pos+1] + sum[2*pos+2];
		}
	}

	void add_range(int start, int end, T dif){
		add_range_rec(start, end, 0, 0, size - 1, dif);
	}

	void add_pos(int i, T sum){
		add_range(i, i, sum);
	}
};

template<typename T>
struct StPer{
	StPer *left, *right;
	int l, r;
	T sum;
 
	StPer(int start, int end): left(NULL), right(NULL), l(start), r(end), sum(0){
		if(l != r){
			int half = l + ((r - l) >> 1);
			left = new StPer(l, half);
			right = new StPer(half+1, r);
		}
	}
	StPer(int start, int end, T val): left(NULL), right(NULL), l(start), r(end), sum(val){}
	StPer(int start, int end, StPer* left, StPer* right): left(left), right(right), l(start), r(end){
		sum = left->sum + right->sum;
	}
 
	T sum_query(int start, int end){
		if(end < l || r < start) return 0;
		if(start <= l && r <= end) return sum;
		else return left->sum_query(start, end) + right->sum_query(start, end);
	}
 
	StPer* update(int pos, T val){
		if(l == r) return new StPer(l, r, sum + val);
		int half = l + ((r - l) >> 1);
		if(pos <= half) return new StPer(l, r, left->update(pos, val), right);
		return new StPer(l, r, left, right->update(pos, val));
	}
};

template<typename T>
struct FenwickTree{
	int N;
	vector<T> bit;

	//build from array in O(n), indexed in 0
	FenwickTree(int N, vector<T> & arr): N(N){
		bit.resize(N);
		for(int i = 0; i < N; ++i){
			bit[i] += arr[i];
			if((i | (i + 1)) < N)
				bit[i | (i + 1)] += bit[i];
		}
	}

	//single element increment
	void update(int pos, T value){
		while(pos < N){
			bit[pos] += value;
			pos |= pos + 1;
		}
	}

	//range query, [0, r]
	T query(int r){
		T res = 0;
		while(r >= 0){
			res += bit[r];
			r = (r & (r + 1)) - 1;
		}
		return res;
	}

	//range query, [l, r]
	T query(int l, int r){
		return query(r) - query(l - 1);
	}
};

struct MOquery{
	int l, r, index, S;
	bool operator<(const MOquery & q) const{
		int c_o = l / S, c_q = q.l / S;
		if(c_o == c_q)
			return r < q.r;
		return c_o < c_q;
	}
};

template<typename T>
struct SQRT{
	int N, S;
	vector<T> A, B;

	SQRT(int N): N(N){
		this->S = sqrt(N + .0) + 1;
		A.assign(N, 0);
		B.assign(S, 0);
	}

	void build(vector<T> & arr){
		A = vector<int>(arr.begin(), arr.end());
		for(int i = 0; i < N; ++i) B[i / S] += A[i];
	}

	//single element update
	void update(int pos, T value){
		int k = pos / S;
		A[pos] = value;
		T res = 0;
		for(int i = k * S, end = min(N, (k + 1) * S) - 1; i <= end; ++i) res += A[i];
		B[k] = res;
	}

	//range query, [l, r]
	T query(int l, int r){
		T res = 0;
		int c_l = l / S, c_r = r / S;
		if(c_l == c_r){
			for(int i = l; i <= r; ++i) res += A[i];
		}else{
			for(int i = l, end = (c_l + 1) * S - 1; i <= end; ++i) res += A[i];
			for(int i = c_l + 1; i <= c_r - 1; ++i) res += B[i];
			for(int i = c_r * S; i <= r; ++i) res += A[i];
		}
		return res;
	}

	//range queries offline using MO's algorithm
	vector<T> MO(vector<MOquery> & queries){
		vector<T> ans(queries.size());
		sort(queries.begin(), queries.end());
		T current = 0;
		int prevL = 0, prevR = -1;
		int i, j;
		for(const MOquery & q : queries){
			for(i = prevL, j = min(prevR, q.l - 1); i <= j; ++i){
				//remove from the left
				current -= A[i];
			}
			for(i = prevL - 1; i >= q.l; --i){
				//add to the left
				current += A[i];
			}
			for(i = max(prevR + 1, q.l); i <= q.r; ++i){
				//add to the right
				current += A[i];
			}
			for(i = prevR; i >= q.r + 1; --i){
				//remove from the right
				current -= A[i];
			}
			prevL = q.l, prevR = q.r;
			ans[q.index] = current;
		}
		return ans;
	}
};

template<typename T>
struct AVLNode{
	AVLNode<T> *left, *right;
	short int height;
	int size;
	T value;

	AVLNode(T value = 0): left(NULL), right(NULL), value(value), height(1), size(1){}

	inline short int balance(){
		return (right ? right->height : 0) - (left ? left->height : 0);
	}

	AVLNode *maxLeftChild(){
		AVLNode *ret = this;
		while(ret->left) ret = ret->left;
		return ret;
	}
};

template<typename T>
struct AVLTree{
	AVLNode<T> *root;

	AVLTree(): root(NULL){}

	inline int nodeSize(AVLNode<T> *& pos){return pos ? pos->size: 0;}

	inline int nodeHeight(AVLNode<T> *& pos){return pos ? pos->height: 0;}

	inline void update(AVLNode<T> *& pos){
		if(!pos) return;
		pos->height = 1 + max(nodeHeight(pos->left), nodeHeight(pos->right));
		pos->size = 1 + nodeSize(pos->left) + nodeSize(pos->right);
	}

	int size(){return nodeSize(root);}

	void leftRotate(AVLNode<T> *& x){
		AVLNode<T> *y = x->right, *t = y->left;
		y->left = x, x->right = t;
		update(x), update(y);
		x = y;
	}

	void rightRotate(AVLNode<T> *& y){
		AVLNode<T> *x = y->left, *t = x->right;
		x->right = y, y->left = t;
		update(y), update(x);
		y = x;
	}

	void updateBalance(AVLNode<T> *& pos){
		if(!pos) return;
		short int bal = pos->balance();
		if(bal > 1){
			if(pos->right->balance() < 0) rightRotate(pos->right);
			leftRotate(pos);
		}else if(bal < -1){
			if(pos->left->balance() > 0) leftRotate(pos->left);
			rightRotate(pos);
		}
	}

	void insert(AVLNode<T> *&pos, T & value){
		if(pos){
			value < pos->value ? insert(pos->left, value) : insert(pos->right, value);
			update(pos), updateBalance(pos);
		}else{
			pos = new AVLNode<T>(value);
		}
	}

	AVLNode<T> *search(T & value){
		AVLNode<T> *pos = root;
		while(pos){
			if(value == pos->value) break;
			pos = (value < pos->value ? pos->left : pos->right);
		}
		return pos;
	}

	void erase(AVLNode<T> *&pos, T & value){
		if(!pos) return;
		if(value < pos->value) erase(pos->left, value);
		else if(value > pos->value) erase(pos->right, value);
		else{
			if(!pos->left) pos = pos->right;
			else if(!pos->right) pos = pos->left;
			else{
				pos->value = pos->right->maxLeftChild()->value;
				erase(pos->right, pos->value);
			}
		}
		update(pos), updateBalance(pos);
	}

	void insert(T value){insert(root, value);}

	void erase(T value){erase(root, value);}

	void updateVal(T old, T New){
		if(search(old))
			erase(old), insert(New);
	}

	T kth(int i){
		assert(0 <= i && i < nodeSize(root));
		AVLNode<T> *pos = root;
		while(i != nodeSize(pos->left)){
			if(i < nodeSize(pos->left)){
				pos = pos->left;
			}else{
				i -= nodeSize(pos->left) + 1;
				pos = pos->right;
			}
		}
		return pos->value;
	}

	int lessThan(T & x){
		int ans = 0;
		AVLNode<T> *pos = root;
		while(pos){
			if(x > pos->value){
				ans += nodeSize(pos->left) + 1;
				pos = pos->right;
			}else{
				pos = pos->left;
			}
		}
		return ans;
	}

	int lessThanOrEqual(T & x){
		int ans = 0;
		AVLNode<T> *pos = root;
		while(pos){
			if(x < pos->value){
				pos = pos->left;
			}else{
				ans += nodeSize(pos->left) + 1;
				pos = pos->right;
			}
		}
		return ans;
	}

	int greaterThan(T & x){
		int ans = 0;
		AVLNode<T> *pos = root;
		while(pos){
			if(x < pos->value){
				ans += nodeSize(pos->right) + 1;
				pos = pos->left;
			}else{
				pos = pos->right;
			}
		}
		return ans;
	}

	int greaterThanOrEqual(T & x){
		int ans = 0;
		AVLNode<T> *pos = root;
		while(pos){
			if(x > pos->value){
				pos = pos->right;
			}else{
				ans += nodeSize(pos->right) + 1;
				pos = pos->left;
			}
		}
		return ans;
	}

	int equalTo(T & x){
		return lessThanOrEqual(x) - lessThan(x);
	}

	void build(AVLNode<T> *& pos, vector<T> & arr, int i, int j){
		if(i > j) return;
		int m = i + ((j - i) >> 1);
		pos = new AVLNode<T>(arr[m]);
		build(pos->left, arr, i, m - 1);
		build(pos->right, arr, m + 1, j);
		update(pos);
	}

	void build(vector<T> & arr){
		build(root, arr, 0, (int)arr.size() - 1);
	}

	void output(AVLNode<T> *pos, vector<T> & arr, int & i){
		if(pos){
			output(pos->left, arr, i);
			arr[++i] = pos->value;
			output(pos->right, arr, i);
		}
	}

	void output(vector<T> & arr){
		int i = -1;
		output(root, arr, i);
	}
};

template<typename T>
struct TreapNode{
	TreapNode<T> *left, *right;
	T value;
	int key, size;

	//fields for queries
	bool rev;
	T sum, add;

	TreapNode(T value = 0): value(value), key(rand()), size(1), left(NULL), right(NULL), sum(value), add(0), rev(false){}
};

template<typename T>
struct Treap{
	TreapNode<T> *root;

	Treap(): root(NULL) {}

	inline int nodeSize(TreapNode<T>* t){return t ? t->size: 0;}

	inline T nodeSum(TreapNode<T>* t){return t ? t->sum : 0;}

	inline void update(TreapNode<T>* &t){
		if(!t) return;
		t->size = 1 + nodeSize(t->left) + nodeSize(t->right);
		t->sum = t->value; //reset node fields
		push(t->left), push(t->right); //push changes to child nodes
		t->sum = t->value + nodeSum(t->left) + nodeSum(t->right); //combine(left,t,t), combine(t,right,t)
	}

	int size(){return nodeSize(root);}

	void merge(TreapNode<T>* &t, TreapNode<T>* t1, TreapNode<T>* t2){
		if(!t1) t = t2;
		else if(!t2) t = t1;
		else if(t1->key > t2->key)
			merge(t1->right, t1->right, t2), t = t1;
		else
			merge(t2->left, t1, t2->left), t = t2;
		update(t);
	}

	void split(TreapNode<T>* t, T & x, TreapNode<T>* &t1, TreapNode<T>* &t2){
		if(!t)
			return void(t1 = t2 = NULL);
		if(x < t->value)
			split(t->left, x, t1, t->left), t2 = t;
		else
			split(t->right, x, t->right, t2), t1 = t;
		update(t);
	}

	void insert(TreapNode<T>* &t, TreapNode<T>* x){
		if(!t) t = x;
		else if(x->key > t->key)
			split(t, x->value, x->left, x->right), t = x;
		else
			insert(x->value < t->value ? t->left : t->right, x);
		update(t);
	}

	TreapNode<T>* search(T & x){
		TreapNode<T> *t = root;
		while(t){
			if(x == t->value) break;
			t = (x < t->value ? t->left : t->right);
		}
		return t;
	}

	void erase(TreapNode<T>* &t, T & x){
		if(!t) return;
		if(t->value == x)
			merge(t, t->left, t->right);
		else
			erase(x < t->value ? t->left : t->right, x);
		update(t);
	}

	void insert(T & x){insert(root, new TreapNode<T>(x));}

	void erase(T & x){erase(root, x);}

	void updateVal(T & old, T & New){
		if(search(old))
			erase(old), insert(New);
	}

	T kth(int i){
		assert(0 <= i && i < nodeSize(root));
		TreapNode<T> *t = root;
		while(i != nodeSize(t->left)){
			if(i < nodeSize(t->left)){
				t = t->left;
			}else{
				i -= nodeSize(t->left) + 1;
				t = t->right;
			}
		}
		return t->value;
	}

	int lessThan(T & x){
		int ans = 0;
		TreapNode<T> *t = root;
		while(t){
			if(x > t->value){
				ans += nodeSize(t->left) + 1;
				t = t->right;
			}else{
				t = t->left;
			}
		}
		return ans;
	}

	//OPERATIONS FOR IMPLICIT TREAP
	inline void push(TreapNode<T>* t){
		if(!t) return;
		//add in range example
		if(t->add){
			t->value += t->add;
			t->sum += t->add * nodeSize(t);
			if(t->left) t->left->add += t->add;
			if(t->right) t->right->add += t->add;
			t->add = 0;
		}
		//reverse range example
		if(t->rev){
			swap(t->left, t->right);
			if(t->left) t->left->rev ^= true;
			if(t->right) t->right->rev ^= true;
			t->rev = false;
		}
	}

	void split2(TreapNode<T>* t, int i, TreapNode<T>* &t1, TreapNode<T>* &t2){
		if(!t)
			return void(t1 = t2 = NULL);
		push(t);
		int curr = nodeSize(t->left);
		if(i <= curr)
			split2(t->left, i, t1, t->left), t2 = t;
		else
			split2(t->right, i - curr - 1, t->right, t2), t1 = t;
		update(t);
	}

	inline int aleatorio(){
		return (rand() << 15) + rand();
	}

	void merge2(TreapNode<T>* &t, TreapNode<T>* t1, TreapNode<T>* t2){
		push(t1), push(t2);
		if(!t1) t = t2;
		else if(!t2) t = t1;
		else if(aleatorio() % (nodeSize(t1) + nodeSize(t2)) < nodeSize(t1))
			merge2(t1->right, t1->right, t2), t = t1;
		else
			merge2(t2->left, t1, t2->left), t = t2;
		update(t);
	}

	//insert the element "x" at position "i"
	void insert_at(T & x, int i){
		if(i > nodeSize(root)) return;
		TreapNode<T> *t1 = NULL, *t2 = NULL;
		split2(root, i, t1, t2);
		merge2(root, t1, new TreapNode<T>(x));
		merge2(root, root, t2);
	}

	//delete element at position "i"
	void erase_at(int i){
		if(i >= nodeSize(root)) return;
		TreapNode<T> *t1 = NULL, *t2 = NULL, *t3 = NULL;
		split2(root, i, t1, t2);
		split2(t2, 1, t2, t3);
		merge2(root, t1, t3);
	}

	void update_at(TreapNode<T>* t, T & x, int i){
		push(t);
		assert(0 <= i && i < nodeSize(t));
		int curr = nodeSize(t->left);
		if(i == curr)
			t->value = x;
		else if(i < curr)
			update_at(t->left, x, i);
		else
			update_at(t->right, x, i - curr - 1);
		update(t);
	}

	T nth(TreapNode<T>* t, int i){
		push(t);
		assert(0 <= i && i < nodeSize(t));
		int curr = nodeSize(t->left);
		if(i == curr)
			return t->value;
		else if(i < curr)
			return nth(t->left, i);
		else
			return nth(t->right, i - curr - 1);
	}

	//update value of element at position "i" with "x"
	void update_at(T & x, int i){update_at(root, x, i);}

	//ith element
	T nth(int i){return nth(root, i);}

	//add "val" in [l, r]
	void add_update(T & val, int l, int r){
		TreapNode<T> *t1 = NULL, *t2 = NULL, *t3 = NULL;
		split2(root, l, t1, t2);
		split2(t2, r - l + 1, t2, t3);
		t2->add += val;
		merge2(root, t1, t2);
		merge2(root, root, t3);
	}

	//reverse [l, r]
	void reverse_update(int l, int r){
		TreapNode<T> *t1 = NULL, *t2 = NULL, *t3 = NULL;
		split2(root, l, t1, t2);
		split2(t2, r - l + 1, t2, t3);
		t2->rev ^= true;
		merge2(root, t1, t2);
		merge2(root, root, t3);
	}

	//rotate [l, r] k times to the right
	void rotate_update(int k, int l, int r){
		TreapNode<T> *t1 = NULL, *t2 = NULL, *t3 = NULL, *t4 = NULL;
		split2(root, l, t1, t2);
		split2(t2, r - l + 1, t2, t3);
		k %= nodeSize(t2);
		split2(t2, nodeSize(t2) - k, t2, t4);
		merge2(root, t1, t4);
		merge2(root, root, t2);
		merge2(root, root, t3);
	}

	//sum query in [l, r]
	T sum_query(int l, int r){
		TreapNode<T> *t1 = NULL, *t2 = NULL, *t3 = NULL;
		split2(root, l, t1, t2);
		split2(t2, r - l + 1, t2, t3);
		T ans = nodeSum(t2);
		merge2(root, t1, t2);
		merge2(root, root, t3);
		return ans;
	}

	void inorder(TreapNode<T>* t){
		if(!t) return;
		push(t);
		inorder(t->left);
		cout << t->value << " ";
		inorder(t->right);
	}

	void inorder(){inorder(root);}
};

template<typename T>
struct SparseTable{
	vector<vector<T>> ST;
	vector<int> logs;
	int K, N;

	SparseTable(vector<T> & arr){
		N = arr.size();
		K = log2(N) + 2;
		ST.assign(K + 1, vector<T>(N));
		logs.assign(N + 1, 0);
		for(int i = 2; i <= N; ++i)
			logs[i] = logs[i >> 1] + 1;
		for(int i = 0; i < N; ++i)
			ST[0][i] = arr[i];
		for(int j = 1; j <= K; ++j)
			for(int i = 0; i + (1 << j) <= N; ++i)
				ST[j][i] = min(ST[j - 1][i], ST[j - 1][i + (1 << (j - 1))]); //put the function accordingly
	}

	T sum(int l, int r){ //non-idempotent functions
		T ans = 0;
		for(int j = K; j >= 0; --j){
			if((1 << j) <= r - l + 1){
				ans += ST[j][l];
				l += 1 << j;
			}
		}
		return ans;
	}

	T minimal(int l, int r){ //idempotent functions
		int j = logs[r - l + 1];
		return min(ST[j][l], ST[j][r - (1 << j) + 1]);
	}
};

//build on O(n log n), queries in O(1) for any operation
template<typename T>
struct DisjointSparseTable{
	vector<vector<T>> left, right;
	int K, N;

	DisjointSparseTable(vector<T> & arr){
		N = arr.size();
		K = log2(N) + 2;
		left.assign(K + 1, vector<T>(N));
		right.assign(K + 1, vector<T>(N));
		for(int j = 0; (1 << j) <= N; ++j){
			int mask = (1 << j) - 1;
			T acum = 0; //neutral element of your operation
			for(int i = 0; i < N; ++i){
				acum += arr[i]; //your operation
				left[j][i] = acum;
				if((i & mask) == mask) acum = 0; //neutral element of your operation
			}
			acum = 0; //neutral element of your operation
			for(int i = N-1; i >= 0; --i){
				acum += arr[i]; //your operation
				right[j][i] = acum;
				if((i & mask) == 0) acum = 0; //neutral element of your operation
			}
		}
	}

	T query(int l, int r){
		if(l == r) return left[0][l];
		int i = 31 - __builtin_clz(l^r);
		return left[i][r] + right[i][l]; //your operation
	}
};

struct WaveletTree{
	int lo, hi;
	WaveletTree *left, *right;
	vector<int> freq;
	vector<int> pref; //just use this if you want sums

	//queries indexed in base 1, complexity for all queries: O(log(max_element))
	//build from [from, to) with non-negative values in range [x, y]
	//you can use vector iterators or array pointers
	WaveletTree(vector<int>::iterator from, vector<int>::iterator to, int x, int y): lo(x), hi(y){
		if(from >= to) return;
		int m = (lo + hi) / 2;
		auto f = [m](int x){return x <= m;};
		freq.reserve(to - from + 1);
		freq.push_back(0);
		pref.reserve(to - from + 1);
		pref.push_back(0);
		for(auto it = from; it != to; ++it){
			freq.push_back(freq.back() + f(*it));
			pref.push_back(pref.back() + *it);
		}
		if(hi != lo){
			auto pivot = stable_partition(from, to, f);
			left = new WaveletTree(from, pivot, lo, m);
			right = new WaveletTree(pivot, to, m + 1, hi);
		}
	}

	//kth element in [l, r]
	int kth(int l, int r, int k){
		if(l > r) return 0;
		if(lo == hi) return lo;
		int lb = freq[l - 1], rb = freq[r];
		int inLeft = rb - lb;
		if(k <= inLeft) return left->kth(lb + 1, rb, k);
		else return right->kth(l - lb, r - rb, k - inLeft);
	}

	//number of elements less than or equal to k in [l, r]
	int lessThanOrEqual(int l, int r, int k){
		if(l > r || k < lo) return 0;
		if(hi <= k) return r - l + 1;
		int lb = freq[l - 1], rb = freq[r];
		return left->lessThanOrEqual(lb + 1, rb, k) + right->lessThanOrEqual(l - lb, r - rb, k);
	}

	//number of elements equal to k in [l, r]
	int equalTo(int l, int r, int k){
		if(l > r || k < lo || k > hi) return 0;
		if(lo == hi) return r - l + 1;
		int lb = freq[l - 1], rb = freq[r];
		int m = (lo + hi) / 2;
		if(k <= m) return left->equalTo(lb + 1, rb, k);
		else return right->equalTo(l - lb, r - rb, k);
	}

	//sum of elements less than or equal to k in [l, r]
	int sum(int l, int r, int k){
		if(l > r || k < lo) return 0;
		if(hi <= k) return pref[r] - pref[l - 1];
		int lb = freq[l - 1], rb = freq[r];
		return left->sum(lb + 1, rb, k) + right->sum(l - lb, r - rb, k);
	}
};

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

template<typename T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;

int main(){
	int t, n, m;
	ordered_set<int> conj;
	while(cin >> t && t != -1){
		cin >> n;
		if(t == 0){ //insert
			conj.insert(n);
		}else if(t == 1){ //search
			if(conj.find(n) != conj.end()) cout << "Found\n";
			else cout << "Not found\n";
		}else if(t == 2){ //delete
			conj.erase(n);
		}else if(t == 3){ //update
			cin >> m;
			if(conj.find(n) != conj.end()){
				conj.erase(n);
				conj.insert(n);
			}
		}else if(t == 4){ //lower bound
			cout << conj.order_of_key(n) << "\n";
		}else if(t == 5){ //get nth element
			auto pos = conj.find_by_order(n);
			if(pos != conj.end()) cout << *pos << "\n";
			else cout << "-1\n";
		}
	}
	return 0;
}

struct HeavyLight{
	int n;
	vector<vector<int>> adj;
	vector<int> parent, level, size, heavy, head, pos, ipos;
	int cur_pos;
	SegmentTree<int> * st;

	HeavyLight(int n, SegmentTree<int> * st): n(n), st(st){
		adj.resize(n), ipos.resize(n);
		parent.resize(n), level.resize(n), size.resize(n);
		heavy.resize(n, -1), head.resize(n), pos.resize(n);
	}

	void dfs(int u){
		size[u] = 1;
		int mx = 0;
		for(int v : adj[u]){
			if(v != parent[u]){
				parent[v] = u;
				level[v] = level[u] + 1;
				dfs(v);
				if(size[v] > mx){
					mx = size[v];
					heavy[u] = v;
				}
				size[u] += size[v];
			}
		}
	}

	void build(int u, int h){
		head[u] = h;
		pos[u] = cur_pos;
		ipos[cur_pos++] = u;
		if(heavy[u] != -1) build(heavy[u], h);
		for(int v : adj[u])
			if(v != parent[u] && v != heavy[u])
				build(v, v);
	}

	void init(int root = 0){
		cur_pos = 0;
		dfs(root);
		build(root, root);
	}

	int query(int a, int b){
		int mx = 0;
		while(head[a] != head[b]){
			if(level[head[a]] > level[head[b]]) swap(a, b);
			mx = max(mx, st->query(pos[head[b]], pos[b]));
			b = parent[head[b]];
		}
		if(level[a] > level[b]) swap(a, b);
		// if(pos[a] + 1 <= pos[b]) for values in edges
		mx = max(mx, st->query(pos[a], pos[b]));
		//LCA at a
		return mx;
	}

	int kth_ancestor(int u, int k){
		while(pos[u] - pos[head[u]] < k){
			k -= pos[u] - pos[head[u]] + 1;
			u = parent[head[u]];
		}
		return ipos[pos[u] - k];
	}
};

/*int main(){
	int n, t, pos, value, l, r;
	cin >> n;
	vector<int> a(n);
	for(int i = 0; i < n; i++) cin >> a[i];
	SegmentTreeDin<int> *st = new SegmentTreeDin<int>(0, n-1, a);
	while(cin >> t && t != -1){
		if(t == 1){ //update single element
			cin >> pos >> value;
			st->add_pos(pos, value);
		}else if(t == 2){ //query
			cin >> l >> r;
			cout << st->sum_query(l, r) << "\n";
		}else if(t == 3){ //update range with element
			cin >> l >> r >> value;
			st->add_range(l, r, value);
		}
	}
	return 0;
}*/

/*int main(){
	int n, t, pos, value, l, r;
	cin >> n;
	vector<int> a(n);
	for(int i = 0; i < n; i++) cin >> a[i];
	SegmentTreeEst<int> *st = new SegmentTreeEst<int>(n, a);
	while(cin >> t && t != -1){
		if(t == 1){ //update single element
			cin >> pos >> value;
			st->add_pos(pos, value);
		}else if(t == 2){ //query
			cin >> l >> r;
			cout << st->sum_query(l, r) << "\n";
		}else if(t == 3){ //update range with element
			cin >> l >> r >> value;
			st->add_range(l, r, value);
		}
	}
	return 0;
}*/

/*int main(){
	int n, q, l, r;
	cin >> n;
	vector<int> a(n);
	for(int i = 0; i < n; i++) cin >> a[i];
	SQRT<int> *s = new SQRT<int>(n);
	s->build(a);
	cin >> q;
	vector<MOquery> queries(q);
	for(int i = 0; i < q; ++i){
		cin >> l >> r;
		queries[i] = {l, r, i, s->S};
	}
	vector<int> ans = s->MO(queries);
	for(int & x : ans){
		cout << x << "\n";
	}
}
*/

/*int main(){
	int t, n, m;
	AVLTree<int> *avl = new AVLTree<int>;
	while(cin >> t && t != -1){
		cin >> n;
		if(t == 0){ //insert
			avl->insert(n);
		}else if(t == 1){ //search
			AVLNode<int> *pos = avl->search(n);
			if(pos) cout << "Found\n";
			else cout << "Not found\n";
		}else if(t == 2){ //delete
			avl->erase(n);
		}else if(t == 3){ //update
			cin >> m;
			avl->updateVal(n, m);
		}else if(t == 4){ //lessThanOrEqual
			cout << avl->lessThan(n) << " " << avl->lessThanOrEqual(n) << " " << avl->greaterThan(n) << " " << avl->greaterThanOrEqual(n) << " " << avl->equalTo(n) << "\n";
		}else if(t == 5){ //get nth element
			cout << avl->kth(n) << "\n";
		}
	}
	return 0;
}*/

/*int main(){
	srand(time(NULL));
	int t, n, m;
	Treap<int> *T = new Treap<int>();
	while(cin >> t && t != -1){
		cin >> n;
		if(t == 0){ //insert
			T->insert(n);
		}else if(t == 1){ //search
			TreapNode<int> *pos = T->search(n);
			if(pos) cout << "Found\n";
			else cout << "Not found\n";
		}else if(t == 2){ //delete
			T->erase(n);
		}else if(t == 3){ //update
			cin >> m;
			T->updateVal(n, m);
		}else if(t == 4){ //lessThan
			cout << T->lessThan(n) << "\n";
		}else if(t == 5){ //get nth element
			cout << T->kth(n) << "\n";
		}
	}
	return 0;
}*/

/*int main(){
	srand(time(NULL));
	int t, n, i, l, r, val, k;
	Treap<int> *T = new Treap<int>();
	while(cin >> t && t != -1){
		if(t == 0){ //insert n at i
			cin >> n >> i;
			T->insert_at(n, i);
		}else if(t == 1){ //delete at i
			cin >> i;
			T->erase_at(i);
		}else if(t == 2){ //update value at i with n
			cin >> n >> i;
			T->update_at(n, i);
		}else if(t == 3){ //get nth element
			cin >> n;
			cout << T->nth(n) << "\n";
		}else if(t == 4){ //add "val" to [l, r]
			cin >> l >> r >> val;
			T->add_update(val, l, r);
		}else if(t == 5){ //get sum in [l, r]
			cin >> l >> r;
			cout << T->sum_query(l, r) << "\n";
		}else if(t == 6){ //reverse [l, r]
			cin >> l >> r;
			T->reverse_update(l, r);
		}else if(t == 7){ //rotate [l, r] k times to the right
			cin >> l >> r >> k;
			T->rotate_update(k, l, r);
		}else if(t == 8){ //inorder trasversal
			T->inorder(); cout << "\n";
		}
	}
	return 0;
}*/

/*int main(){
	int n, l, r;
	cin >> n;
	vector<int> arr(n);
	for(int i = 0; i < n; ++i)
		cin >> arr[i];
	SparseTable<int> table(arr);
	while(cin >> l && l != -1){
		cin >> r;
		cout << table.minimal(l, r) << "\n";
	}
	return 0;
}*/

/*int main(){
	int n, l, r;
	cin >> n;
	vector<int> arr(n);
	for(int i = 0; i < n; ++i)
		cin >> arr[i];
	DisjointSparseTable<int> table(arr);
	while(cin >> l && l != -1){
		cin >> r;
		cout << table.query(l, r) << "\n";
	}
	return 0;
}*/

/*int main(){
	int n, t, l, r, k;
	cin >> n;
	vector<int> arr(n);
	for(int i = 0; i < n; ++i){
		cin >> arr[i];
	}
	WaveletTree w(begin(arr), end(arr), *min_element(begin(arr), end(arr)), *max_element(begin(arr), end(arr)));
	while(cin >> t && t != -1){
		cin >> l >> r >> k;
		if(t == 0){ //kth smallest
			cout << w.kth(l, r, k) << "\n";
		}else if(t == 1){ //less than or equal to k
			cout << w.lessThanOrEqual(l, r, k) << "\n";
		}else if(t == 2){ //equal to k
			cout << w.equalTo(l, r, k) << "\n";
		}else if(t == 3){ //sum of elements less than or equal to k
			cout << w.sum(l, r, k) << "\n";
		}
	}
}*/