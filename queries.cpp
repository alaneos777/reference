#include <bits/stdc++.h>
using namespace std;

template<typename T>
struct SegmentTree{
	int N;
	vector<T> ST;

	SegmentTree(int N){
		this->N = N;
		ST.assign(N << 1, 0);
	}

	void build(vector<T> & arr){
		for(int i = 0; i < N; ++i)
			ST[N + i] = arr[i];
		for(int i = N - 1; i > 0; --i)
			ST[i] = ST[i << 1] + ST[i << 1 | 1];
	}

	//single element update in pos
	void update(int pos, T value){
		ST[pos += N] = value;
		while(pos >>= 1)
			ST[pos] = ST[pos << 1] + ST[pos << 1 | 1];
	}

	//single element update in [l, r]
	void update(int l, int r, T value){
		for(int i = l; i <= r; ++i)
			ST[N + i] = value;
		for(int i = r; i > l; --i)
			ST[i] = ST[i << 1] + ST[i << 1 | 1];
	}

	//range query, [l, r]
	T query(int l, int r){
		++r;
		T res = 0;
		for(l += N, r += N; l < r; l >>= 1, r >>= 1) {
			if(l & 1) res += ST[l++];
			if(r & 1) res += ST[--r];
		}
		return res;
	}
};

template<typename T>
struct FenwickTree{
	int N;
	vector<T> bit;

	FenwickTree(int N){
		this->N = N;
		bit.assign(N, 0);
	}

	void build(vector<T> & arr){
		for(int i = 0; i < arr.size(); ++i){
			update(i, arr[i]);
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

template<typename T>
struct SQRT{
	int N, S;
	vector<T> A, B;

	SQRT(int N){
		this->N = N;
		this->S = sqrt(N + .0) + 1;
		A.assign(N, 0);
		B.assign(S, 0);
	}

	void build(vector<T> & arr){
		A = vector<int>(arr.begin(), arr.end());
		for(int i = 0; i < N; ++i){
			B[i / S] += A[i];
		}
	}

	//single element update
	void update(int pos, T value){
		int k = pos / S;
		A[pos] = value;
		T res = 0;
		for(int i = k * S, end = min(N, (k + 1) * S) - 1; i <= end; ++i){
			res += A[i];
		}
		B[k] = res;
	}

	//range query, [l, r]
	T query(int l, int r){
		T res = 0;
		int c_l = l / S, c_r = r / S;
		if(c_l == c_r){
			for(int i = l; i <= r; ++i)
				res += A[i];
		}else{
			for(int i = l, end = (c_l + 1) * S - 1; i <= end; ++i){
				res += A[i];
			}
			for(int i = c_l + 1; i <= c_r - 1; ++i){
				res += B[i];
			}
			for(int i = c_r * S; i <= r; ++i){
				res += A[i];
			}
		}
		return res;
	}
};

template<typename T>
struct AVLNode
{
	AVLNode<T> *left;
	AVLNode<T> *right;
	short int height;
	int size;
	T value;

	AVLNode(T value){
		left = right = NULL;
		this->value = value;
		height = 1, size = 0;
	}

	short int balance(){
		return (right ? right->height : 0) - (left ? left->height : 0);
	}

	void updateHeight(){
		height = 1 + max(left ? left->height : 0, right ? right->height : 0);
	}

	AVLNode *maxLeftChild(){
		AVLNode *ret = this;
		while(ret->left) ret = ret->left;
		return ret;
	}

	AVLNode *maxRightChild(){
		AVLNode *ret = this;
		while(ret->right) ret = ret->right;
		return ret;
	}
};

template<typename T>
struct AVLTree
{
	AVLNode<T> *root;
	int size;

	AVLTree(){
		root = NULL;
		size = 0;
	}

	int nodeSize(AVLNode<T> *& pos){
		return pos ? pos->size + 1: 0;
	}

	void leftRotate(AVLNode<T> *& x){
		AVLNode<T> *y = x->right;
		AVLNode<T> *t = y->left;
		y->left = x;
		x->right = t;
		x->updateHeight();
		y->updateHeight();
		int size = nodeSize(t);
		x->size = x->size - (y->size + 1) + size;
		y->size = y->size - size + (x->size + 1);
		x = y;
	}

	void rightRotate(AVLNode<T> *& y){
		AVLNode<T> *x = y->left;
		AVLNode<T> *t = x->right;
		x->right = y;
		y->left = t;
		y->updateHeight();
		x->updateHeight();
		int size = nodeSize(t);
		y->size = y->size - (x->size + 1) + size;
		x->size = x->size - size + (y->size + 1);
		y = x;
	}

	void updateBalance(AVLNode<T> *& pos){
		short int bal = pos->balance();
		if(bal > 1){
			if(pos->right->balance() < 0){
				rightRotate(pos->right);
				leftRotate(pos);
			}else{
				leftRotate(pos);
			}
		}else if(bal < -1){
			if(pos->left->balance() > 0){
				leftRotate(pos->left);
				rightRotate(pos);
			}else{
				rightRotate(pos);
			}
		}
	}

	void build(AVLNode<T> *& pos, const vector<T> & arr, int i, int j){
		if(i > j) return;
		int m = i + ((j - i) >> 1);
		pos = new AVLNode<T>(arr[m]);
		build(pos->left, arr, i, m - 1);
		build(pos->right, arr, m + 1, j);
		pos->size = j - i;
		pos->updateHeight();
	}

	void build(const vector<T> & arr){
		size = arr.size();
		build(root, arr, 0, size - 1);
	}

	void insert(AVLNode<T> *&pos, const T & value){
		if(pos){
			value < pos->value ? insert(pos->left, value) : insert(pos->right, value);
			++pos->size;
			pos->updateHeight();
			updateBalance(pos);
		}else{
			pos = new AVLNode<T>(value);
			++size;
		}
	}

	void insert(T value){
		insert(root, value);
	}

	AVLNode<T> *search(const T & value){
		AVLNode<T> *pos = root;
		while(pos){
			if(value == pos->value) break;
			pos = (value < pos->value ? pos->left : pos->right);
		}
		return pos;
	}

	bool erase(AVLNode<T> *&pos, const T & value){
		AVLNode<T> *tmp, *next;
		if(!pos) return false;
		bool success = false;
		if(value < pos->value){
			success = erase(pos->left, value);
			if(success) --pos->size;
		}else if(value > pos->value){
			success = erase(pos->right, value);
			if(success) --pos->size;
		}else{
			success = true;
			if(!pos->left){
				pos = pos->right;
			}else if(!pos->right){
				pos = pos->left;
			}else{
				next = pos->right->maxLeftChild(); //pos->left->maxRightChild();
				pos->value = next->value;
				erase(pos->right, pos->value);     //erase(pos->left, pos->value);
				--pos->size;
			}
		}
		if(pos && success){
			pos->updateHeight();
			updateBalance(pos);
		}
		return success;
	}

	void erase(T value){
		if(erase(root, value)) --size;
	}

	void update(T old, T New){
		erase(old);
		insert(New);
	}

	int lessThan(const T & x){
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

	int lessThanOrEqual(const T & x){
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

	int greaterThan(const T & x){
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

	int greaterThanOrEqual(const T & x){
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

	int equalTo(const T & x){
		return lessThanOrEqual(x) - lessThan(x);
	}

	T index(int i){
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

	void inorden(AVLNode<T> *pos){
		if(pos){
			inorden(pos->left);
			cout << pos->value << " " << pos->height << " " << pos->balance() << " " << pos->size << "\n";
			inorden(pos->right);
		}
	}
};

int main(){
	int n, q, t, pos, value, l, r;
	cin >> n;
	vector<int> a(n);
	for(int i = 0; i < n; i++) cin >> a[i];
	SegmentTree<int> *st = new SegmentTree<int>(n);
	st->build(a);
	cin >> q;
	while(q--){
		cin >> t;
		if(t == 1){ //update single element
			cin >> pos >> value;
			st->update(pos, value);
		}else if(t == 2){ //query
			cin >> l >> r;
			cout << st->query(l, r) << "\n";
		}else if(t == 3){ //update range with element
			cin >> l >> r >> value;
			st->update(l, r, value);
		}
	}
	return 0;
}

/*int main(){
	int q, t, n, m;
	AVLTree<int> *avl = new AVLTree<int>;
	cin >> q;
	while(q--){
		cin >> t;
		cin >> n;
		if(t == 0){ //insert
			avl->insert(n);
			cout << "Inorden:\n"; avl->inorden(avl->root); cout << "Size: " << avl->size << "\n\n";
		}else if(t == 1){ //search
			AVLNode<int> *pos = avl->search(n);
			if(pos) cout << "Found\n";
			else cout << "Not found\n";
		}else if(t == 2){ //delete
			avl->erase(n);
			cout << "Inorden:\n"; avl->inorden(avl->root); cout << "Size: " << avl->size << "\n\n";
		}else if(t == 3){ //update
			cin >> m;
			avl->update(n, m);
		}else if(t == 4){ //lessThanOrEqual
			cout << avl->lessThan(n) << " " << avl->lessThanOrEqual(n) << " " << avl->greaterThan(n) << " " << avl->greaterThanOrEqual(n) << " " << avl->equalTo(n) << "\n";
		}else if(t == 5){ //get nth element
			cout << avl->index(n) << "\n";
		}
	}
	return 0;
}*/