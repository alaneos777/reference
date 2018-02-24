#include <bits/stdc++.h>
using namespace std;

//Selection Sort
void selectionSort(vector<int> & arr, int a, int b){
	for(int i = a; i < b; i++){
		int minPos = i;
		for(int j = i + 1; j <= b; j++){
			if(arr[j] < arr[minPos]) minPos = j;
		}
		swap(arr[i], arr[minPos]);
	}
}

//Insertion Sort
void insertionSort(vector<int> & arr, int a, int b){
	for(int i = a + 1; i <= b; i++){
		int v = arr[i];
		int j = i;
		while(j > a && arr[j - 1] > v){
			arr[j] = arr[j - 1];
			j--;
		}
		arr[j] = v;
	}
}

//Bubble Sort
void bubbleSort(vector<int> & arr, int a, int b){
	bool change = true;
	int s = 0;
	while(change){
		change = false;
		for(int p = a; p < b - s; p++){
			if(arr[p] > arr[p + 1]){
				swap(arr[p], arr[p + 1]);
				change = true;
			}
		}
		s++;
	}
}

//Shell sort
void shellSort(vector<int> & arr, int a, int b){
	vector<int> gaps = {701, 301, 132, 57, 23, 10, 4, 1};
	for(int gap: gaps){
		for(int i = a + gap; i <= b; i++){
			int v = arr[i];
			int j = i;
			while(j >= a + gap && arr[j - gap] > v){
				arr[j] = arr[j - gap];
				j -= gap;
			}
			arr[j] = v;
		}
	}
}

//Merge Sort
void merge(vector<int> & arr, int a, int b, int c){
	int i = a, j = b + 1, p = 0;
	vector<int> tmp(c - a + 1);
	while(i <= b && j <= c)
		if(arr[i] <= arr[j]) tmp[p++] = arr[i++];
		else tmp[p++] = arr[j++];
	while(i <= b) tmp[p++] = arr[i++];
	while(j <= c) tmp[p++] = arr[j++];
	for(p = a; p <= c; p++){
		arr[p] = tmp[p - a];
	}
}

void mergeSort(vector<int> & arr, int a, int b){
	if(a < b){
		int m = a + ((b - a) >> 1);
		mergeSort(arr, a, m);
		mergeSort(arr, m + 1, b);
		merge(arr, a, m, b);
	}
}

//Quick Sort
int partition(vector<int> & arr, int a, int b){
	int p = b;
	a--;
	int pivot = arr[p];
	while(a < b){
		while(arr[++a] < pivot);
		while(b > 0 && arr[--b] > pivot);
		if(a < b) swap(arr[a], arr[b]);
	}
	swap(arr[a], arr[p]);
	return a;
}

void quickSort(vector<int> & arr, int a, int b){
	if(a < b){
		int m = partition(arr, a, b);
		quickSort(arr, a, m - 1);
		quickSort(arr, m + 1, b);
	}
}

//Heap sort
int left(int i, int a){
	return ((i - a) << 1) + 1 + a;
}
int right(int i, int a){
	return ((i - a) << 1) + 2 + a;
}
int parent(int i, int a){
	return ((i - a - 1) >> 1) + a;
}

void heapify(vector<int> & arr, int i, int a, int b){
	int l = left(i, a), r = right(i, a), largest = i;
	if(a <= l && l <= b && arr[l] > arr[largest]) largest = l;
	if(a <= r && r <= b && arr[r] > arr[largest]) largest = r;
	if(largest != i){
		swap(arr[i], arr[largest]);
		heapify(arr, largest, a, b);
	}
}

void makeHeap(vector<int> & arr, int a, int b){
	for(int i = parent(b, a); i >= a; i--){
		heapify(arr, i, a, b);
	}
}

void heapSort(vector<int> & arr, int a, int b){
	makeHeap(arr, a, b);
	while(a < b){
		swap(arr[a], arr[b]);
		heapify(arr, a, a, --b);
	}
}

//Intro Sort
void introSortRec(vector<int> & arr, int a, int b, int maxDepth){
	int size = b - a;
	if(size < 16){
		insertionSort(arr, a, b);
	}else if(maxDepth == 0){
		heapSort(arr, a, b);
	}else{
		int m = partition(arr, a, b);
		introSortRec(arr, a, m - 1, maxDepth - 1);
		introSortRec(arr, m + 1, b, maxDepth - 1);
	}
}

void introSort(vector<int> & arr, int a, int b){
	int maxDepth = 2 * log(b - a);
	introSortRec(arr, a, b, maxDepth);
}

//Bucket Sort
int maxValue(vector<int> & arr, int a, int b){
	int maxN = 0;
	for(int i = a; i <= b; i++){
		if(arr[i] > maxN) maxN = arr[i];
	}
	return maxN;
}

void bucketSort(vector<int> & arr, int a, int b){
	int maxN = maxValue(arr, a, b);
	vector<int> bucket(maxN + 1);
	for(int i = a; i <= b; i++){
		bucket[arr[i]]++;
	}
	int p = 0;
	for(int i = 0; i <= maxN; i++){
		for(int j = 0; j < bucket[i]; j++){
			arr[a + p++] = i;
		}
	}
}

//Counting Sort (only non-negative values)
void countingSort(vector<int> & arr, int a, int b, int exp = 0){
	int maxN;
	if(exp > 0) maxN = 9;
	else maxN = maxValue(arr, a, b);
	vector<int> counting(maxN + 1);
	for(int i = a; i <= b; i++){
		int pos = arr[i];
		if(exp > 0) pos = (pos / exp) % 10;
		counting[pos]++;
	}
	for(int i = 1; i <= maxN; i++){
		counting[i] += counting[i - 1];
	}
	vector<int> ans(b - a + 1);
	for(int i = b; i >= a; i--){
		int pos = arr[i];
		if(exp > 0) pos = (pos / exp) % 10;
		ans[counting[pos] - 1] = arr[i];
		counting[pos]--;
	}
	for(int i = a; i <= b; i++){
		arr[i] = ans[i - a];
	}
}

//Radix Sort (only non-negative values)
void radixSort(vector<int> & arr, int a, int b){
	int maxN = maxValue(arr, a, b);
	int exp = 1;
	while(maxN > 0){
		countingSort(arr, a, b, exp);
		exp *= 10;
		maxN /= 10;
	}
}

int random(int minimo, int maximo){
	return rand() % (maximo - minimo + 1) + minimo;
}

int main(){srand(time(0));
	int N = 2e7;
	clock_t begin, end;
	vector<int> arr(N), test(N);
	for(int i = 0; i < arr.size(); i++) arr[i] = random(0, (int)arr.size());

	test = arr;
	begin = clock();
	mergeSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "MergeSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	quickSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "QuickSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	heapSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "HeapSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	introSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "IntroSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	sort(test.begin(), test.end());
	end = clock();
	cout << fixed << setprecision(4) << "STL Sort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	bucketSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "BucketSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	countingSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "CountingSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";

	test = arr;
	begin = clock();
	radixSort(test, 0, test.size() - 1);
	end = clock();
	cout << fixed << setprecision(4) << "RadixSort:\t" << double(end - begin) / CLOCKS_PER_SEC << "s\n";
	return 0;
}