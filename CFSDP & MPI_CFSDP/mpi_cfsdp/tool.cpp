//#include "tool.h"
//
//
//
//template<typename T>			//此种实现方式要比下面的交换方式  快n多
//void find_KNN(vector<T> &vec, int left, int right, int n) {
//	T key;
//	int lindex;
//	int rindex;
//	int lens = right;
//
//	while (right - left > n) {
//		key = vec[left];
//		lindex = left;
//		rindex = right - 1;
//		while (lindex < rindex) {
//			while (lindex < rindex && key<vec[rindex]) rindex--;
//			if (lindex < rindex)  vec[lindex++] = vec[rindex];//find  不是正常退出
//			while (lindex < rindex&&key>vec[lindex]) lindex++;
//			if (lindex < rindex) 	vec[rindex--] = vec[lindex];
//		}
//		vec[lindex] = key;
//		right = lindex;
//	}
//	if (right - left == n || right - left == n - 1) return;
//	else {
//		find_KNN(vec, right + 1, lens, n - (right - left + 1));
//	}
//}
///*
//
//template<typename RandomAccessIterator, typename Comp>
//void find_KNN(RandomAccessIterator left, RandomAccessIterator right, int n, Comp comp) {
//RandomAccessIterator pivot;
//RandomAccessIterator lindex;
//RandomAccessIterator rindex;
//int lens = right;
//
//while (right - left > n) {
//pivot = left;
//lindex = left;
//rindex = right - 1;
//while (lindex < rindex) {
//while (lindex < rindex && comp(*pivot, *rindex)) rindex--;
//if (lindex < rindex)  iter_swap(lindex++, rindex);//find  不是正常退出
//while (lindex < rindex&&!comp(*pivot, *lindex))lindex++;
//if (lindex < rindex) 	iter_swap(rindex--, lindex);
//}
//right = lindex;
//}
//if (right - left == n || right - left == n - 1) return;
//else {
//find_KNN(right + 1, left + lens, n - (right - left + 1), comp);
//}
//
//}*/
//template<typename RandomAccessIterator>
//void find_KNN(RandomAccessIterator left, RandomAccessIterator right, int n) {
//	RandomAccessIterator pivot;
//	RandomAccessIterator lindex;
//	RandomAccessIterator rindex;
//	int lens = right - left;
//
//	while (right - left > n) {
//		pivot = left;
//		lindex = left;
//		rindex = right - 1;
//		while (lindex < rindex) {
//			while (lindex < rindex && *pivot < *rindex) rindex--;
//			if (lindex < rindex)  iter_swap(lindex++, rindex);//find  不是正常退出 
//			pivot = rindex;
//			while (lindex < rindex&&*pivot>*lindex)lindex++;
//			if (lindex < rindex) 	iter_swap(rindex--, lindex);
//			pivot = lindex;
//		}
//		right = lindex;
//	}
//	if (right - left == n || right - left == n - 1) return;
//	else {
//		find_KNN(right + 1, left + lens, n - (right - left + 1));
//	}
//
//
//}
//
///*
//void main() {
//
//	clock_t s1, s2, f1, f2, t1, t2;
//	s1 = s2 = f1 = f2 = t1 = t2 = 0;
//
//	vector<int>	a, b;
//	int size = 2000;
//	a.reserve(size);
//	b.reserve(size);
//	for (int i = 0; i < size; i++) {
//		a.clear();
//		b.clear();
//		srand(clock());
//		int s, db;
//		for (int i = 0; i < size; i++) {
//			s = rand() % 50000;
//			db = s;
//			a.push_back(s);
//			b.push_back(db);
//
//		}
//		a[0] *= 10;
//		b[0] *= 10;
//		
//		s1 = clock();
//		find_KNN(a, 0, a.size(), 20);
//		// sort(a.begin(), a.end(), less<int>());
//		f1 = clock();
//		t1 += f1 - s1;
//		
//		s2 = clock();
//		sort(b.begin(), b.end());
//		f2 = clock();
//		t2 += f2 - s2;
//
//	}
//	for (int i = 0; i < 20; i++) {
//		cout << a[i] << "  ";
//	}
//	cout << endl;
//	cout << endl;
//	for (int i = 0; i < 20; i++) {
//		cout << b[i] << "  ";
//	}
//
//
//		cout << endl;
//
//		cout << "sort============" << t1 << "'s\n";
//		cout << "KNN============" << t2 << "'s\n";
//
//	}*/
//	
//
//
//
//
//
//
//
//
//
//
