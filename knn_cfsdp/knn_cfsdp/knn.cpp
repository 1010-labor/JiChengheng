#include "func.h"

bool distCmp(Dist a, Dist b) {
	return a.distance < b.distance;
}

void sortDist(Dist* dm_line, int left, int right) {
	if (right > left) {
		Dist key = dm_line[left];
		int lindex = left;
		int rindex = right;
		while (lindex < rindex) {
			while (lindex < rindex && key.distance <= dm_line[rindex].distance) rindex--;
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			while (lindex < rindex&&key.distance >= dm_line[lindex].distance) lindex++;
			if (lindex < rindex) 	dm_line[rindex--] = dm_line[lindex];
		}
		if (lindex != rindex)
			cout << "";
		dm_line[lindex] = key;
		sortDist(dm_line, left, lindex - 1);
		sortDist(dm_line, rindex + 1, right);
	} else {
		return;
	}

}


void findKnn(Dist* dm_line, int left, int right, int k) {

	int current_find = right - left + 1;
	Dist key;
	int lindex;
	int rindex = right;
	int key_history;

	while (current_find > k + 1) {

		lindex = left;
		key_history = rindex;
		key = dm_line[left];
		while (lindex < rindex) {
			while (lindex < rindex && key.distance <= dm_line[rindex].distance) rindex--;
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			while (lindex < rindex&&key.distance >= dm_line[lindex].distance) lindex++;
			if (lindex < rindex) 	dm_line[rindex--] = dm_line[lindex];
		}
		assert(lindex == rindex);
		dm_line[lindex] = key;
		current_find = rindex - left + 1;
	}
	if (current_find == k + 1 || current_find == k) return;
	else {
		findKnn(dm_line, rindex + 1, key_history, k - current_find);
	}

}

void readData(double* x, double* y, int num, string filePath) {
	fstream file(filePath);
	char delimiter;
	for (int i = 0; i < num; i++) {
		file >> x[i];
		file >> y[i] >> delimiter;
	}
}

void getRhoAndSm(double* x, double* y, int num, int k, int shortest, double* rho, Dist** sm) {
	Dist* dline = new Dist[num];
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < num; j++) {
			dline[j] = Dist(i, pow((pow((x[i] - x[j]), 2), pow((y[i] - y[j]), 2)), 0.5));
		}
		findKnn(dline, 0, num - 1, k + 1);
		double tmp_sum = 0;
		for (int j = 0; j < k; j++) {
			tmp_sum += dline[j].distance;
		}
		rho[i] = 1.0 / tmp_sum;
		sort(dline, dline + k, distCmp);
		memcpy(sm[i], dline + 1, shortest*sizeof(Dist));
	}

	delete[] dline;

}

void getDetaAndNeighbor(double* x, double* y, int num, int shortest, Dist** sm, int* seq_of_rho, double* rho, double* deta, int* neighbor) {
	int first = seq_of_rho[0];
	for (int i = 0; i < num; i++) {
		int index_this = seq_of_rho[i];
		double rho_this = rho[index_this];
		int neighbor_this = first;
		double deta_this = pow((pow((x[i] - x[first]), 2), pow((y[i] - y[first]), 2)), 0.5);
		bool findInKnn = false;

		int tmp_with;
		double tmp_dist;
		for (int j = 0; j < shortest; j++) {
			tmp_with = ((*sm + i) + j)->with;
			tmp_dist = ((*sm + i) + j)->distance;
			if (tmp_dist < deta_this && rho_this < rho[tmp_with]);
			{
				deta[index_this] = tmp_dist;
				neighbor[index_this] = tmp_with;
				findInKnn = true;
				break;
			}
		}
		if (!findInKnn) {
			for (int j = 1; j < i; j++) {
				tmp_with = seq_of_rho[j];
				tmp_dist = pow((pow((x[i] - x[j]), 2), pow((y[i] - y[j]), 2)), 0.5);
				if (tmp_dist < deta_this && rho_this < rho[tmp_with]) {
					deta_this = tmp_dist;
					neighbor_this = tmp_with;
				}
			}
			deta[index_this] = deta_this;
			neighbor[index_this] = neighbor_this;
		}
	}
}


