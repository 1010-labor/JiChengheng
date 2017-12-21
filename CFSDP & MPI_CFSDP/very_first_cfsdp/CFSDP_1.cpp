/*
#include "CFSDP_1.h"
#define DATAFILE "d:\\data\\data.dat"
#define RESFILE "d:\\test\\result.dat"

void findKnn(double* dm_line, int left, int right, int k) {

	int current_find = right - left + 1;
	double key;
	int lindex;
	int rindex = right;
	int key_history;


	while (current_find > k + 1) {

		lindex = left;//rindex实指最后一个元素
		key_history = rindex;
		//记录left作为哨兵
		key = dm_line[left];


		while (lindex < rindex) {
			//逆向寻找比哨兵小的元素
			while (lindex < rindex && key <= dm_line[rindex]) rindex--;
			//找到比哨兵小的元素 移动到左索引位置，左索引前进一位
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			//正向寻找比哨兵大的元素
			while (lindex < rindex&&key>= dm_line[lindex]) lindex++;
			//找到比哨兵大的元素 移动到右索引位置，右索引后退一位
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
DistFunction distfunc = EuclideanDistance;
void main() {
	double max_density, min_density;
	double max_deta, min_deta;


	ofstream timelist("d:\\time.txt");
	clock_t start, finish;
	clock_t start_0;
	double totaltime;
	start_0 = clock();
	start = clock();

	vector<vector<double>> data_vec;
	ReadData(DATAFILE, data_vec);
	int data_num = data_vec.size();

	finish = clock();
	cout << "read data:\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	timelist << "read data:\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();



	double** dm = new double*[data_vec.size()];
	for (int i = 0; i < data_num; i++) {
		*(dm + i) = new double[data_num - i];

	}
	start = clock();
	CreateDM(data_vec, dm, distfunc);

	finish = clock();
	cout << "creat dm :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	timelist << "creat dm :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	double dc = FindCutDist(dm, data_num);
	cout << endl << dc << endl;

	finish = clock();
	cout << "find dc :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	timelist << "find dc :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	double* density_vec=new double[data_num];
	int* order=new int[data_num];
	DensityVec(dm, data_num, density_vec, dc);
	finish = clock();
	cout << "compute density :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	timelist << "compute  density :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();
	SortDensity(density_vec, order,data_num);

	finish = clock();
	cout << "sort density :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	timelist << "sort density :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	vector<double> deta(data_num);
	vector<int> neighbor(data_num);
	DetaAndNeighbor(dm, order, deta, neighbor,min_deta);


	max_deta = deta[0];
	max_density = density_vec[order[0]];
	min_density= density_vec[order[data_num-1]];

	finish = clock();
	cout << "compute deta :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	timelist << "compute deta :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	

	vector<double> gama(data_num);
	GamaVec(density_vec, deta, min_deta, min_density, max_deta, max_density, gama);

	finish = clock();
	timelist << "gama :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	vector<int> result(data_num, 0);
	FindInitialCenters(gama, 4, result);
	AssignClusters(order, neighbor, result);

	finish = clock();
	timelist << "result :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	ofstream output(RESFILE);
	for (int i = 0; i < data_num; i++) {
		output << i << " " << result[i] << " " << 0 << endl;
	}
	output.close();

	finish = clock();
	timelist << "output :\t" << (double)(finish - start) / CLOCKS_PER_SEC << endl;
	start = clock();

	for (int i = 0; i < data_num; i++) {
		delete[] * (dm + i);
	}
	finish = clock();
	cout << "total :\t" << (double)(finish - start_0) / CLOCKS_PER_SEC << endl;
	start = clock();

	timelist.close();







}
void SetMatrixData(double** dist_matrix, int i, int j, double val) {
	int row = i < j ? i : j;
	int col = i > j ? i : j;
	col -= row;
	*(*(dist_matrix + row) + col) = val;
}
double GetMatrixData(double** dist_matrix, int i, int j) {
	int row = i < j ? i : j;
	int col = i > j ? i : j;
	col -= row;
	return *(*(dist_matrix + row) + col);
}
double EuclideanDistance(const vector<double>& vec1, const vector<double>& vec2) {

	double res = 0.0;
	vector<double>::const_iterator iter1 = vec1.begin();
	vector<double>::const_iterator iter2 = vec2.begin();
	while (iter1 != vec1.end() && iter2 != vec2.end()) {
		res += (*iter1 - *iter2)*(*iter1 - *iter2);
		++iter1;
		++iter2;
	}
	return sqrt(res);
}
void ReadData(string file_name, vector<vector<double>> &data_set) {
	stringstream ss;
	ifstream file(file_name);
	if (!file) {
		cout << "file open failed!";
		return;
	}
	string line;
	double value;
	while (getline(file, line)) {
		ss.clear();
		ss << line;
		vector<double> record;
		while (ss >> value)	record.push_back(value);
		data_set.push_back(record);
	}
	file.close();
}
void CreateDM(const vector<vector<double>>&  data_vec, double** dist_matrix, DistFunction  dist_func) {
	int s = data_vec.size();
	for (int i = 0; i < s; i++) {
		SetMatrixData(dist_matrix, i, i, 0);
		for (int j = i + 1; j < s; j++) {
			SetMatrixData(dist_matrix, i, j, dist_func(data_vec[i], data_vec[j]));
		}

	}


}
double FindCutDist(double** dist_matrix, int data_num) {
	int nElem = data_num*(data_num - 1) / 2;
	double* dist=new double[nElem];

	for (int i = 0, k = 0; i < data_num - 1; ++i)
		for (int j = i + 1; j < data_num; ++j) {
			dist[k++] = GetMatrixData(dist_matrix, i, j);
		}
	int pos = int(round(0.0175*nElem));//position of d_c
	findKnn(dist, 0, nElem, pos);
	return dist[pos - 1];

}
void DensityVec(double** dist_matrix, int data_num, double* density_vec, double cut_dist) {
	for (int i = 0; i < data_num; i++) {
		density_vec[i] = -1;
		for (int j = 0; j < data_num; j++) {
			if (GetMatrixData(dist_matrix, i, j) <= cut_dist) {
				density_vec[i]++;
			}
		}
	}

}
void SortDensity(double* density_vec, int* order,int data_num) {
	for (int t = 0; t < data_num; ++t)
		order[t] = t;
	for (int i = 0; i < data_num - 1; ++i) {
		int max = i;
		for (int j = i + 1; j < data_num; ++j)
			if (density_vec[ order[max]] < density_vec [ order[j]])
				max = j;
		std::swap(order[i], order[max]);
	}

}
void DetaAndNeighbor(double** dist_matrix, int* order, vector<double> &deta, vector<int> &neighbor, double &min_deta) {
	int Num = deta.size();
	double globalMax = GetMatrixData(dist_matrix, 0, 0);
	min_deta = GetMatrixData(dist_matrix, order[1], order[0]);
	double min;
	double buf;
	for (int i = 1; i < Num; ++i) {
		min = GetMatrixData(dist_matrix, order[i], order[0]);
		buf = 0;
		neighbor[order[i]] = order[0];
		for (int j = 0; j < i; ++j) {
			buf = GetMatrixData(dist_matrix, order[i], order[j]);
			if (buf > globalMax) globalMax = buf;
			if (buf < min) {
				min = buf;
				neighbor[order[i]] = order[j];
			}
		}
		deta[order[i]] = min;
		if (min < min_deta) min_deta = min;

	}
	deta[order[0]] = globalMax;
}
void GamaVec(double* density_vec, const  vector<double> &deta, double min_deta, double min_density, double max_deta, double max_density, vector<double> &gama) {
	int Num = deta.size();
	double divisor = (max_density - max_deta)*(max_deta - min_deta);
	for (int i = 0; i < Num; i++)
		gama[i] = (density_vec[i] - min_density)*(deta[i] - min_deta) / divisor;



}
void FindInitialCenters(vector<double> &gama, int n_clus, vector<int>& res_vec) {
	int first_id = 1;
	for (int s = 0; s < n_clus; ++s) {
		int max = 0;
		for (int j = 0; j < gama.size(); ++j) {
			if (gama[max] < gama[j])
				max = j;
		}
		gama[max] = 0;


		res_vec[max] = first_id++;
	}

}
void AssignClusters(int* order, const vector<int> &neiborgh, vector<int>& res_vec) {
	int Num = neiborgh.size();
	for (int i = 0; i < Num; i++) {
		if (res_vec[order[i]] == 0)

			res_vec[order[i]] = res_vec[neiborgh[order[i]]];
	}
}

*/
