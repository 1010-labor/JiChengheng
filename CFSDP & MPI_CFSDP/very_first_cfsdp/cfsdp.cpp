/*
#include "cfsdp.h"
using namespace std;
#define CLUSNUM 4
#define DATAFILE "d:\\data\\data.dat"

void assignClusters(const int* order, const int* neighbor,
	int Num, const vector<int>& vec, int* res) {

	memset(res, 0, Num*sizeof(int));
	for (int sz = 1; sz <= vec.size(); ++sz)
		*(res + vec[sz - 1]) = sz;
	/ *cout << neighbor[order[1]] <<  "  mi du di er  da de lin ju "<<order[0]<<endl;* /
	for (int t = 0; t < Num; ++t) {
		if (res[order[t]] == 0)
			res[order[t]] = res[neighbor[order[t]]];
	}

}

void findKnn(dist* dm_line, int left, int right, int k) {

	int current_find = right - left + 1;
	dist key;
	int lindex;
	int rindex = right;
	int key_history;


	while (current_find > k + 1) {

		lindex = left;//rindexʵָ���һ��Ԫ��
		key_history = rindex;
		//��¼left��Ϊ�ڱ�
		key = dm_line[left];


		while (lindex < rindex) {
			//����Ѱ�ұ��ڱ�С��Ԫ��
			while (lindex < rindex && key.distance <= dm_line[rindex].distance) rindex--;
			//�ҵ����ڱ�С��Ԫ�� �ƶ���������λ�ã�������ǰ��һλ
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			//����Ѱ�ұ��ڱ����Ԫ��
			while (lindex < rindex&&key.distance >= dm_line[lindex].distance) lindex++;
			//�ҵ����ڱ����Ԫ�� �ƶ���������λ�ã�����������һλ
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


//calculate the distance between two vectors with Euclidean metric
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

void readData(string dataFile, vector<vector<double>>& dataVec) {
	stringstream ss;
	ifstream ifs(dataFile);
	double val;
	if (ifs) {
		dataVec.clear();
		string line;
		vector<double> subvec;

		while (getline(ifs, line)) {
			subvec.clear();
			ss.clear();
			ss << line;
			while (ss >> val) {
				subvec.push_back(val);

			}
			dataVec.push_back(subvec);



		}
	} else
		cerr << dataFile << " doesn't exist!\n";
}

void sortedOrder(const double* rho, int Num, int* order) {
	for (int t = 0; t < Num; ++t)
		order[t] = t;
	for (int i = 0; i < Num - 1; ++i) {
		int max = i;
		for (int j = i + 1; j < Num; ++j)
			if (*(rho + order[max]) < *(rho + order[j]))
				max = j;
		std::swap(order[i], order[max]);
	}
}


void enlargeDataSet(vector<vector<double>>& dataVec, int times, double** newset) {
	int d = 0;
	for (vector<vector<double>>::iterator it = dataVec.begin(); it != dataVec.end(); it++, d++) {
		for (int i = 0; i < times; i++) {

			double val = 0.0;
			int index = 0;
			double p;
			for (vector<double > ::iterator it2 = (*it).begin(); it2 != (*it).end(); it2++, index++) {
				p = 1.0 + rand() % 100 / (double)100 * (rand() % 55 / (double)100);
				val = (*it2)*p;
				newset[d*times + i][index] = val;
			}

		}
	}

}


int main(int argc, char* argv[]) {

	clock_t start, finish;
	double totaltime;
	start = clock();


	//Read dataset
	vector<vector<double>> dataVec;
	readData(DATAFILE, dataVec);
	int dataNum = dataVec.size();


	//�����̿�ʼ�������Ӧ���������ľ������  ����һ�������ܶ�

	dist** partDistMatrix = new dist*[dataNum];
	for (int j = 0; j < dataNum; j++) {
		partDistMatrix[j] = new dist[dataNum];
	}
	for (int i = 0; i < dataNum; i++) {
		for (int j = 0; j < dataNum; j++) {
			(*(partDistMatrix + i) + j)->with = j;
			(*(partDistMatrix + i) + j)->distance = EuclideanDistance(dataVec[i], dataVec[j]);
		}
	}




	//knn��ʽ�����ܶ�
	double* density = new double[dataNum];
	int K = (int)((double)dataNum*0.015);              //�����ܶ�ʱ��ȡ���ھ���


	for (int i = 0; i < dataNum; i++) {
		findKnn(*(partDistMatrix + i), 0, dataNum - 1, K);
		double tmp = 0;
		for (int j = 0; j < K; j++) {
			tmp += (*(partDistMatrix + i) + j)->distance;
		}
		density[i] = dataNum / tmp;

	}





	//�����ܶ�����
	int* order = new int[dataNum];
	sortedOrder(density, dataNum, order);






	//�����ܶ�����  �����̼�����������������detaֵ���Լ��ܶȴ�������������ھ�   ��ͳ��deta�ľֲ��ֲ�����
	double* deta = new double[dataNum];
	int* neighbor = new int[dataNum];
	double deta_max = DBL_MIN;
	double deta_min = DBL_MAX;
	double deta_range = 0.0;

	deta[order[0]] = 0;
	for (int i = 0; i < dataNum; i++) {
		double  min_dist = EuclideanDistance(dataVec[i], dataVec[order[0]]);;
		int  min_with = order[0];
		int tmp_with;
		double tmp_dist;
		double density_of_i = density[i];
		bool   find_in_knn = false;
		for (int j = 0; j < K; j++) {
			tmp_with = (*(partDistMatrix + i) + j)->with;
			tmp_dist = (*(partDistMatrix + i) + j)->distance;
			if (tmp_dist<min_dist && density[tmp_with]>density_of_i) {
				min_dist = tmp_dist;
				min_with = tmp_with;
				find_in_knn = true;
			}

			if (!find_in_knn) {
				for (int k = 0; order[k] != i; k++) {
					tmp_dist = EuclideanDistance(dataVec[i], dataVec[order[k]]);
					if (tmp_dist < min_dist) {
						min_dist = tmp_dist;
						min_with = order[k];
					}
				}

			}

			deta[i] = min_dist;
			neighbor[i] = min_with;

			deta_max = deta_max < min_dist ? min_dist : deta_max;
			deta_min = deta_min > min_dist ? min_dist : deta_min;
		}

	} 
	deta_range = deta_max - deta_min;


	//�����̼������������������gamaֵ
	double density_max = density[order[0]];
	double density_min = density[order[dataNum - 1]];
	double density_range = density_max - density_min;

	double* gama = new double[dataNum];
	memset(gama, 0, sizeof(double)*dataNum);
	
		deta[order[0]] = deta_max;
		for (int i = 0; i < dataNum; i++) {
			gama[i] = (density[i] - density_min) * (deta[i] - deta_min) / density_range*deta_range;
		}
	



	//������ͬ��gama  ѡȡ��������   ����
	
		int* order2 = new int[dataNum];
		for (int j = 0; j < dataNum; j++) {
			order2[j] = j;
		}
		sortedOrder(gama, dataNum, order2);
		int* clus = new int[dataNum];
		memset(clus, 0, sizeof(int)*dataNum);
		vector<int> centres;
		for (int i = 0; i < CLUSNUM; i++)	centres.push_back(order2[i]);

		assignClusters(order, neighbor, dataNum, centres, clus);
		ofstream out("d:\\result");
		for (int i = 0; i < dataNum; i++) {

			out << i << " " << *(clus + i) << " 0" << endl;
		}
		out.close();
	
	finish = clock();
	cout << (double)(finish - start) / CLOCKS_PER_SEC << " seconds costs by the procedure !" << endl;

}




*/