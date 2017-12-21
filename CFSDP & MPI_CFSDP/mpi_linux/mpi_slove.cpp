#include "mpi_solve.h"
#include "iomanip"
#pragma comment (lib, "msmpi.lib")  
using namespace std;
#define SHORTEST 10
#define CLUSNUM 4
#define DATAFILE "f:\\DATA\\data.dat"
#define RES "e:\\RESULT_MPI.DAT"

void assignClusters(const int* order, const int* neighbor,
	int Num, const vector<int>& vec, int* res) {

	memset(res, 0, Num*sizeof(int));
	for (int sz = 1; sz <= vec.size(); ++sz)
		*(res + vec[sz - 1]) = sz;
	/*cout << neighbor[order[1]] <<  "  mi du di er  da de lin ju "<<order[0]<<endl;*/
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

		lindex = left;//rindex实指最后一个元素
		key_history = rindex;
		//记录left作为哨兵
		key = dm_line[left];


		while (lindex < rindex) {
			//逆向寻找比哨兵小的元素
			while (lindex < rindex && key.distance <= dm_line[rindex].distance) rindex--;
			//找到比哨兵小的元素 移动到左索引位置，左索引前进一位
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			//正向寻找比哨兵大的元素
			while (lindex < rindex&&key.distance >= dm_line[lindex].distance) lindex++;
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
double EuclideanDistance2(double* vec1, double* vec2, int dim) {

	double res = 0.0;
	while (dim--) {
		res += (vec2[dim] - vec1[dim])*(vec2[dim] - vec1[dim]);

	}
	return sqrt(res);
}
void readData(string dataFile, vector<vector<double> >& dataVec) {
	stringstream ss;
	ifstream ifs(dataFile.c_str());
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


void enlargeDataSet(vector<vector<double> >& dataVec, int times, double** newset) {
	int d = 0;
	for (vector<vector<double> >::iterator it = dataVec.begin(); it != dataVec.end(); it++, d++) {
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

	clock_t start, finish, start0;
	double totaltime;
	start = clock();
	start0 = clock();
	int id, size;
	ofstream timeres("E:\\time_res.dat");

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	//Read dataset
	FILE * fp = NULL; //文件指针。
	int c, lc = 0; //c为文件当前字符，lc为上一个字符，供结尾判断用。
	int line = 0; //行数统计
	fp = fopen(DATAFILE, "r");//以只读方式打开文件。
	while ((c = fgetc(fp)) != EOF) //逐个读入字符直到文件结尾
	{
		if (c == '\n') line++; //统计行数。
		lc = c; //保存上一字符。
	}
	if (lc != '\n') line++;//处理末行
	fseek(fp, 0L, SEEK_SET);

	int dim = 2;
	int dataNum = line;
	double** dataVec = new double*[dataNum];
	for (int i = 0; i < dataNum; i++) {
		dataVec[i] = new double[dim];

		(fscanf(fp, "%lf", (*(dataVec + i))));
		fgetc(fp);//过滤掉分隔符
		(fscanf(fp, "%lf", (*(dataVec + i) + 1)));
		fgetc(fp);//过滤掉分隔符


	}
	
	int partSize = ceil(dataNum / size);


	//各进程开始计算其对应部分样本的距离矩阵  并进一步计算密度

	dist* dist_line = new dist[dataNum];
	dist** shortest_matrix = new dist*[partSize];
	for (int j = 0; j < partSize; j++) {
		shortest_matrix[j] = new dist[SHORTEST];
	}

	

	double* density = new double[dataNum];
	int K = (int)((double)dataNum*0.015);



	if (id != size - 1) {
		for (int i = id*partSize; i < (id + 1)*partSize; i++) {
			for (int j = 0; j < dataNum; j++) {
				(dist_line + j)->with = j;
				(dist_line + j)->distance = EuclideanDistance2(dataVec[i], dataVec[j],2);
			}
			findKnn(dist_line, 0, dataNum - 1, K);
			findKnn(dist_line, 0, K, SHORTEST);
			double tmp = 0;
			for (int j = 0; j < K; j++) {
				tmp += (dist_line + j)->distance;
				
			}
			for (int j = 0; j < SHORTEST; j++) {
				shortest_matrix[i%partSize][j] = dist_line[j];
			}
			
			density[i] = dataNum / tmp;

		}
	} else {
		for (int i = id*partSize; i < (id + 1)*partSize; i++) {
			for (int j = 0; j < dataNum; j++) {
				(dist_line + j)->with = j;
				(dist_line + j)->distance = EuclideanDistance2(dataVec[i], dataVec[j], 2);
			}
			findKnn(dist_line, 0, dataNum - 1, K);
			double tmp = 0;
			for (int j = 0; j < K; j++) {
				tmp += (dist_line + j)->distance;

			}
			for (int j = 0; j < SHORTEST; j++) {
				shortest_matrix[i%partSize][j] = dist_line[j];
			}

			density[i] = dataNum / tmp;

		}
	}

	
	if (id == size - 1) {
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "DENSITY & SHORTEST MATRIX    -" << endl;
		start = clock();
		cout<<"ok here"<<endl;
	}
	

	//同步密度  计算密度排序  
	int* order = new int[dataNum];

	if (id < size - 1) {
		MPI_Send((density + id*partSize), partSize, MPI_DOUBLE, size - 1, id*partSize, MPI_COMM_WORLD);
	}
	if (id == size - 1) {
		MPI_Status status;
		for (int i = 0; i < size - 1; i++) {
			MPI_Recv(density + i*partSize, partSize, MPI_DOUBLE, i, i*partSize, MPI_COMM_WORLD, &status);
		}
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "  -DENSITY synch    " << endl;
		start = clock();
		for (int i = 0; i < dataNum; i++) {
			order[i] = i;
		}
		sortedOrder(density, dataNum, order);


	}

	//同步密度和密度排序

	MPI_Bcast(density, dataNum, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);
	MPI_Bcast(order, dataNum, MPI_INT, size - 1, MPI_COMM_WORLD);

	if (id == size - 1) {
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "  -DENSITY & ORDER bcast    " << endl;
		start = clock();
	}


	//根据密度排序  各进程计算其所负责样本的deta值、以及密度大于自身的最近邻居   并统计deta的局部分布区间
	double* deta = new double[dataNum];
	int* neighbor = new int[dataNum];
	double deta_max = 0;
	double deta_min = 100000;
	double other_deta_max = 0;
	double other_deta_min = 10000;
	double deta_range = 0.0;
	deta[order[0]] = 0;
	if (id != size - 1) {
		for (int i = id*partSize; i < (id + 1)*partSize; i++) {
			double  min_dist = EuclideanDistance2(dataVec[i], dataVec[order[0]],2);;
			int  min_with = order[0];
			int tmp_with;
			double tmp_dist;
			double density_of_i = density[i];
			bool   find_in_knn = false;
			for (int j = 0; j < SHORTEST; j++) {
				tmp_with = (*(shortest_matrix + i%partSize) + j)->with;
				tmp_dist = (*(shortest_matrix + i % partSize) + j)->distance;
				if (tmp_dist<min_dist && density[tmp_with]>density_of_i) {
					min_dist = tmp_dist;
					min_with = tmp_with;
					find_in_knn = true;
				}
			}
			if (!find_in_knn) {
				for (int k = 0; order[k] != i; k++) {
					tmp_dist = EuclideanDistance2(dataVec[i], dataVec[order[k]],2);
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

	} else {
		for (int i = id*partSize; i < dataNum; i++) {
			double  min_dist = EuclideanDistance2(dataVec[i], dataVec[order[0]],2);;
			int  min_with = order[0];
			int tmp_with;
			double tmp_dist;
			double density_of_i = density[i];
			bool   find_in_knn = false;
			for (int j = 0; j <SHORTEST; j++) {
				tmp_with = (*(shortest_matrix + i % partSize) + j)->with;
				tmp_dist = (*(shortest_matrix + i % partSize) + j)->distance;
				if (tmp_dist<min_dist && density[tmp_with]>density_of_i) {
					min_dist = tmp_dist;
					min_with = tmp_with;
					find_in_knn = true;
				}
			}
			if (!find_in_knn) {
				for (int k = 0; order[k] != i; k++) {
					tmp_dist = EuclideanDistance2(dataVec[i], dataVec[order[k]],2);
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

	MPI_Reduce(&deta_max, &other_deta_max, 1, MPI_DOUBLE, MPI_MAX, size - 1, MPI_COMM_WORLD);
	MPI_Reduce(&deta_min, &other_deta_min, 1, MPI_DOUBLE, MPI_MIN, size - 1, MPI_COMM_WORLD);

	//统计deta的全局分布区间，全局最小deta值  并广播     
	if (id == size - 1) {
		deta_max = deta_max < other_deta_max ? other_deta_max : deta_max;
		deta_min = deta_min > other_deta_min ? other_deta_min : deta_min;

		//修正密度最大样本的deta值
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "  -DETA    " << endl;
		start = clock();

	}
	MPI_Bcast(&deta_max, 1, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);
	MPI_Bcast(&deta_min, 1, MPI_DOUBLE, size - 1, MPI_COMM_WORLD);
	deta[order[0]] = deta_max;
	deta_range = deta_max - deta_min;
	if (id == size - 1) {
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "  -DETA synch    " << endl;
		start = clock();
	}

	//各进程计算其所负责样本点的gama值  
	double density_max = density[order[0]];
	double density_min = density[order[dataNum - 1]];
	double density_range = density_max - density_min;

	double* gama = new double[dataNum];
	memset(gama, 0, sizeof(double)*dataNum);

	if (id != size - 1) {

		for (int i = id*partSize; i < (id + 1)*partSize; i++) {
			gama[i] = (density[i] - density_min) * (deta[i] - deta_min) / density_range*deta_range;
		}
	} else {
		for (int i = id*partSize; i < dataNum; i++) {
			gama[i] = (density[i] - density_min) * (deta[i] - deta_min) / density_range*deta_range;
		}
	}

	if (id == size - 1) {
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "  GAMA    " << endl;
		start = clock();
	}
	//主进程同步gama  选取聚类中心   归类
	if (id < size - 1) {
		MPI_Send((gama + id*partSize), partSize, MPI_DOUBLE, size - 1, id*partSize, MPI_COMM_WORLD);
		MPI_Send((neighbor + id*partSize), partSize, MPI_INT, size - 1, id*partSize, MPI_COMM_WORLD);
	}
	if (id == size - 1) {
		MPI_Status status;
		for (int i = 0; i < size - 1; i++) {
			MPI_Recv(gama + i*partSize, partSize, MPI_DOUBLE, i, i*partSize, MPI_COMM_WORLD, &status);
			MPI_Recv(neighbor + i*partSize, partSize, MPI_INT, i, i*partSize, MPI_COMM_WORLD, &status);
		}
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << endl;
		start = clock();
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
		ofstream out(RES);
		for (int i = 0; i < dataNum; i++) {

			out << i << " " << *(clus + i) << " 0" << endl;
		}
		out.close();
		finish = clock();
		timeres << setw(20) << (double)(finish - start) / CLOCKS_PER_SEC << "  CLUSTERING    " << endl;
		start = clock();
		finish = clock();
		timeres << setw(20) << (double)(finish - start0) / CLOCKS_PER_SEC << "  CLUSTERING    " << endl;
		start = clock();
	}



	MPI_Finalize();
	finish = clock();

}




