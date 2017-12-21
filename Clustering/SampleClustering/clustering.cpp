#include "clustering.h"
#include "makesample.h"

void clustering(const vector<Record >& vec, int nClus, double tau,  int* clus) {
	int Num = vec.size();
	double** matrix = new double*[Num];
	for (int i = 0; i < Num; ++i)
		*(matrix + i) = new double[Num - i]();
	//calculate the two-dimensional distance matrix
	
	distanceMatrix(vec, matrix, eucDistance);

	double* rho = new double[Num];

	double radius = searchRadius(matrix, Num, tau);
	density(matrix, Num, radius, rho);

	//get delta
	double* delta = new double[Num];
	int* neighbor = new int[Num];
	int* order = new int[Num];
	getDelta(matrix, Num, rho, delta, neighbor, order);

	vector<int> clustersVec;
	findInitialCenters(rho, delta, Num, nClus, clustersVec);
	//save the index of samples treated as cluster centers into file
	//assign cluster centers to samples
	assignClusters(vec,order, neighbor, Num, clustersVec, clus);
	for (int ii = 0; ii < Num; ++ii)	delete[] * (matrix + ii);
	delete[] matrix;
	delete[] rho;
	delete[] delta;
	delete[] neighbor;
	delete[] order;

}

void density(double** matrix, int Num, double radius,  double* rho) {
	double dist, val;
	memset(rho, 0, sizeof(double)*Num);//reset values in res
	vector<double> vec(Num, 0);
	
		for (int i = 0; i < Num; ++i)
			for (int j = i + 1; j < Num; ++j) {
				dist = getMatrixData(matrix, i, j);
				if (dist < radius) {
					*(rho + i) += 1;
					*(rho + j) += 1;
				}
			}
	
}

double getMatrixData(double** matrix, int i, int j) {
	int row = i < j ? i : j;
	int col = i > j ? i : j;
	col -= row;
	return *(*(matrix + row) + col);
}

void setMatrixData(double** matrix, int i, int j, double val) {
	int row = i < j ? i : j;
	int col = i > j ? i : j;
	col -= row;
	*(*(matrix + row) + col) = val;
}

void distanceMatrix(const vector<Record>& vec, double** matrix,
	double(*metricfun)(const Record&, const Record&)) {
	
	size_t sz = vec.size();
	double dist = 0.0;

	for (size_t i = 0; i < sz; ++i) {
		setMatrixData(matrix, i, i, 0);
		for (size_t j = i + 1; j < sz; ++j) {
			dist = metricfun(vec[i], vec[j]);
			setMatrixData(matrix, i, j, dist);
		}

	}

	//matdis.close();


}

double searchRadius(double** matrix, int Num, double tau) {
	int nElem = Num*(Num - 1) / 2;
	vector<double> dist;
	dist.reserve(nElem);
	int cnt = 0;
	for (int i = 0; i < Num - 1; ++i)
		for (int j = i + 1; j < Num; ++j)
			dist.push_back(getMatrixData(matrix, i, j));

	int pos = int(round(tau*nElem));//position of d_c
	nth_element(dist.begin(), dist.begin() + pos - 1, dist.end());
	return dist[pos - 1];
}

void sortByDensity(const double* rho, int Num, int* order) {
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


void getDelta(double** matrix, int Num, const double* rho, double* delta,
	int* neighbor, int* order) {
	sortByDensity(rho, Num, order);
	double globalMax = getMatrixData(matrix, 0, 0);
	for (int i = 0; i < Num; ++i) {
		double min = getMatrixData(matrix, order[i], order[0]);
		double buf = 0;
		*(neighbor + order[i]) = order[0];
		for (int j = 0; j < i; ++j) {
			buf = getMatrixData(matrix, order[i], order[j]);
			if (buf > globalMax) globalMax = buf;
			if (buf < min) {
				min = buf;
				*(neighbor + order[i]) = order[j];
			}
		}
		*(delta + order[i]) = min;
	}
	*(delta + order[0]) = globalMax;
}

int findInitialCenters(const double* rho, const double* delta,
	int Num, int nClus, vector<int>& vec) {
	if (NULL == rho || NULL == delta)
		return -1;
	vec.clear();
	//scale delta and rho into the range of [0,1]
	double rho_min, rho_max;
	rho_min = rho_max = *rho;
	for (int i = 1; i < Num; ++i) {
		if (*(rho + i) > rho_max) rho_max = *(rho + i);
		else if (*(rho + i) < rho_min) rho_min = *(rho + i);
	}
	double rho_range = rho_max - rho_min;

	double delta_min, delta_max;
	delta_min = delta_max = *delta;
	for (int ii = 1; ii < Num; ++ii) {
		if (*(delta + ii) > delta_max) delta_max = *(delta + ii);
		else if (*(delta + ii) < delta_min) delta_min = *(delta + ii);
	}
	double delta_range = delta_max - delta_min;

	double *pgamma = new double[Num];
	for (int t = 0; t < Num; ++t)
		*(pgamma + t) = (*(rho + t) - rho_min)*(*(delta + t) - delta_min) / (rho_range*delta_range);

	for (int s = 0; s < nClus; ++s) {
		int max = s;
		for (int j = s + 1; j < Num; ++j)
			if (*(pgamma + max) < *(pgamma + j)) max = j;
		std::swap(*(pgamma + s), *(pgamma + max));
		vec.push_back(max);
	}
	delete[] pgamma;
	return nClus;
}

void assignClusters(const vector<Record> &data_vec, const int* order, const int* neighbor,
	int Num, const vector<int>& centre, int* clus) {
	for (int j = 0; j < centre.size(); ++j)
		*(clus + data_vec[centre[j]].index) = j+1;
	for (int t = 0; t < Num; ++t)
		if (*(clus + data_vec[*(order + t)].index) == 0)//waiting for assignment
			*(clus + data_vec[*(order + t)].index) = *(clus + data_vec[*(neighbor + *(order + t))].index);
}


