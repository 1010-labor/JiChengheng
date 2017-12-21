#pragma once
#include "st.h"



//由对称矩阵获取数据
double getMatrixData(double** matrix, int i, int j);

//设定矩阵指定位置处的值
void setMatrixData(double** matrix, int i, int j, double val);

//计算二维距离矩阵
void distanceMatrix(const vector<Record>& vec, double** matrix,
	double(*metricfun)(const Record&, const Record&));

//寻找合适搜索半径
double searchRadius(double** matrix, int Num, double tau);

//计算各个样本密度
void density(double** matrix, int Num, double radius,  double* rho);

//计算单个样本的deta值
void getDelta(double** matrix, int Num, const double* dens, double* delta,
	int* neighbor, int* order);

void sortByDensity(const double* dens, int Num, int* index);

//find initial nClus cluster centers
int findInitialCenters(const double* dens, const double* delta, int Num,
	int nClus, vector<int>& vec);

//将各个样本聚集到相应的聚类中心
void assignClusters(const vector<Record> &data_vec,const int* order, const int* neighbor, int Num,
	const vector<int>& vec, int* res);

//algorithm of clustering 
//聚类算法
void clustering(const vector<Record >& vec, int nClus, double cut_tau, int* clus);



