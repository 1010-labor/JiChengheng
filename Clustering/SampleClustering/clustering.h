#pragma once
#include "st.h"



//�ɶԳƾ����ȡ����
double getMatrixData(double** matrix, int i, int j);

//�趨����ָ��λ�ô���ֵ
void setMatrixData(double** matrix, int i, int j, double val);

//�����ά�������
void distanceMatrix(const vector<Record>& vec, double** matrix,
	double(*metricfun)(const Record&, const Record&));

//Ѱ�Һ��������뾶
double searchRadius(double** matrix, int Num, double tau);

//������������ܶ�
void density(double** matrix, int Num, double radius,  double* rho);

//���㵥��������detaֵ
void getDelta(double** matrix, int Num, const double* dens, double* delta,
	int* neighbor, int* order);

void sortByDensity(const double* dens, int Num, int* index);

//find initial nClus cluster centers
int findInitialCenters(const double* dens, const double* delta, int Num,
	int nClus, vector<int>& vec);

//�����������ۼ�����Ӧ�ľ�������
void assignClusters(const vector<Record> &data_vec,const int* order, const int* neighbor, int Num,
	const vector<int>& vec, int* res);

//algorithm of clustering 
//�����㷨
void clustering(const vector<Record >& vec, int nClus, double cut_tau, int* clus);



