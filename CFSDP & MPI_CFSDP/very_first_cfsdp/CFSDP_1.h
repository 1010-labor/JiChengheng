#pragma once
#include "sta.h"
using namespace std;

typedef double (*DistFunction)(const vector<double>&, const vector<double>&);



double EuclideanDistance(const vector<double>& vec1, const vector<double>& vec2);
void SetMatrixData(double** dist_matrix, int i, int j, double val);
double GetMatrixData(double** dist_matrix, int i, int j);

//read data from a file
void ReadData(string file_name,vector<vector<double>> &data_set);

void CreateDM(const vector<vector<double>>& data_vec, double** dist_matrix, DistFunction  dist_func);

double FindCutDist(double** dist_matrix,int data_num);

void DensityVec(double** dist_matrix,int data_num,double* density_vec,double cut_dist);

void SortDensity(double* density_vec,int* order,int data_num);


void DetaAndNeighbor(double** dist_matrix, int * order, vector<double> &deta,vector<int> &neiborgh, double &min_deta);


void GamaVec(double *density_vec, const vector<double> &deta, double min_deta, double min_density, double max_deta, double max_density,  vector<double> &gama);

void FindInitialCenters(vector<double> &gama, int n_clus, vector<int>& res_vec);

void AssignClusters(int* order, const vector<int> &neiborgh, vector<int>& res_vec);







