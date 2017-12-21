#pragma once

#include "st.h"

struct Dist {
	int with;
	double distance;
public:
	Dist(int w, double d) :with(w), distance(d) {}
	Dist() {}

};

bool distCmp(Dist a, Dist b);



void findKnn(Dist* dm_line, int left, int right, int k);
void readData(double* x, double* y, int num, string filePath);
void sortDist(Dist* dm_line, int left, int right);
void getRhoAndSm(double* x, double* y, int num, int k, int shortest, double* rho, Dist** sm);
void getDetaAndNeighbor(double* x, double* y, int num, int shortest, Dist** sm, int* seq_of_rho, double* rho, double* deta, int* neighbor);
