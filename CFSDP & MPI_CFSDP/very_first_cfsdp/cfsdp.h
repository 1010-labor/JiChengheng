#pragma  once

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <string>
#include <cstring>
#include <iterator>
#include <fstream>
#include <sstream>
#include <windows.h>
#include <time.h>
#include "assert.h"
#include <limits>
using namespace std;


typedef double(*DistFunction)(const vector<double>&, const vector<double>&);

bool readData(string dataFile, bool withLable, vector<vector<double>>& dataVec, vector<int>& lableVec);
double EuclideanDistance(const vector<double>& vec1, const vector<double>& vec2);

struct dist {
	int with;
	double distance;
};