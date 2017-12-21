#include "stdfx.h"
#include <mpi.h>
using namespace std;


typedef double (*DistFunction)(const vector<double>&, const vector<double>&);

bool readData(string dataFile, bool withLable, vector<vector<double>>& dataVec, vector<int>& lableVec);
double EuclideanDistance(const vector<double>& vec1, const vector<double>& vec2);

struct dist {
	int with;
	double distance;
};
