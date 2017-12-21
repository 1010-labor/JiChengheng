#pragma once
#pragma once

#include "st.h"



#define  K 40


using namespace std;
void computeDensity(dist** dm, int data_num, dens* density);
double eucDistance(const Record& p1, const Record& p2);
void createDM(const vector<Record>& vec, dist** dm,
	double(*distFunction)(const Record&, const Record&));
void findKnn(dist* dm_line, int left, int right, int k);
void computeDensity(dens* dens_vec, dist** dm, int data_num, int k);
int divisionData(string dataFile, vector<Record>* partVec, int partsNum);
void partMakeSample(vector<Record> &sample_vec, map<int, int> &del_source_map, vector<Record> &part_data, int k);
void makeSample(vector<Record>* part_data_arr, int parts_num,
	vector<Record> &sample_vec, map<int, int>& del_source,int k);
int recordComp(Record &rc1, Record &rc2);
class CompDens {
private:
	dens* dsp;
public:
	CompDens(dens * a) {
		this->dsp = a;
	}
public:bool operator ()(int a, int b) {
	return (dsp + a)->density > (dsp + b)->density;
}

};
