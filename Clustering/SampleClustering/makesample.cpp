#include "st.h"
#include "makesample.h"
using namespace std;

//比较两个样本的一维大小
int recordComp(Record &rc1, Record &rc2) {
	return rc1._x < rc2._x;
}

//此函数计算两个点的欧氏距离
double eucDistance(const Record& p1, const Record& p2) {
	return sqrt(pow((p1._x - p2._x), 2) + pow((p1._y - p2._y), 2));
}

//此函创建距离矩阵，根据函数距离函数
void createDM(const vector<Record>& vec, dist** dm,
	double(*distFunction)(const Record&, const Record&)) {
	int data_num = vec.size();
	for (int i = 0; i < data_num; i++)
	{
		dist tmp;
		tmp.with = i;
		tmp.dis =0;
		dm[i][i] = tmp;

		for (int j = i+1; j < data_num; j++) {
					
			tmp.with = j;
			tmp.dis = distFunction(vec[i], vec[j]);
			dm[j][i] = dm[i][j] = tmp;
		}
	}






}

int distComp(dist &a, dist &b) {
	return a.dis < b.dis;


}
//此函数计算密度序列
void computeDensity(dens* dens_vec, dist** dm, int data_num, int k) {

	for (int i = 0; i < data_num; i++)
	{
		findKnn(*(dm + i), 0, data_num -1, k+1 );
		sort(*(dm + i), *(dm + i) + k + 1, distComp);
	}

	for (int i = 0; i < data_num; i++) {
		(dens_vec + i)->index = i;
		double tmp = 0;
		for (int j = 1; j <= k; j++) {
			tmp += 1 / (*(dm + i) + j)->dis;
		}
		(dens_vec + i)->density = tmp;
	}


}

void findKnn(dist* dm_line, int left, int right, int k) {

	int current_find = right - left + 1;
	dist key;
	int lindex;
	int rindex = right;
	int key_history;


	while (current_find > k+1) {

		lindex = left;//rindex实指最后一个元素
		key_history = rindex;
		//记录left作为哨兵
		key.with = dm_line[left].with;
		key.dis = dm_line[left].dis;
		
							
		while (lindex < rindex) {
			//逆向寻找比哨兵小的元素
			while (lindex < rindex && key.dis <= (dm_line[rindex]).dis) rindex--;
			//找到比哨兵小的元素 移动到左索引位置，左索引前进一位
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			//正向寻找比哨兵大的元素
			while (lindex < rindex&&key.dis>=(dm_line[lindex]).dis) lindex++;
			//找到比哨兵大的元素 移动到右索引位置，右索引后退一位
			if (lindex < rindex) 	dm_line[rindex--] = dm_line[lindex];
		}
		assert(lindex == rindex);
		dm_line[lindex] = key;
		current_find = rindex - left+1;
	}
	if (current_find == k+1|| current_find == k) return;
	else {
		findKnn(dm_line, rindex+1, key_history, k - current_find);
	}

}




//此函数对每个子集进行抽样，返回子集的抽样，以及子集映射表


//此函数对划分好的各部分数据进行抽样
void makeSample(vector<Record>* part_data_arr, int parts_num, vector<Record> &sample_vec, map<int, int>& del_source,int k) {

	for (int i = 0; i < parts_num; i++)
	{
		partMakeSample(sample_vec, del_source, part_data_arr[i], k);

	}
}



//此函数对二维数据集进行带边界划分
int divisionData(string dataFile, vector<Record>* partVec, int partsNum) {
	ifstream ifs(dataFile);
	int dataNum=0;
	if (ifs) {
		list<Record> dataVec; //you sort fang fa
		stringstream ss;
		double val;
		string line;
		int index = 0;
		//读文件，填充dataVec
		while (getline(ifs, line)) {
			Record rc;
			rc.index = index++;
			ss.clear();
			ss << line;
			ss >> val;
			rc._x = val;
			ss >> val;
			rc._y = val;

			dataVec.push_back(rc);
		}
		//一维排序
		dataVec.sort(recordComp);

		dataNum= dataVec.size();
		int border = dataNum /40;
		int partsize = ceil(dataNum / partsNum);
		list<Record>::iterator it;
		int i = 0;
		for (int j = 0; j < partsNum; j++)
		{
			
			for (i = 0, it = dataVec.begin(); i < j* partsize - border; i++, it++);
			for (; i < (j+1)*partsize + border&&i < dataNum; i++, it++)
			{
				
				Record temp;
				temp.isBorder = i < (j + 1)*partsize&&i >= j*partsize ? 0 :1;
				temp.index = it->index;
				temp._x = it->_x;
				temp._y = it->_y;
				partVec[j].push_back(temp);

			}
		}
		
		
	}
	return dataNum;
}

void partMakeSample(vector<Record> &sample_vec, map<int, int> &del_source_map, vector<Record> &part_data, int k) {
	int part_data_num = part_data.size();

	//计算距离矩阵
	dist** dm;
	dm = new dist*[part_data_num];
	for (int i = 0; i < part_data_num; i++) *(dm + i) = new dist[part_data_num];
	createDM(part_data, dm, eucDistance);
	//计算密度并排序
	dens* dens_vec = new dens[part_data_num];
	computeDensity(dens_vec, dm, part_data_num, K);
	CompDens comp(dens_vec);
	int* densDecend = new int[part_data_num];
	for (int i = 0; i < part_data_num; i++)
	{
		densDecend[i] = i;
	}
	sort(densDecend, densDecend + part_data_num, comp);
	//用densDecend存储密度降序下标

	//抽样
	list<int> HH;
	for (int i = 0; i < part_data_num; i++)
	{
		HH.push_back(densDecend[i]);
	}
	delete[] densDecend;
	list<int>::iterator  it = HH.begin();

	while (it != HH.end())
	{
		if (!part_data[*it].isBorder) {
			sample_vec.push_back(part_data[*it]);
			int with;
			for (int i = 1; i <= k; i++)
			{
				with = (*(dm + *it) + i)->with;
				if ((dens_vec + *it)->density >= (dens_vec + with)->density)
				{
					list<int>::iterator  it_foward = it;
					while (++it_foward != HH.end()) {
						if (*it_foward == with) {
							if (!part_data[*it_foward].isBorder) {
								int s = part_data[with].index;
								int p = part_data[*it].index;
								del_source_map.insert(make_pair(s, p));

							}
							HH.remove(with);
							break;
						}
					}
				}
			}
			it++;
		}
		else
		{
			it = HH.erase(it);
		}

	}


	for (int i = 0; i < part_data_num; i++) delete[] *(dm + i);
	delete[] dens_vec;

};

