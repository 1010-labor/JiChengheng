#include "st.h"
#include "makesample.h"
using namespace std;

//�Ƚ�����������һά��С
int recordComp(Record &rc1, Record &rc2) {
	return rc1._x < rc2._x;
}

//�˺��������������ŷ�Ͼ���
double eucDistance(const Record& p1, const Record& p2) {
	return sqrt(pow((p1._x - p2._x), 2) + pow((p1._y - p2._y), 2));
}

//�˺�����������󣬸��ݺ������뺯��
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
//�˺��������ܶ�����
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

		lindex = left;//rindexʵָ���һ��Ԫ��
		key_history = rindex;
		//��¼left��Ϊ�ڱ�
		key.with = dm_line[left].with;
		key.dis = dm_line[left].dis;
		
							
		while (lindex < rindex) {
			//����Ѱ�ұ��ڱ�С��Ԫ��
			while (lindex < rindex && key.dis <= (dm_line[rindex]).dis) rindex--;
			//�ҵ����ڱ�С��Ԫ�� �ƶ���������λ�ã�������ǰ��һλ
			if (lindex < rindex)  dm_line[lindex++] = dm_line[rindex];
			//����Ѱ�ұ��ڱ����Ԫ��
			while (lindex < rindex&&key.dis>=(dm_line[lindex]).dis) lindex++;
			//�ҵ����ڱ����Ԫ�� �ƶ���������λ�ã�����������һλ
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




//�˺�����ÿ���Ӽ����г����������Ӽ��ĳ������Լ��Ӽ�ӳ���


//�˺����Ի��ֺõĸ��������ݽ��г���
void makeSample(vector<Record>* part_data_arr, int parts_num, vector<Record> &sample_vec, map<int, int>& del_source,int k) {

	for (int i = 0; i < parts_num; i++)
	{
		partMakeSample(sample_vec, del_source, part_data_arr[i], k);

	}
}



//�˺����Զ�ά���ݼ����д��߽绮��
int divisionData(string dataFile, vector<Record>* partVec, int partsNum) {
	ifstream ifs(dataFile);
	int dataNum=0;
	if (ifs) {
		list<Record> dataVec; //you sort fang fa
		stringstream ss;
		double val;
		string line;
		int index = 0;
		//���ļ������dataVec
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
		//һά����
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

	//����������
	dist** dm;
	dm = new dist*[part_data_num];
	for (int i = 0; i < part_data_num; i++) *(dm + i) = new dist[part_data_num];
	createDM(part_data, dm, eucDistance);
	//�����ܶȲ�����
	dens* dens_vec = new dens[part_data_num];
	computeDensity(dens_vec, dm, part_data_num, K);
	CompDens comp(dens_vec);
	int* densDecend = new int[part_data_num];
	for (int i = 0; i < part_data_num; i++)
	{
		densDecend[i] = i;
	}
	sort(densDecend, densDecend + part_data_num, comp);
	//��densDecend�洢�ܶȽ����±�

	//����
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

