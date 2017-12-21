#include "makesample.h"
#include "clustering.h"
#define  DATA "d:\\part_4.dat"
//#define  DATA "D:\\data.dat"
#define CLUSNUM 4
#define PARTNUM 8
#define KK 10
void indKnn(vector<int>& line, int left, int right, int k);

void  sampleClustering() {
	clock_t start, finish;
	double totaltime;
	double t1, t2;
	start = clock();
	

	vector<Record> part_vecs[PARTNUM];
	int data_num = divisionData(DATA, part_vecs, PARTNUM);
	
	cout << "part_num     =====    " <<PARTNUM<< endl;
	cout << "origin size  =====    " << data_num << endl;
	vector<Record> sample;
	map<int, int> mm;
	int k = KK;
	cout << "k            =====    " << k<<endl;
	makeSample(part_vecs, PARTNUM, sample, mm, k);

	finish = clock();

	t1 = double(finish - start) / CLOCKS_PER_SEC;
	cout << "t1  ============    " <<t1 << "  sencods"<<endl;
	cout <<  "sample size====    " << sample.size() << endl;
	start = clock();

	int* clus = new int[data_num];
	memset(clus, 0, sizeof(int)*data_num);
	clustering(sample, CLUSNUM, 0.02, clus);

	ofstream sam_file("d:\\sample.txt");
	ofstream sam_res("d:\\sample_res.txt");


	for (int i = 0; i < data_num; i++)
	{
		if (clus[i] == 0)
		{
			
			clus[i] = clus[(mm.find(i))->second];
		}
	}
	for (int i = 0; i < data_num; i++)
	{
		sam_res << i << " " << clus[i] << " " << 0 << endl;
	}
	sam_file.close();
	sam_res.close();


	finish = clock();
	t2 = double(finish - start) / CLOCKS_PER_SEC;
	cout << "t2  ============    " << t2 << "  sencods" << endl;

	cout << "total time =====    " << t1 + t2 << endl;



}
void directClustering() {
	clock_t start, finish;
	double totaltime;
	start = clock();


	vector<Record> part_vecs[4];
	int data_num = divisionData(DATA, part_vecs, 4);
	part_vecs[1];


	vector<Record> sample;
	map<int, int> mm;
	int k = 6;
	makeSample(part_vecs, CLUSNUM, sample, mm, k);

	finish = clock();
	cout << " make sample time ====" << (finish - start) / CLOCKS_PER_SEC << "  sencods";
	start = clock();

	int* clus = new int[data_num];
	memset(clus, 0, sizeof(int)*data_num);
	clustering(sample, 7, 0.01, clus);

	ofstream sam_file("d:\\sample.txt");
	ofstream sam_res("d:\\sample_res.txt");


	for (int i = 0; i < data_num; i++)
	{
		if (clus[i] == 0)
		{
			clus[i] = clus[(mm.find(i))->second];
		}
	}
	for (int i = 0; i < data_num; i++)
	{
		sam_res << i << " " << clus[i] << " " << 0 << endl;
	}
	sam_file.close();
	sam_res.close();

	finish = clock();
	cout << "time ====" << (finish - start) / CLOCKS_PER_SEC << "  sencods";

}
void main() {

	sampleClustering();


	/*
		set<int> sample_index;
		vector<Record>::iterator it = sample.begin();
		while (it != sample.end())
		{

			sample_index.insert(it->index);
			it++;
		}
		map<int, int>::iterator it2=mm.begin();
		while (it2!=mm.end())
		{

			assert(sample_index.find(it2->first)==sample_index.end());
			assert(sample_index.find(it2->second) != sample_index.end());
			it2++;
		}*/
		/*
			vector<int> line;
			line.reserve(50000);
			srand(time(0));
			int nn0=0;
			for (int i = 0; i < 50000;i++)
			{
				line.push_back(rand() % 5000);
				if (line[i]==0)
				{
					nn0++;

				}
			}
			cout << nn0<<"              find "<<endl;

			indKnn(line, 0, 49999, nn0);

			//sort(line.begin(), line.end());

			int K_max = line[0];
			for (int i = 0; i < nn0; i++)
			{
				cout << line[i]<<" "<<endl;
				K_max = line[i] > K_max ? line[i] : K_max;

			}
			for (int i =nn0; i < 50000; i++)
			{
				if (line[i]<K_max)
				{
					cout << "fuck"<<line[i]<<"   "<<K_max<<endl;
				}

			}


			finish = clock();
			totaltime = (double)(finish - start) ;
			cout << "\n此程序的运行时间为" << totaltime << "毫秒！" << endl;*/

}

void indKnn(vector<int>& line, int left, int right, int k) {

	int current_find = right - left + 1;
	int key;
	int lindex;
	int rindex = right;
	int key_history;

	while (current_find > k + 1) {

		lindex = left;//rindex实指最后一个元素
		key_history = rindex;
		//记录left作为哨兵
		key = line[left];


		while (lindex < rindex) {
			//逆向寻找比哨兵小的元素
			while (lindex < rindex && key <= line[rindex]) rindex--;
			//找到比哨兵小的元素 移动到左索引位置，左索引前进一位
			if (lindex < rindex)  line[lindex++] = line[rindex];
			//正向寻找比哨兵大的元素
			while (lindex < rindex&&key >= line[lindex]) lindex++;
			//找到比哨兵大的元素 移动到右索引位置，右索引后退一位
			if (lindex < rindex) 	line[rindex--] = line[lindex];
		}
		assert(lindex == rindex);
		line[lindex] = key;
		current_find = rindex - left + 1;
	}
	if (current_find == k + 1 || current_find == k) return;
	else {
		indKnn(line, rindex + 1, key_history, k - current_find);
	}

}