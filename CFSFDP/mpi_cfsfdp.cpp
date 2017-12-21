/*************************************************************************
    > File Name: mpi_cfsfdp.cpp
    > Author: Zhang Jiafeng
    > Mail: 761369516@qq.com 
    > Created Time: 2017-05-03
 ************************************************************************/

#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<mpi.h>
#include<algorithm>
#include<math.h>
#include<time.h>
//后面待_t的都是数据类型，它是一个与机器相关的unsigned类型，其大小足以保证存储内存中对象的大小。经常用来存储数组大小，数组下表，因为它足够大，大于等于地址线的位数，所以一个指针可以被安全地放进为size_t类型
size_t dataNum, dim, K, partSize;
double threshold = 0.10;

using namespace std;

struct dataIndex
{
    int index;
    double data;
};

void readData(const char * filename, vector<vector<double> > & dataVec)
{
    stringstream ss;
    double val;
    ifstream in(filename);
    if(in)
    {
        dataVec.clear();
        string line;
        vector<double> subvec;
        while(getline(in,line))
        {
            subvec.clear();
            ss.clear();
            ss << line;
            while(ss >> val)
                subvec.push_back(val);
            dataVec.push_back(subvec);
        }
    }
    else
        cerr << filename << "doesn't exist!\n";
}

double euclideanDistance(const vector<double> &vec1, const vector<double> &vec2)
{
    double sum = 0.0;
    if(vec1.size() != vec2.size())
    {
        cout << "euclideanDistance function error: two parameter not equal"<<endl;
        exit(1);
    }
    for(size_t i=0; i < vec1.size(); i++)
    {
        sum += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }
    return sqrt(sum);
}

void KNN(vector<dataIndex> &vec1)
{
    //头文件: <algorithm>，对给定范围[first,last)内的元素进行重新布置.方法是,nth位置的元素放置的值就是把所有元素排序后在nth位置的值.把所有不大于nth的值放到nth的前面,把所有不小于nth的值放到nth后面.
    nth_element(vec1.begin(), vec1.begin()+K, vec1.end(),[] (const dataIndex &a, const dataIndex &b)
                { return a.data < b.data;});//部分排序，前K小的元素放在前面，并不保证前K各元素有序。

}

void getDensity(const vector<vector<double> > &dataVec, vector<double> &partDensity, vector<vector<dataIndex> > &KNN_part_table, int rank)
{
    vector<dataIndex>  distVec(dataNum-1);
    for(size_t i = 0; i < partSize; i++)
    {    
        size_t jpre = 0;
        double tmp = 0.0;
        for(size_t j = 0; j < dataNum; j++)
        {
            if(rank*partSize+i == j)//防止对本身求距离，jpre是在rank*partSize+i==j之后修正的j
                continue;       
            distVec[jpre].data = euclideanDistance(dataVec[rank*partSize+i],dataVec[j]);
            distVec[jpre].index = j;
            jpre++;
        }
        KNN(distVec);
        for(size_t t = 0; t < K; t++)
        {
            tmp += distVec[t].data;
        }
        nth_element(distVec.begin(), distVec.begin()+10, distVec.begin()+K,[] (const dataIndex &a, const dataIndex &b) { return a.data < b.data;});
        partDensity[i] = dataNum/tmp;
        for(size_t m = 0; m < 10; m++)
        {
            KNN_part_table[i][m].data = distVec[m].data;
            KNN_part_table[i][m].index = distVec[m].index;
        }
    }
}

void getdelta(const vector<double> &density, const vector<vector<double> >&dataVec, const vector<vector<dataIndex> > &shortTable, vector<double> &partDelta,vector<int> &partdeltaIndex, vector<int> &densityIndex, int rank)
{
    //对density增广后进行排序
    vector<dataIndex> enlargeDensity(dataNum);
    for(size_t t = 0; t < dataNum; t++)
    {
        enlargeDensity[t].index = t;
        enlargeDensity[t].data = density[t];
    }
    sort(enlargeDensity.begin(),enlargeDensity.end(),[](const dataIndex &a, const dataIndex &b){return a.data > b.data;});
    for(size_t i = 0; i < dataNum; i++)
        densityIndex[i] = enlargeDensity[i].index;
    for(size_t k = 0; k < partSize; k++)
    {
        size_t i = rank * partSize + k;
        double tmp = 0.0;
        bool flag = false;
        double max = 0.0;
        int index = 0;
        for(size_t j = 0; j < 10; j++)
        {
            if(density[i] < density[shortTable[k][j].index])
            {
                if(tmp == 0.0 || tmp > shortTable[k][j].data)
                {  
                    tmp = shortTable[k][j].data;
                    index = shortTable[k][j].index;
                }
                flag = true;
            }
        }
        if(!flag)
        {
            for(size_t m = 0; density[i] < enlargeDensity[m].data; m++)
            {
                double tmpDist = euclideanDistance(dataVec[i],dataVec[enlargeDensity[m].index]);
                if(tmp == 0.0 || tmp > tmpDist)
                {
                    tmp = tmpDist;
                    index = enlargeDensity[m].index;
                }
                flag = true;
            }
        }
        if(!flag)
        {
            vector<double> vec1(dataNum);
            for(size_t n = 0; n < dataNum; n++)
            {
                vec1[n] = euclideanDistance(dataVec[i],dataVec[n]);
            }
            tmp = *max_element(vec1.begin(), vec1.end());
            index = i;
            flag = true;
        }
        partDelta[k] = tmp;
        partdeltaIndex[k] = index;
    }
}

void assignCluster(const vector<int> &densityIndex, const vector<int> &deltaIndex, const vector<double> &gama, vector<int> &clusterId)
{
    vector<dataIndex> enlargeGama(dataNum);
    for(size_t i = 0; i < dataNum; i++)
    {
        enlargeGama[i].data = gama[i];
        enlargeGama[i].index = i;
    }
    //对gama排序，找出大于阈值的点，即判断有多少个聚类
    sort(enlargeGama.begin(),enlargeGama.end(),[](const dataIndex &a, const dataIndex &b) { return a.data > b.data;});
    for(size_t i = 0; i < dataNum; i++)
        if(enlargeGama[i].data > threshold)
            clusterId[enlargeGama[i].index] = i+1;
        else
            break;
    for(size_t i = 0; i< dataNum; i++)
    {
        size_t index = densityIndex[i];
        if(clusterId[index] == 0 )
            clusterId[index] = clusterId[deltaIndex[index]];
    }
}

int main(int argc, char **argv)
{
    int32_t myid,size;
    const char* filename = "data6.dat";
    const char* RES = "result.dat";
    vector<vector<double> > dataVec;
    readData(filename, dataVec);
    dataNum = dataVec.size();
    dim = dataVec[0].size();
    K =(int) (dataNum * 0.015);
    sort(dataVec.begin(), dataVec.end(),[](const vector<double> &a, const vector<double> &b){ return a[0] < b[0];});//将数据排序，以便数据划分
    clock_t start,end;
    start = clock();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    partSize = ceil((double)dataNum/size);
//计算密度
    vector<double> density(dataNum);
    vector<double> partDensity(partSize);
    vector<vector<dataIndex> > partShortTable(partSize, vector<dataIndex> (10));
    if(myid != 0)
    {
        getDensity(dataVec, partDensity, partShortTable, myid);
        MPI_Send(&partDensity[0], partDensity.size(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); 
    }
    if(myid == 0)
    {
        getDensity(dataVec, density, partShortTable, myid);
        MPI_Status status;
        for(size_t i = 1; i < size; i++)
        {
            MPI_Recv(&density[i*partSize], partSize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
        }
    }
//广播密度
    MPI_Bcast(&density[0], dataNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
//计算deta值   
    vector<double> delta(dataNum);
    vector<int> deltaIndex(dataNum);
    vector<double> partDelta(partSize);
    vector<int> partdeltaIndex(partSize);
    vector<int> densityIndex(dataNum);
    if(myid !=0)
    {
        getdelta(density, dataVec, partShortTable,partDelta, partdeltaIndex, densityIndex, myid);
        MPI_Send(&partDelta[0], partSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD); 
        MPI_Send(&partdeltaIndex[0], partSize, MPI_INT, 0, 1, MPI_COMM_WORLD); 
    }
    if(myid == 0)
    {
        getdelta(density, dataVec, partShortTable, delta, deltaIndex,densityIndex, myid);
        MPI_Status status;
        for(size_t i = 1; i < size; i++)
        {
            MPI_Recv(&delta[i*partSize], partSize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&deltaIndex[i*partSize], partSize, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
        }
    }
//广播delta值
    MPI_Bcast(&delta[0], dataNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&deltaIndex[0], dataNum, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
//计算gama    
//gama[i]=((density[i]-density[min])*(delta[i]-delta[min]))/((density[max]-density[min])*(delta[max]-delta[min]))
    vector<double> gama(dataNum);
    vector<double> partGama(partSize);
    double min_den,max_den,min_del,max_del;
    min_den = density[densityIndex[dataNum-1]];
    max_den = density[densityIndex[0]];
    min_del = *min_element(delta.begin(),delta.end());
    max_del = *max_element(delta.begin(),delta.end());
    if(myid != 0)
    {
        for(size_t i = 0; i < partSize; i++)
        {
            size_t j = myid * partSize + i;
            partGama[i] = ((density[j]-min_den)*(delta[j]-min_del))/((max_den-min_den)*(max_del-min_del));
        }
        MPI_Send(&partGama[0], partSize, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    if(myid == 0)
    {
        for(size_t i = 0; i < partSize; i++)
        {
            size_t j = myid * partSize + i;
            gama[i] = ((density[j]-min_den)*(delta[j]-min_del))/((max_den-min_den)*(max_del-min_del));
        }
        MPI_Status status;
        for(size_t i = 1; i < size; i++)
            MPI_Recv(&gama[i*partSize], partSize, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
        
    }
//广播gama值
//    MPI_Bcast(&gama[0], dataNum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 0)
    {
        vector<int> clusterId(dataNum);
        assignCluster(densityIndex, deltaIndex, gama, clusterId);
        ofstream out(RES);
        for(size_t i = 0; i < dataNum; i++)
        {
            for(size_t j = 0; j < dim; j++)
                out << dataVec[i][j] << "\t";
            out<<clusterId[i]<<"\n";
        }
        end = clock();
        cout<<"runtime: "<<(double)(end-start)/CLOCKS_PER_SEC<<endl;
    }
    MPI_Finalize();
}

