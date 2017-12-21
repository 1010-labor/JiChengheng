#include "st.h"
#include "func.h"

#define FILEPATH "d:\\data.dat"
#define DATANUM 3200
#define K 32
#define SHORTEST 30
double rho[DATANUM];
using namespace std;
bool cmpp(int a, int b) {
	return rho[a] > rho[b];

}

int main() {

	double x[DATANUM], y[DATANUM];
	readData(x, y, DATANUM, FILEPATH);
	Dist** sm = new Dist*[DATANUM];
	for (int i = 0; i < DATANUM; i++) {
		*(sm + i) = new Dist[SHORTEST];
	}
	
	getRhoAndSm(x, y, DATANUM, K, SHORTEST, rho, sm);

	fstream f("d:\\rho.txt");
	for (int i = 0; i < DATANUM; i++) {
		f << setw(10) << x[i] << setw(10) << y[i] <<" "<< rho[i]*100<<endl;
	}
	f.flush();
	f.close();

	int seq_of_rho[DATANUM];
	for (int i = 0; i < DATANUM; i++) {
		seq_of_rho[i] = i;
	}
	sort(seq_of_rho, seq_of_rho +DATANUM, cmpp);
	cout << "------------------------------------" << endl;
	for (int i = 0; i < 50;i++)
	{
		cout << rho[seq_of_rho[i]] << endl;
	}
	int* neighbor = new int[DATANUM];
	double* deta = new double[DATANUM];
	getDetaAndNeighbor(x, y, DATANUM, SHORTEST, sm, seq_of_rho, rho, deta, neighbor);

	cout << "------------------------------------" << endl;
	for (int i = 0; i < 50; i++) {
		cout << neighbor[i] << " ";
		if (i % 100 == 0) {
			cout << endl;
		}

	}
}
