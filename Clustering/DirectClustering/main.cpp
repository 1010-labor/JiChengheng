#include "drc.h"
#define  DATA "d:\\data.dat"
#define CLUSNUM 4

struct Node {
	int data;
	Node* next;
};

bool checkLoop(Node* head) {
	Node *a, *b;
	a = b = head; a = a->next;
	if((*a).next == NULL);
	while (true)
	{
		if (a->next == NULL) {
			return false;
		}
	}


}

void main(){
	int N;
	int sum, remain;
	sum = remain = 0;
	cin >> N;
	vector<int> arrr;
	arrr.push_back(1);
	int *arr = new int[N];
	for (int i = 0; i < N;i++)
	{
		cin >> arr[i];
		sum += arr[i];
	}
	remain = sum;
	
	sort(arr, arr + N);
	for (int i = N-1; i >1;i--)
	{
		remain -= arr[i];
		sum += remain;
	}
	cout<<sum;

}