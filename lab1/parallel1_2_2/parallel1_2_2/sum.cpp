#include<iostream>
#include<vector>
using namespace std;

int main() {
	int n;
	cin >> n;

	long long int sum = 0;
	long long int sum1 = 0;
	long long int sum2 = 0;
	vector<int> a(n);

	for (int i = 0; i < n; i++) {
		a[i] = i;
	}

	//cache优化  多链路式

	for (int i = 0; i < n; i += 2) {
		sum1 += a[i];
		sum2 += a[i + 1];
	}

	sum = sum1 + sum2;

	cout << sum;

	return 0;
}