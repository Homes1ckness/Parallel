#include<iostream>
#include<vector>
using namespace std;

int main() {
	int n;
	cin >> n;
	long long int sum = 0;
	vector<int> a(n);

	for (int i = 0; i < n; i++) {
		a[i] = i;
	}

	//Æ½·²Ëã·¨
	for (int i = 0; i < n; i++) {
		sum += a[i];
	}

	cout << sum;

	return 0;
}