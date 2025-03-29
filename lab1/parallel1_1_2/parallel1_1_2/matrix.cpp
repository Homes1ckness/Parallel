#include<iostream>
#include<vector>
using namespace std;

int main() {
	long long int n;
	cin >> n;
	vector<vector<int>> matrix(n, vector<int>(n));
	vector<int> res(n);
	vector<int> v(n);



	//初始化
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = i + j;
		}
	}

	for (int i = 0; i < n; i++) {
		v[i] = i;
	}

	//cache优化
	for (int i = 0; i < n; i++) {
		res[i] = 0;
	}

	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			res[i] += matrix[j][i] * v[j];
		}
	}



	return 0;
}