#include<iostream>
#include<vector>
using namespace std;

int main() {
	int n;
	cin >> n;
	vector<vector<int>> matrix(n,vector<int>(n));
	vector<int> res(n);
	vector<int> v(n);
	
	

	//��ʼ��
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = i + j;
		}
	}

	for (int i = 0; i < n; i++) {
		v[i] = i;
	}

	//ƽ���㷨
	for (int i = 0; i < n; i++) {
		res[i] = 0;
		for (int j = 0; j < n; j++) {
			res[i] += matrix[j][i] * v[j];
		}
	}



	return 0;
}