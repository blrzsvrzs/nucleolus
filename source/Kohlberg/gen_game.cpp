#include "gen_game.h"

void type1(vector<double> &v, unsigned int &s, unsigned short int &n) {
	unsigned short int range = 201;
	unsigned short int base = 100 * (n - 2);
	v[s] = rand() % range + base;
	base = 1;
	for (unsigned int i = 0; i < s; i++) {
		unsigned short int S = 1;
		int j = log2(i + 1);
		unsigned int k = i + 1 - pow(2, j);
		while (k > 0) {
			j = log2(k);
			S++;
			k -= pow(2, j);
		}
		if (S == 1)
			v[i] = 0;
		else {
			range = 100 * S;
			v[i] = rand() % range + base;
		}
	}
}

void type2(vector<double> &v, unsigned int &s, unsigned short int &n) {
	unsigned short int range = 50 * n;
	unsigned int j = pow(2, 3);
	v[2] = rand() % range + 1;
	for (unsigned int i = 4; i < s + 1; i++) {
		if (i < (j - 1))
			v[i] = rand() % range + 1;
		else
			j += j;
	}
}

void type3(vector<bool> &v, unsigned int &s, unsigned short int &n) {
	v[s] = true;
	unsigned short int rnd = rand() % 10;
	for (unsigned short int i = 0; i < n; i++) {
		if (rnd < 9)
			v[s - pow(2, i)] = true;
		rnd = rand() % 10;
		for (unsigned short int j = i + 1; j < n; j++) {
			if (rnd < 9)
				v[s - pow(2, i) - pow(2, j)] = true;
			rnd = rand() % 10;
		}
	}
}

void type4(vector<double> &v, unsigned int &s, unsigned short int &n) {
	unsigned int j = pow(2, 3);
	v[2] = rand() % n + 1;
	for (unsigned int i = 4; i < s + 1; i++) {
		if (i < (j - 1))
			v[i] = rand() % n + 1;
		else
			j += j;
	}
}

void type5(vector<bool> &v, unsigned int &s, unsigned short int &n) {
	vector<double> weights(n, 1);
	for (unsigned short int i = 0; i < 5; i++) {
		weights[i] = floor((n - 3) / 2);
	}
	double q = 4 * weights[0] + n - 4;
	wvg(n, weights, q, v, s);
}

void wvg(unsigned short int &n, vector<double> &weights, double &q, vector<bool> &v, unsigned int &s) {
	double sumw = 0;
	for (unsigned short int i = 0; i < n; i++) {
		if (weights[i] > 0)
			sumw += weights[i];
	}
	if (sumw > q) {
		sumw = 0;
		for (unsigned int k = 0; k < s + 1; k++) {
			unsigned int i = k + 1;
			int j = log2(i);
			if (weights[j] != 0)
				sumw += weights[j];
			//A[k][j] = true;
			i -= pow(2, j);
			while (i > 0) {
				j = log2(i);
				if (weights[j] != 0)
					sumw += weights[j];
				//A[k][j] = true;
				i -= pow(2, j);
			}
			if (sumw >= q)
				v[k] = true;
			sumw = 0;
		}
	}
}

void type6(vector<bool> &v, unsigned int &s, unsigned short int &n) {
	vector<double> weights(n, 0);
	double sumw = 0;
	for (unsigned short int i = 0; i < n; i++) {
		weights[i] = rand() % n + 1;
		sumw += weights[i];
	}
	double q = sumw / 2 + 1;
	wvg(n, weights, q, v, s);
}

void type7(vector<double> &v, unsigned short int &n) {
	unsigned int k = pow(2, n - 3) - n + 2;
	vector<unsigned short int> sizes(n - 3, 0);
	vector<unsigned short int> counts(n - 4, 0);
	for (unsigned short int i = 1; i < n - 3; i++) {
		sizes[i] = nchoosek(n - 3, i + 1);
	}
	unsigned int s = pow(2, n - 3) - 1;
	for (unsigned int i = 0; i < s; i++) {
		unsigned short int S = 1;
		int j = log2(i + 1);
		unsigned int l = i + 1 - pow(2, j);
		while (l > 0) {
			j = log2(l);
			S++;
			l -= pow(2, j);
		}
		if (S == 1)
			v[i] = k + 3;
		else {
			v[i] = sizes[S - 2] + counts[S - 2] + 1;
			counts[S - 2]++;
		}
	}
	v[pow(2, n - 3) - 1] = k + 1; v[pow(2, n - 2) - 1] = k + 1; v[pow(2, n - 1) - 1] = k + 1;
	v[pow(2, n - 3) + pow(2, n - 2) - 1] = k + 2; v[pow(2, n - 2) + pow(2, n - 1) - 1] = k + 2; v[pow(2, n - 1) + pow(2, n - 3) - 1] = k + 2;
	v[pow(2, n - 3) + pow(2, n - 2) + pow(2, n - 1) - 1] = k + 3;
	v[pow(2, n) - 2] = n - 2;
}

double nchoosek(unsigned short int n, unsigned short int k) {
	if (k == 0)
		return 1;
	return (n * nchoosek(n - 1, k - 1) / k);
}