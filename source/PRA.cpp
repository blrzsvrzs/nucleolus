/*
*    Nucleolus
*    PRA.cpp
*    Purpose: finding the nucleolus of a cooperative game using
*             Potters et al. (1996) - Computing the nucleolus by solving
*			  a prolonged simplex algorithm
*
*    @author Marton Benedek
*    @version 1.2 16/08/2019
*
*    This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program. If not, see:
*    <https://github.com/blrzsvrzs/nucleolus>.
*/

#include "PRA.h"
#include "gen_game.h"

int main() {
	unsigned short int n = 0;
	unsigned short int type = 0;
	unsigned int seed = 0;
	bool disp = false;
	bool memo = false;
	bool nlsu = false;
	ifstream inp;
	cout << "Reading the input...";
	inp.open("input.txt");
	inp >> n >> type >> seed >> disp >> memo >> nlsu;
	inp.close();
	cout << "done!" << endl;
	if (seed == 0)
		seed = GetTickCount();
	srand(seed);
	unsigned int s = pow(2, n) - 2;
	vector<double> x(n, 0);
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	unsigned short int iter = 0;
	unsigned int piv = 0;
	double t = 0;
	if (type == 3 || type == 5) {
		cout << "Generating game...";
		vector<bool> v(s + 1, 0);
		if (type == 3)
			type3(v, s, n);
		else if (type == 5)
			type5(v, s, n);
		cout << "done!" << endl;
		cout << "Running PRA..." << endl;
		double t1 = cpuTime();
		PRA_sg(v, n, s, unsettled, disp, prec, iter, piv, x, t, t1);
	}
	else {
		vector<double> v(s + 1, 0);
		if (type == 0) {
			cout << "Loading game...";
			inp.open("v.txt");
			for (unsigned int i = 0; i < s + 1; i++)
				inp >> v[i];
			inp.close();
			cout << "done!" << endl;
		}
		else {
			cout << "Generating game...";
			if (type == 1)
				type1(v, s, n);
			else if (type == 2)
				type2(v, s, n);
			else if (type == 4)
				type4(v, s, n);
			cout << "done!" << endl;
		}
		cout << "Running PRA..." << endl;
		double t1 = cpuTime();
		PRA(v, n, s, unsettled, disp, prec, iter, piv, x, t, t1);
	}
	ofstream res;
	res.open("results.txt", ofstream::out | ofstream::trunc);
	res << seed << endl << t << endl << iter << endl << piv << endl;
	for (unsigned int i = 0; i < n; i++)
		res << fixed << setprecision(17) << x[i] << endl;
	res.close();
	cout << "Press 0 then Enter to quit: ";
	double quit;
	cin >> quit;
	cin.get();
	return 0;
}

void PRA(vector<double> &v, unsigned short int &n, unsigned int &s, vector<bool> &unsettled, bool &disp, double &prec, unsigned short int &iter, unsigned int &piv, vector<double> &x, double &t, double &t1) {
	// find c
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	double c = 0;
	for (unsigned int i = 0; i < s; i++)
		if (v[i] > c)
			c = v[i];
	// create tables
	vector<vector<double>> A(s + 1, vector<double>(n, 0));
	vector<vector<double>> B(s + 1, vector<double>(s, 0));
	vector<double> RHS(s + 1, 0);
	for (unsigned short int i = 0; i < n; i++)
		A[0][i] = 1;
	RHS[0] = v[s];
	vector<unsigned int> basis(s + 1, n + s);
	basis[0] = n - 1;
	vector<bool> basis_ind(n + s + 2, true);
	basis_ind[n + s] = false;
	unsigned short int pb = 0;
	vector<bool> settled_rows(s + 1, false);
	vector<bool> a(n, 0);
	unsigned int settled = 0;
	for (unsigned int i = 0; i < s / 2; i++) {
		de2bi(i, a, n);
		for (unsigned short int j = 0; j < n; j++) {
			if (i == 0)
				basis_ind[j] = false;
			if (a[j]) {
				A[i + 1][j] = -1;
				A[s - i][j] = 1;
			}
		}
		B[i + 1][i] = 1;
		B[i + 1 + s / 2][i + s / 2] = 1;
		RHS[i + 1] = c - v[i];
		RHS[i + 1 + s / 2] = c - v[i + s / 2] + v[s];
	}
	basis_ind[n - 1] = true;
	if (disp && n <= 5) {
		cout << "Starting table:" << endl;
		for (unsigned int i = 0; i < s + 1; i++) {
			for (unsigned int j = 0; j < n; j++) {
				if (A[i][j] >= 0)
					cout << " ";
				cout << A[i][j] << " ";
			}
			for (unsigned int j = 0; j < s; j++) {
				if (B[i][j] >= 0)
					cout << " ";
				cout << B[i][j] << " ";
			}
			if (i == 0)
				cout << " " << 0 << " ";
			else
				cout << " " << 1 << " ";
			cout << " " << RHS[i] << endl;
		}
		cin.get();
	}
	// find pivot row in t-column: find min of RHS
	unsigned int row = 0;
	double min_ratio = DBL_MAX;
	unsigned int min_idx = INT_MAX;
	for (unsigned int i = 1; i < s + 1; i++) {
		if (RHS[i] < prec) {
			row = i;
			break;
		}
		else {
			if (min_ratio > RHS[i]) {
				min_ratio = RHS[i];
				row = i;
			}
		}
	}
	basis_ind[n + row - 1] = false;
	// pivot on t-column: subtract row from every row of A, B and RHS, except first one
	for (unsigned int i = 1; i < s + 1; i++) {
		if (i != row) {
			basis[i] -= s - i + 1;
			for (unsigned short int j = 0; j < n; j++)
				if (A[row][j] > prec || A[row][j] < -prec)
					A[i][j] -= A[row][j];
			B[i][row - 1] = -1;
			if (RHS[i] > prec)
				RHS[i] -= RHS[row];
		}
		else
			basis[i] = n + s + 1;
	}
	piv++;
	if (disp && n <= 5) {
		cout << "New table:" << endl;
		for (unsigned int i = 0; i < s + 1; i++) {
			for (unsigned int j = 0; j < n; j++) {
				if (A[i][j] >= 0)
					cout << " ";
				cout << A[i][j] << " ";
			}
			for (unsigned int j = 0; j < s; j++) {
				if (B[i][j] >= 0)
					cout << " ";
				cout << B[i][j] << " ";
			}
			cout << " " << RHS[i] << endl;
		}
		cout << "Basis:" << endl;
		for (unsigned int i = 0; i < s + 1; i++)
			cout << basis[i] << " ";
		cout << endl;
		cout << "Basis idx:" << endl;
		for (unsigned int i = 0; i < n + s + 2; i++)
			cout << basis_ind[i] << " ";
		cout << endl;
		cin.get();
	}
	solve(row, A, B, RHS, s, n, prec, min_ratio, disp, piv, basis, basis_ind, unsettled, pb, settled_rows, min_idx);
	double eps = RHS[row] - c;
	update(s, unsettled, B, A, row, n, prec, pb, RHS, settled_rows, min_ratio, basis, basis_ind, settled, disp, rank, Arref, J, a);
	if (rank == n)
		settled = s;
	if (disp) {
		cout << "Least core solution:" << endl;
		min_ratio = 0;
		row = s + pb + 1;
		for (unsigned int j = 0; j < n; j++) {
			if (!basis_ind[j]) {
				for (unsigned int i = 0; i < s + pb + 1; i++) {
					if (!settled_rows[i]) {
						if (A[i][j] > prec || A[i][j] < -prec) {
							if (min_ratio == 0 && A[i][j] < 1 + prec && A[i][j] > 1 - prec) {
								row = i;
								min_ratio += 1;
							}
							else {
								if (min_ratio == 1)
									row = s + pb + 1;
							}
						}
					}
				}
				if (row < s + pb + 1) {
					x[j] = RHS[row];
					row = s + pb + 1;
				}
				min_ratio = 0;
			}
		}
		for (unsigned int i = 0; i < s + pb + 1; i++) {
			if (basis[i] < n)
				x[basis[i]] = RHS[i];
		}
		for (unsigned int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << eps << endl;
		if (n <= 5) {
			cout << "New table:" << endl;
			for (unsigned int i = 0; i < s + pb + 1; i++) {
				if (!settled_rows[i]) {
					for (unsigned int j = 0; j < n; j++) {
						if (A[i][j] >= 0)
							cout << " ";
						cout << A[i][j] << " ";
					}
					for (unsigned int j = 0; j < s; j++) {
						if (unsettled[j]) {
							if (B[i][j] >= 0)
								cout << " ";
							cout << B[i][j] << " ";
						}
					}
					cout << " " << RHS[i] << endl;
				}
			}
			cout << "Basis:" << endl;
			for (unsigned int i = 0; i < s + 1 + pb; i++)
				if (!settled_rows[i])
					cout << basis[i] << " ";
			cout << endl;
			cout << "Basis idx:" << endl;
			for (unsigned int i = 0; i < n + s + 2; i++)
				cout << basis_ind[i] << " ";
			cout << endl;
			cin.get();
		}
	}
	iter++;
	while (settled < s)
		iteration(s, pb, settled_rows, unsettled, B, prec, A, row, min_ratio, RHS, n, basis_ind, basis, piv, disp, eps, iter, settled, x, rank, Arref, J, a, min_idx);
	min_ratio = 0;
	row = s + pb + 1;
	for (unsigned int j = 0; j < n; j++) {
		if (!basis_ind[j]) {
			for (unsigned int i = 0; i < s + pb + 1; i++) {
				if (!settled_rows[i]) {
					if (A[i][j] > prec || A[i][j] < -prec) {
						if (min_ratio == 0 && A[i][j] < 1 + prec && A[i][j] > 1 - prec) {
							row = i;
							min_ratio += 1;
						}
						else {
							if (min_ratio == 1)
								row = s + pb + 1;
						}
					}
				}
			}
			if (row < s + pb + 1) {
				x[j] = RHS[row];
				row = s + pb + 1;
			}
			min_ratio = 0;
		}
	}
	for (unsigned int i = 0; i < s + pb + 1; i++) {
		if (basis[i] < n)
			x[basis[i]] = RHS[i];
	}
	t = cpuTime() - t1;
	cout << "PRA finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	return;
}

void PRA_sg(vector<bool> &v, unsigned short int &n, unsigned int &s, vector<bool> &unsettled, bool &disp, double &prec, unsigned short int &iter, unsigned int &piv, vector<double> &x, double &t, double &t1) {
	// find c
	double c = 0;
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	for (unsigned int i = 0; i < s; i++)
		if (v[i] > c)
			c = v[i];
	// create tables
	vector<vector<double>> A(s + 1, vector<double>(n, 0));
	vector<vector<double>> B(s + 1, vector<double>(s, 0));
	vector<double> RHS(s + 1, 0);
	for (unsigned short int i = 0; i < n; i++)
		A[0][i] = 1;
	RHS[0] = v[s];
	vector<unsigned int> basis(s + 1, n + s);
	basis[0] = n - 1;
	vector<bool> basis_ind(n + s + 2, true);
	basis_ind[n + s] = false;
	unsigned short int pb = 0;
	vector<bool> settled_rows(s + 1, false);
	vector<bool> a(n, 0);
	unsigned int settled = 0;
	for (unsigned int i = 0; i < s / 2; i++) {
		de2bi(i, a, n);
		for (unsigned short int j = 0; j < n; j++) {
			if (i == 0)
				basis_ind[j] = false;
			if (a[j]) {
				A[i + 1][j] = -1;
				A[s - i][j] = 1;
			}
		}
		B[i + 1][i] = 1;
		B[i + 1 + s / 2][i + s / 2] = 1;
		RHS[i + 1] = c - v[i];
		RHS[i + 1 + s / 2] = c - v[i + s / 2] + v[s];
	}
	basis_ind[n - 1] = true;
	if (disp && n <= 5) {
		cout << "Starting table:" << endl;
		for (unsigned int i = 0; i < s + 1; i++) {
			for (unsigned int j = 0; j < n; j++) {
				if (A[i][j] >= 0)
					cout << " ";
				cout << A[i][j] << " ";
			}
			for (unsigned int j = 0; j < s; j++) {
				if (B[i][j] >= 0)
					cout << " ";
				cout << B[i][j] << " ";
			}
			if (i == 0)
				cout << " " << 0 << " ";
			else
				cout << " " << 1 << " ";
			cout << " " << RHS[i] << endl;
		}
		cin.get();
	}
	// find pivot row in t-column: find min of RHS
	unsigned int row = 0;
	double min_ratio = DBL_MAX;
	unsigned int min_idx = INT_MAX;
	for (unsigned int i = 1; i < s + 1; i++) {
		if (RHS[i] < prec) {
			row = i;
			break;
		}
		else {
			if (min_ratio > RHS[i]) {
				min_ratio = RHS[i];
				row = i;
			}
		}
	}
	basis_ind[n + row - 1] = false;
	// pivot on t-column: subtract row from every row of A, B and RHS, except first one
	for (unsigned int i = 1; i < s + 1; i++) {
		if (i != row) {
			basis[i] -= s - i + 1;
			for (unsigned short int j = 0; j < n; j++)
				if (A[row][j] > prec || A[row][j] < -prec)
					A[i][j] -= A[row][j];
			B[i][row - 1] = -1;
			if (RHS[i] > prec)
				RHS[i] -= RHS[row];
		}
		else
			basis[i] = n + s + 1;
	}
	piv++;
	if (disp && n <= 5) {
		cout << "New table:" << endl;
		for (unsigned int i = 0; i < s + 1; i++) {
			for (unsigned int j = 0; j < n; j++) {
				if (A[i][j] >= 0)
					cout << " ";
				cout << A[i][j] << " ";
			}
			for (unsigned int j = 0; j < s; j++) {
				if (B[i][j] >= 0)
					cout << " ";
				cout << B[i][j] << " ";
			}
			cout << " " << RHS[i] << endl;
		}
		cout << "Basis:" << endl;
		for (unsigned int i = 0; i < s + 1; i++)
			cout << basis[i] << " ";
		cout << endl;
		cout << "Basis idx:" << endl;
		for (unsigned int i = 0; i < n + s + 2; i++)
			cout << basis_ind[i] << " ";
		cout << endl;
		cin.get();
	}
	solve(row, A, B, RHS, s, n, prec, min_ratio, disp, piv, basis, basis_ind, unsettled, pb, settled_rows, min_idx);
	double eps = RHS[row] - c;
	update(s, unsettled, B, A, row, n, prec, pb, RHS, settled_rows, min_ratio, basis, basis_ind, settled, disp, rank, Arref, J, a);
	if (rank == n)
		settled = s;
	if (disp) {
		cout << "Least core solution:" << endl;
		min_ratio = 0;
		row = s + pb + 1;
		for (unsigned int j = 0; j < n; j++) {
			if (!basis_ind[j]) {
				for (unsigned int i = 0; i < s + pb + 1; i++) {
					if (!settled_rows[i]) {
						if (A[i][j] > prec || A[i][j] < -prec) {
							if (min_ratio == 0 && A[i][j] < 1 + prec && A[i][j] > 1 - prec) {
								row = i;
								min_ratio += 1;
							}
							else {
								if (min_ratio == 1)
									row = s + pb + 1;
							}
						}
					}
				}
				if (row < s + pb + 1) {
					x[j] = RHS[row];
					row = s + pb + 1;
				}
				min_ratio = 0;
			}
		}
		for (unsigned int i = 0; i < s + pb + 1; i++) {
			if (basis[i] < n)
				x[basis[i]] = RHS[i];
		}
		for (unsigned int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << eps << endl;
		if (n <= 5) {
			cout << "New table:" << endl;
			for (unsigned int i = 0; i < s + pb + 1; i++) {
				if (!settled_rows[i]) {
					for (unsigned int j = 0; j < n; j++) {
						if (A[i][j] >= 0)
							cout << " ";
						cout << A[i][j] << " ";
					}
					for (unsigned int j = 0; j < s; j++) {
						if (unsettled[j]) {
							if (B[i][j] >= 0)
								cout << " ";
							cout << B[i][j] << " ";
						}
					}
					cout << " " << RHS[i] << endl;
				}
			}
			cout << "Basis:" << endl;
			for (unsigned int i = 0; i < s + 1 + pb; i++)
				if (!settled_rows[i])
					cout << basis[i] << " ";
			cout << endl;
			cout << "Basis idx:" << endl;
			for (unsigned int i = 0; i < n + s + 2; i++)
				cout << basis_ind[i] << " ";
			cout << endl;
			cin.get();
		}
	}
	iter++;
	while (settled < s)
		iteration(s, pb, settled_rows, unsettled, B, prec, A, row, min_ratio, RHS, n, basis_ind, basis, piv, disp, eps, iter, settled, x, rank, Arref, J, a, min_idx);
	min_ratio = 0;
	row = s + pb + 1;
	for (unsigned int j = 0; j < n; j++) {
		if (!basis_ind[j]) {
			for (unsigned int i = 0; i < s + pb + 1; i++) {
				if (!settled_rows[i]) {
					if (A[i][j] > prec || A[i][j] < -prec) {
						if (min_ratio == 0 && A[i][j] < 1 + prec && A[i][j] > 1 - prec) {
							row = i;
							min_ratio += 1;
						}
						else {
							if (min_ratio == 1)
								row = s + pb + 1;
						}
					}
				}
			}
			if (row < s + pb + 1) {
				x[j] = RHS[row];
				row = s + pb + 1;
			}
			min_ratio = 0;
		}
	}
	for (unsigned int i = 0; i < s + pb + 1; i++) {
		if (basis[i] < n)
			x[basis[i]] = RHS[i];
	}
	t = cpuTime() - t1;
	cout << "PRA finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	return;
}

void iteration(unsigned int &s, unsigned short int &pb, vector<bool> &settled_rows, vector<bool> &unsettled, vector<vector<double>> &B, double &prec, vector<vector<double>> &A, unsigned int &row, double &min_ratio, vector<double> &RHS, unsigned short int &n, vector<bool> &basis_ind, vector<unsigned int> &basis, unsigned int &piv, bool &disp, double &eps, unsigned short int &iter, unsigned int &settled, vector<double> &x, unsigned short int &rank, vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &a, unsigned int &min_idx) {
	min_idx = INT_MAX;
	vector<double> t(s + 1 + pb, 0);
	for (unsigned int i = 0; i < s + 1 + pb; i++) {
		if (!settled_rows[i]) {
			for (unsigned int j = 0; j < s; j++) {
				if (unsettled[j]) {
					if (B[i][j] > prec || B[i][j] < -prec)
						t[i] += B[i][j];
				}
			}
		}
	}
	row = 0;
	min_ratio = DBL_MAX;
	for (unsigned int i = 0; i < s + 1 + pb; i++) {
		if (!settled_rows[i]) {
			if (t[i] > prec) {
				if (min_ratio >= RHS[i] / t[i]) {
					if (min_ratio > RHS[i] / t[i]) {
						min_ratio = RHS[i] / t[i];
						min_idx = basis[i];
						row = i;
					}
					else {
						if (basis[i] < min_idx) {
							min_idx = basis[i];
							row = i;
						}
					}
				}
			}
		}
	}
	if (t[row] > 1 + prec || t[row] < 1 - prec) {
		for (unsigned int j = 0; j < n; j++) {
			if (A[row][j] > prec || A[row][j] < -prec)
				A[row][j] /= t[row];
		}
		for (unsigned int j = 0; j < s; j++) {
			if (unsettled[j])
				if (B[row][j] > prec || B[row][j] < -prec)
					B[row][j] /= t[row];
		}
		if (RHS[row] > prec || RHS[row] < -prec)
			RHS[row] /= t[row];
		t[row] = 1;
	}
	for (unsigned int i = 0; i < s + 1 + pb; i++) {
		if (!settled_rows[i]) {
			if (i != row) {
				if (t[i] > prec || t[i] < -prec) {
					for (unsigned int j = 0; j < n; j++)
						if (A[row][j] > prec || A[row][j] < -prec)
							A[i][j] -= A[row][j] * t[i];
					for (unsigned int j = 0; j < s; j++)
						if (unsettled[j])
							if (B[row][j] > prec || B[row][j] < -prec)
								B[i][j] -= B[row][j] * t[i];
					if (RHS[row] > prec || RHS[row] < -prec)
						RHS[i] -= t[i] * RHS[row];
				}
			}
		}
	}
	if (disp)
		cout << "t enters, " << basis[row] << " leaves the basis" << endl << endl;
	basis_ind[basis[row]] = false;
	basis_ind[n + s + 1] = true;
	basis[row] = n + s + 1;
	piv++;
	if (disp && n <= 5) {
		cout << "New table:" << endl;
		for (unsigned int i = 0; i < s + 1; i++) {
			for (unsigned int j = 0; j < n; j++) {
				if (A[i][j] >= 0)
					cout << " ";
				cout << A[i][j] << " ";
			}
			for (unsigned int j = 0; j < s; j++) {
				if (B[i][j] >= 0)
					cout << " ";
				cout << B[i][j] << " ";
			}
			cout << " " << RHS[i] << endl;
		}
		cout << "Basis:" << endl;
		for (unsigned int i = 0; i < s + 1; i++)
			cout << basis[i] << " ";
		cout << endl;
		cout << "Basis idx:" << endl;
		for (unsigned int i = 0; i < n + s + 2; i++)
			cout << basis_ind[i] << " ";
		cout << endl;
		cin.get();
	}
	solve(row, A, B, RHS, s, n, prec, min_ratio, disp, piv, basis, basis_ind, unsettled, pb, settled_rows, min_idx);
	eps += RHS[row];
	update(s, unsettled, B, A, row, n, prec, pb, RHS, settled_rows, min_ratio, basis, basis_ind, settled, disp, rank, Arref, J, a);
	if (rank == n)
		settled = s;
	if (disp) {
		cout << "New solution:" << endl;
		min_ratio = 0;
		row = s + pb + 1;
		for (unsigned int j = 0; j < n; j++) {
			if (!basis_ind[j]) {
				for (unsigned int i = 0; i < s + pb + 1; i++) {
					if (!settled_rows[i]) {
						if (A[i][j] > prec || A[i][j] < -prec) {
							if (min_ratio == 0 && A[i][j] < 1 + prec && A[i][j] > 1 - prec) {
								row = i;
								min_ratio += 1;
							}
							else {
								if (min_ratio == 1)
									row = s + pb + 1;
							}
						}
					}
				}
				if (row < s + pb + 1) {
					x[j] = RHS[row];
					row = s + pb + 1;
				}
				min_ratio = 0;
			}
		}
		for (unsigned int i = 0; i < s + pb + 1; i++) {
			if (basis[i] < n)
				x[basis[i]] = RHS[i];
		}
		for (unsigned int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "New epsilon: " << eps << endl;
	}
	if (disp && n <= 5) {
		cout << "New table:" << endl;
		for (unsigned int i = 0; i < s + pb + 1; i++) {
			if (!settled_rows[i]) {
				for (unsigned int j = 0; j < n; j++) {
					if (A[i][j] >= 0)
						cout << " ";
					cout << A[i][j] << " ";
				}
				for (unsigned int j = 0; j < s; j++) {
					if (unsettled[j]) {
						if (B[i][j] >= 0)
							cout << " ";
						cout << B[i][j] << " ";
					}
				}
				cout << " " << RHS[i] << endl;
			}
		}
		cout << "Basis:" << endl;
		for (unsigned int i = 0; i < s + 1 + pb; i++)
			if (!settled_rows[i])
				cout << basis[i] << " ";
		cout << endl;
		cout << "Basis idx:" << endl;
		for (unsigned int i = 0; i < n + s + 2; i++)
			cout << basis_ind[i] << " ";
		cout << endl;
		cin.get();
	}
	iter++;
	return;
}

void update(unsigned int &s, vector<bool> &unsettled, vector<vector<double>> &B, vector<vector<double>> &A, unsigned int &row, unsigned short int &n, double &prec, unsigned short int &pb, vector<double> &RHS, vector<bool> &settled_rows, double &min_ratio, vector<unsigned int> &basis, vector<bool> &basis_ind, unsigned int &settled, bool &disp, unsigned short int &rank, vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &a) {
	for (unsigned int j = 0; j < s; j++) {
		if (unsettled[j]) {
			if (B[row][j] > prec) {
				if (disp)
					cout << j << " is settled";
				de2bi(j, a, n);
				if (binrank(Arref, J, a, n)) {
					rank++;
					if (disp)
						cout << ", rank increased to " << rank;
					if (rank == n) {
						if (disp)
							cout << endl;
						return;
					}
					rowechform(Arref, J, a, n, rank);
				}
				unsettled[j] = false;
				settled++;
				if (disp)
					cout << endl;
			}
		}
	}
	for (unsigned short int j = 0; j < n; j++) {
		if (A[row][j] > prec) {
			for (unsigned int i = 0; i < s + 1 + pb; i++)
				A[i][j] = 0;
			A.push_back(vector<double>(n, 0));
			pb++;
			A[s + pb][j] = 1;
			RHS.push_back(0);
			B.push_back(vector<double>(s, 0));
			basis.push_back(j);
			basis_ind[j] = true;
			settled_rows.push_back(false);
			if (disp)
				cout << "x_" << j << " = 0" << endl;
		}
	}
	settled_rows[row] = true;
	min_ratio = 0;
	unsigned int j = 0;
	if (disp && n <= 5)
		cout << "Searching for elementary rows: " << endl;
	for (unsigned int i = 0; i < s + 1 + pb; i++) {
		if (disp && n <= 5)
			cout << "Row " << i << ": ";
		if (basis[i] >= n) {
			if (!settled_rows[i]) {
				while (j < n) {
					if (A[i][j] > prec || A[i][j] < -prec)
						break;
					j++;
				}
				if (j == n) {
					if (disp && n <= 5)
						cout << "A is okay... ";
					while (j < n + s) {
						if (unsettled[j - n]) {
							if (B[i][j - n] > prec || B[i][j - n] < -prec) {
								min_ratio += 1;
								row = j - n;
								if (disp && n <= 5)
									cout << "non-zero in B at col " << j - n << "... ";
							}
							if (min_ratio == 2) {
								if (disp && n <= 5)
									cout << "that's the second one!" << endl;
								break;
							}
						}
						j++;
					}
					if (j == n + s && min_ratio == 1) {
						if (disp)
							cout << row << " is settled";
						de2bi(row, a, n);
						if (binrank(Arref, J, a, n)) {
							rank++;
							if (disp)
								cout << ", rank increased to " << rank;
							if (rank == n) {
								if (disp)
									cout << endl;
								return;
							}
							rowechform(Arref, J, a, n, rank);
						}
						unsettled[row] = false;
						settled_rows[i] = true;
						settled++;
						if (disp)
							cout << endl;
					}
					min_ratio = 0;
				}
				j = 0;
			}
			else {
				if (disp && n <= 5)
					cout << "settled row" << endl;
			}
		}
		else {
			if (disp && n <= 5)
				cout << "basis row for x_" << basis[i] << endl;
		}
	}
	return;
}

void solve(unsigned int &J, vector<vector<double>> &A, vector<vector<double>> &B, vector<double> &RHS, unsigned int &s, unsigned short int &n, double &prec, double &min_ratio, bool &disp, unsigned int &piv, vector<unsigned int> &basis, vector<bool> &basis_ind, vector<bool> &unsettled, unsigned short int &pb, vector<bool> &settled_rows, unsigned int &min_idx) {
	bool opt = false;
	unsigned int col = INT_MAX;
	unsigned int row = INT_MAX;
	min_ratio = DBL_MAX;
	min_idx = INT_MAX;
	while (!opt) {
		for (unsigned short int j = 0; j < n; j++){
			if (!basis_ind[j]) {
				if (A[J][j] < -prec) {
					col = j;
					break;
				}
			}
		}
		if (col == INT_MAX) {
			for (unsigned short int j = 0; j < s; j++) {
				if (unsettled[j]) {
					if (!basis_ind[n + j]) {
						if (B[J][j] < -prec) {
							col = n + j;
							break;
						}
					}
				}
			}
		}
		if (col == INT_MAX)
			break;
		else { // find pivot row for col
			if (col < n) {
				for (unsigned short int i = 0; i < s + 1 + pb; i++) {
					if (!settled_rows[i]) {
						if (A[i][col] > prec) {
							if (RHS[i] / A[i][col] <= min_ratio) {
								if (RHS[i] / A[i][col] < min_ratio) {
									min_ratio = RHS[i] / A[i][col];
									min_idx = basis[i];
									row = i;
								}
								else {
									if (min_idx > basis[i]) {
										min_idx = basis[i];
										row = i;
									}
								}
							}
						}
					}
				}
			}
			else {
				for (unsigned short int i = 0; i < s + 1 + pb; i++) {
					if (!settled_rows[i]) {
						if (B[i][col - n] > prec) {
							if (RHS[i] / B[i][col - n] <= min_ratio) {
								if (RHS[i] / B[i][col - n] < min_ratio) {
									min_ratio = RHS[i] / B[i][col - n];
									min_idx = basis[i];
									row = i;
								}
								else {
									if (min_idx > basis[i]) {
										min_idx = basis[i];
										row = i;
									}
								}
							}
						}
					}
				}
			}
			// pivot on row and col
			if (disp) {
				if (col < n)
					cout << "pivot row " << row << " (with ratio " << RHS[row] / A[row][col] << "); column " << col << " (with t-coeff " << A[J][col] << ")" << endl;
				else
					cout << "pivot row " << row << " (with ratio " << RHS[row] / B[row][col - n] << "); column " << col << " (with t-coeff " << B[J][col - n] << ")" << endl;
				if (col < n)
					cout << "x_" << col << " enters, ";
				else
					cout << "y_" << col - n << " enters, ";
				if (basis[row] < n)
					cout << "x_" << basis[row] << " leaves the basis" << endl;
				else
					cout << "y_" << basis[row] - n << " leaves the basis" << endl;
				//cin.get();
			}
			pivot(n, col, row, A, B, prec, s, RHS, unsettled, pb, settled_rows);
			basis_ind[col] = true;
			basis_ind[basis[row]] = false;
			basis[row] = col;
			piv++;
			if (disp && n <= 5) {
				cout << "New table:" << endl;
				for (unsigned int i = 0; i < s + 1 + pb; i++) {
					if (!settled_rows[i]) {
						for (unsigned int j = 0; j < n; j++) {
							if (A[i][j] >= 0)
								cout << " ";
							cout << A[i][j] << " ";
						}
						for (unsigned int j = 0; j < s; j++) {
							if (unsettled[j]) {
								if (B[i][j] >= 0)
									cout << " ";
								cout << B[i][j] << " ";
							}
						}
						cout << " " << RHS[i] << endl;
					}
				}
				cout << "Basis:" << endl;
				for (unsigned int i = 0; i < s + 1 + pb; i++)
					if (!settled_rows[i])
						cout << basis[i] << " ";
				cout << endl;
				cout << "Basis idx:" << endl;
				for (unsigned int i = 0; i < n + s + 2; i++)
					cout << basis_ind[i] << " ";
				cout << endl;
				cin.get();
			}
			col = INT_MAX;
			row = INT_MAX;
			min_ratio = DBL_MAX;
			min_idx = INT_MAX;
		}
	}
	return;
}

void pivot(unsigned short int &n, unsigned int &col, unsigned int &row, vector<vector<double>> &A, vector<vector<double>> &B, double &prec, unsigned int &s, vector<double> &RHS, vector<bool> &unsettled, unsigned short int &pb, vector<bool> &settled_rows) {
	if (col < n) {
		if (A[row][col] > 1 + prec || A[row][col] < 1 - prec) {
			for (unsigned short int j = 0; j < n; j++)
				if (j != col)
					if (A[row][j] > prec || A[row][j] < -prec)
						A[row][j] /= A[row][col];
			for (unsigned short int j = 0; j < s; j++)
				if (unsettled[j])
					if (B[row][j] > prec || B[row][j] < -prec)
						B[row][j] /= A[row][col];
			RHS[row] /= A[row][col];
			A[row][col] = 1;
		}
		for (unsigned int i = 0; i < s + 1 + pb; i++) {
			if (!settled_rows[i]) {
				if (i != row) {
					for (unsigned short int j = 0; j < n; j++) {
						if (j != col) {
							if ((A[row][j] > prec || A[row][j] < -prec) && (A[i][col] > prec || A[i][col] < -prec))
								A[i][j] -= A[row][j] * A[i][col];
						}
					}
				}
			}
		}
		for (unsigned int i = 0; i < s + 1 + pb; i++) {
			if (!settled_rows[i]) {
				if (i != row) {
					for (unsigned short int j = 0; j < s; j++) {
						if (unsettled[j])
							if ((B[row][j] > prec || B[row][j] < -prec) && (A[i][col] > prec || A[i][col] < -prec))
								B[i][j] -= B[row][j] * A[i][col];
					}
					if ((RHS[row] > prec || RHS[row] < -prec) && (A[i][col] > prec || A[i][col] < -prec))
						RHS[i] -= RHS[row] * A[i][col];
					A[i][col] = 0;
				}
			}
		}
	}
	else {
		if (B[row][col - n] > 1 + prec || B[row][col - n] < 1 - prec) {
			for (unsigned short int j = 0; j < n; j++)
				if (A[row][j] > prec || A[row][j] < -prec)
					A[row][j] /= B[row][col - n];
			for (unsigned short int j = 0; j < s; j++)
				if (unsettled[j])
					if (j != col - n)
						if (B[row][j] > prec || B[row][j] < -prec)
							B[row][j] /= B[row][col - n];
			RHS[row] /= B[row][col - n];
			B[row][col - n] = 1;
		}
		for (unsigned int i = 0; i < s + 1 + pb; i++) {
			if (!settled_rows[i]) {
				if (i != row) {
					for (unsigned short int j = 0; j < n; j++) {
						if ((A[row][j] > prec || A[row][j] < -prec) && (B[i][col - n] > prec || B[i][col - n] < -prec))
							A[i][j] -= A[row][j] * B[i][col - n];
					}
				}
			}
		}
		for (unsigned int i = 0; i < s + 1 + pb; i++) {
			if (!settled_rows[i]) {
				if (i != row) {
					for (unsigned short int j = 0; j < s; j++) {
						if (unsettled[j]) {
							if (j != col - n) {
								if ((B[row][j] > prec || B[row][j] < -prec) && (B[i][col - n] > prec || B[i][col - n] < -prec))
									B[i][j] -= B[row][j] * B[i][col - n];
							}
						}
					}
					if ((RHS[row] > prec || RHS[row] < -prec) && (B[i][col - n] > prec || B[i][col - n] < -prec))
						RHS[i] -= RHS[row] * B[i][col - n];
					B[i][col - n] = 0;
				}
			}
		}
	}
	return;
}