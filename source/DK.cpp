/*
*    Nucleolus
*    DK.cpp
*    Purpose: finding the nucleolus of a cooperative game using
*             Derks and Kuipers (1997) - Implementing the simplex method
*             for computing the prenucleolus of transferable utility games
*             https://www.researchgate.net/publication/265080719
*
*    @author Marton Benedek
*    @version 1.0 18/12/2018
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

#include "DK.h"
#include "gen_game.h"

int mainzy() {
	unsigned short int n = 0;
	unsigned short int type = 0;
	unsigned int seed = 0;
	bool disp = false;
	bool memo = false;
	ifstream inp;
	cout << "Reading the input...";
	inp.open("input.txt");
	inp >> n >> type >> seed >> disp >> memo;
	inp.close();
	cout << "done!" << endl;
	if (seed == 0)
		seed = GetTickCount();
	srand(seed);
	unsigned int s = pow(2, n) - 2;
	vector<double> x(n, 0);
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	double prec = pow(10, -6);
	vector<double> excess(s, 0);
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
		cout << "Running DK..." << endl;
		double t1 = cpuTime();
		for (unsigned short int i = 0; i < n; i++) {
			singleton_bounds[i] = v[pow(2, i) - 1];
			impu += singleton_bounds[i];
		}
		x = singleton_bounds;
		for (unsigned short int i = 0; i < n; i++)
			x[i] += (v[s] - impu) / n;
		if (memo) {
			vector<bool> a(n, false);
			excess_init_sg_mem(excess, unsettled, a, x, v, s, n);
			DK_mem(disp, n, s, prec, singleton_bounds, unsettled, a, iter, piv, t, x, excess, t1);
		}
		else {
			vector<vector<bool>> A(s + 1, vector<bool>(n, false));
			A_mx(A, n, s);
			excess_init_sg(excess, unsettled, A, x, v, s, n);
			DK(disp, n, s, prec, singleton_bounds, unsettled, A, iter, piv, t, x, excess, t1);
		}
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
		cout << "Running DK..." << endl;
		double t1 = cpuTime();
		for (unsigned short int i = 0; i < n; i++) {
			singleton_bounds[i] = v[pow(2, i) - 1];
			impu += singleton_bounds[i];
		}
		x = singleton_bounds;
		for (unsigned short int i = 0; i < n; i++)
			x[i] += (v[s] - impu) / n;
		if (memo) {
			vector<bool> a(n, false);
			excess_init_mem(excess, unsettled, a, x, v, s, n);
			DK_mem(disp, n, s, prec, singleton_bounds, unsettled, a, iter, piv, t, x, excess, t1);
		}
		else {
			vector<vector<bool>> A(s + 1, vector<bool>(n, false));
			A_mx(A, n, s);
			excess_init(excess, unsettled, A, x, v, s, n);
			DK(disp, n, s, prec, singleton_bounds, unsettled, A, iter, piv, t, x, excess, t1);
		}
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

void DK(bool &disp, unsigned short int &n, unsigned int &s, double &prec, vector<double> &singleton_bounds, vector<bool> &unsettled, vector<vector<bool>> &A, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x, vector<double> &excess, double &t1) {
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<vector<bool>> Asettled(n, vector<bool>(n, 0));
	Asettled[0] = vector<bool>(n, true);
	if (disp) {
		cout << "Starting point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	vector<double> d(n, 0);
	double epsi = 0;
	while (rank < n)
		iteration(epsi, s, excess, prec, n, A, Arref, J, unsettled, rank, d, x, disp, Asettled, iter, piv, unsettled_p, singleton_bounds);
	t = cpuTime() - t1;
	cout << "DK finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;

}

void iteration(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<vector<bool>>&A, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, vector<double> &singleton_bounds) {
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int l = 0;
	while (true) {
		if (unsettled[l] && abs(excess[l] - epsi) < prec) {
			T_coord.push_back(l);
			T[l] = true;
			break;
		}
		l++;
	}
	if (disp)
		cout << "Initial tight coalition: " << l + 1 << endl;
	bool pivo = true;
	unsigned int t_size = 1;
	vector<vector<bool>> Atight(1, vector<bool>(n, false));
	Atight[0] = A[l];
	unsigned short int t2_size = 0;
	vector<bool> T2(n, false);
	vector<vector<bool>> Atight2(0, vector<bool>(n, false));
	vector<unsigned int> T2_coord(0, 0);
	while (pivo)
		pivot(epsi, s, excess, prec, n, A, Arref, J, unsettled, rank, d, x, disp, Asettled, pivo, t_size, Atight, T, T_coord, iter, piv, unsettled_p, singleton_bounds, Atight2, t2_size, T2, T2_coord);
}

void pivot(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<vector<bool>>&A, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, bool &pivo, unsigned int &t_size, vector<vector<bool>> &Atight, vector<bool> &T, vector<unsigned int> &T_coord, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, vector<double> &singleton_bounds, vector<vector<bool>> &Atight2, unsigned short int &t2_size, vector<bool> &T2, vector<unsigned int> &T2_coord) {
	IloEnv env;
	IloModel model(env);
	IloNumVarArray lambda(env, t_size + t2_size + rank, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	IloExpr q(env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(env);
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				p += lambda[j];
			if (i == 0)
				q += lambda[j];
			if (j < rank && Asettled[j][i] == true)
				p += lambda[j + t_size];
		}
		if (rank > t_size) {
			for (unsigned short int j = t_size; j < rank; j++) {
				if (Asettled[j][i] == true)
					p += lambda[j + t_size];
			}
		}
		for (unsigned int j = 0; j < t2_size; j++) {
			if (i == 0)
				lambda[j + t_size + rank].setLB(0);
			if (Atight2[j][i])
				p += lambda[j + t_size + rank];
		}
		IloConstraint r = (p == 0);
		model.add(r);
	}
	IloConstraint r = (q == 1);
	model.add(r);
	IloObjective OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	IloCplex sr(model);
	sr.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sr.setOut(env.getNullStream());
	if (disp)
		cout << endl << "   ---===   CHECKING DUAL FEASIBILITY   ===---   " << endl << endl;
	bool feas = sr.solve();
	if (disp)
		cout << "dual feasibility: " << feas << endl;
	if (!feas) {
		piv++;
		env.end();
		if (disp)
			cout << endl << "   ---===   SOLVING IMPROVING DIRECTION LP   ===---   " << endl << endl;
		imprdir(d, n, t_size, Atight, rank, Asettled, disp, t2_size, Atight2);
		if (disp)
			cout << endl << "   ---===   IMPROVING DIRECTION OBTAINED   ===---   " << endl << endl;
		if (disp) {
			cout << "Improving direction:" << endl;
			for (unsigned short int i = 0; i < n; i++)
				cout << d[i] << "    ";
			cout << endl;
		}
		unsigned short int l = 0;
		for (unsigned short int l = 0; l < t2_size; l++) {
			if (T2[log2(T2_coord[l] + 1)] && d[log2(T2_coord[l] + 1)] > prec) {
				if (disp)
					cout << "COALITION LEAVING THE TIGHT SET: " << T2_coord[l] + 1 << endl;
				T2[log2(T2_coord[l] + 1)] = false;
				T2_coord.erase(T2_coord.begin() + l);
				Atight2.erase(Atight2.begin() + l);
				t2_size--;
				l--;
			}
		}
		if (disp)
			cout << endl << "   ---===   COMPUTING STEP SIZE   ===---   " << endl << endl;
		step(T, unsettled, s, A, epsi, excess, d, n, x, disp, T_coord, prec, Atight, t_size, unsettled_p, singleton_bounds, Atight2, t2_size, T2, T2_coord);
	}
	else {
		vector<double> u(t_size, 0);
		for (unsigned short int i = 0; i < t_size; i++)
			u[i] = sr.getValue(lambda[i]);
		vector<double> u_impu(t2_size, 0);
		for (unsigned short int i = 0; i < t2_size; i++)
			u_impu[i] = sr.getValue(lambda[i + t_size + rank]);
		env.end();
		double m1;
		vec_min(m1, u);
		if (m1 < -prec) {
			if (disp) {
				cout << endl << "   ---===   FEASIBLE SYSTEM WITH NEGATIVE DUAL VARIABLE   ===---   " << endl << endl;
				cout << "T_coord: ";
				for (unsigned short int l = 0; l < t_size; l++)
					cout << T_coord[l] << "  ";
				cout << endl << "u: ";
				for (unsigned short int l = 0; l < t_size; l++)
					cout << u[l] << "  ";
				cout << endl;
			}
			unsigned int l_coord = 0;
			unsigned int min_coord = s + 1;
			for (unsigned int l = 0; l < t_size; l++) {
				if (u[l] < -prec) {
					if (T_coord[l] < min_coord) {
						min_coord = T_coord[l];
						l_coord = l;
					}
				}
			}
			if (disp)
				cout << "Coalition leaving the tight set: " << min_coord + 1 << endl;
			T[min_coord] = false;
			T_coord.erase(T_coord.begin() + l_coord);
			Atight.erase(Atight.begin() + l_coord);
			t_size--;
			piv++;
			if (disp)
				cout << endl << "   ---===   SOLVING IMPROVING DIRECTION LP   ===---   " << endl << endl;
			imprdir(d, n, t_size, Atight, rank, Asettled, disp, t2_size, Atight2);
			if (disp)
				cout << endl << "   ---===   IMPROVING DIRECTION OBTAINED   ===---   " << endl << endl;
			if (disp) {
				cout << "Improving direction:" << endl;
				for (unsigned short int i = 0; i < n; i++) {
					cout << d[i] << "    ";
				}
				cout << endl;
			}
			if (disp)
				cout << endl << "   ---===   COMPUTING STEP SIZE   ===---   " << endl << endl;
			step(T, unsettled, s, A, epsi, excess, d, n, x, disp, T_coord, prec, Atight, t_size, unsettled_p, singleton_bounds, Atight2, t2_size, T2, T2_coord);
		}
		else {
			pivo = false;
			iter++;
			if (disp)
				cout << endl << "   ---===   DUAL FEASIBLE POINT FOUND!   ===---   " << endl << endl;
			for (unsigned int i = 0; i < t_size; i++) {
				if (u[i] > prec) {
					if (binrank(Arref, J, Atight[i], n)) {
						rank++;
						if (disp)
							cout << "Rank increased to " << rank << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i] << ") getting settled." << endl;
						if (rank == n) {
							if (disp)
								cout << "Rank condition satisfied!" << endl;
							return;
						}
						Asettled[rank - 1] = Atight[i];
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						rowechform(Arref, J, Atight[i], n, rank);
					}
					else {
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
					}
				}
			}
			for (unsigned int i = 0; i < t2_size; i++) {
				if (u_impu[i] > prec) {
					if (binrank(Arref, J, Atight2[i], n)) {
						rank++;
						if (disp)
							cout << "Rank increased to " << rank << " with " << T2_coord[i] + 1 << " (and " << s - T2_coord[i] << ") getting settled." << endl;
						if (rank == n) {
							if (disp)
								cout << "Rank condition satisfied!" << endl;
							return;
						}
						Asettled[rank - 1] = Atight2[i];
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						rowechform(Arref, J, Atight2[i], n, rank);
					}
					else {
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						if (disp)
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
					}
				}
			}
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					if (!(binrank(Arref, J, A[i], n))) {
						unsettled[i] = false;
						unsettled[s - 1 - i] = false;
					}
				}
			}
			if (disp)
				cout << "Rank increased to: " << rank << endl;
			for (unsigned short int i = 0; i < n; i++)
				if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
					unsettled_p[i] = false;
		}
	}
}

void step(vector<bool> &T, vector<bool> &unsettled, unsigned int &s, vector<vector<bool>> &A, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, bool &disp, vector<unsigned int> &T_coord, double &prec, vector<vector<bool>> &Atight, unsigned int &t_size, vector<bool> &unsettled_p, vector<double> &singleton_bounds, vector<vector<bool>> &Atight2, unsigned short int &t2_size, vector<bool> &T2, vector<unsigned int> &T2_coord) {
	double alpha1 = DBL_MAX;
	double alpha2 = DBL_MAX;
	double alpha = DBL_MAX;
	double Ad;
	unsigned int coord1;
	unsigned int coord2;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha1) {
				alpha1 = (epsi - excess[i]) / (Ad - 1);
				coord1 = i;
				if (alpha1 <= 0)
					break;
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] && !T2[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha2) {
				alpha2 = (singleton_bounds[i] - x[i]) / d[i];
				coord2 = i;
				if (alpha2 <= 0)
					break;
			}
		}
	}
	if (alpha1 <= alpha2) {
		T[coord1] = true;
		T_coord.push_back(coord1);
		alpha = alpha1;
		Atight.push_back(A[T_coord[t_size]]);
		t_size++;
		epsi += alpha;
		if (disp) {
			cout << "Coalition entering the tight set: " << coord1 + 1 << endl;
			cout << "Step size: " << alpha << endl;
			cout << "Epsilon: " << epsi << endl;
		}
	}
	else {
		T2[coord2] = true;
		Atight2.push_back(vector<bool>(n, false));
		Atight2[t2_size][coord2] = true;
		T2_coord.push_back(pow(2, coord2) - 1);
		t2_size++;
		alpha = alpha2;
		epsi += alpha;
		if (disp) {
			cout << "Coalition entering the tight set: " << pow(2, coord2) << endl;
			cout << "Step size: " << alpha << endl;
			cout << "Epsilon: " << epsi << endl;
		}
	}
	if (disp)
		cout << endl << "   ---===   STEP SIZE OBTAINED   ===---   " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New x point: " << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			excess[i] += alpha*Ad;
		}
	}
}

void stepDK(vector<bool> &T, vector<bool> &unsettled, unsigned int &s, vector<vector<bool>> &A, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, bool &disp, vector<unsigned int> &T_coord, double &prec, vector<vector<bool>> &Atight, unsigned int &t_size, vector<bool> &unsettled_p, vector<double> &singleton_bounds, vector<vector<bool>> &Atight2, unsigned short int &t2_size, vector<bool> &T2, vector<unsigned int> &T2_coord) {
	double alpha1 = DBL_MAX;
	double alpha2 = DBL_MAX;
	double alpha = DBL_MAX;
	double Ad;
	unsigned int coord1;
	unsigned int coord2;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha1) {
				alpha1 = (epsi - excess[i]) / (Ad - 1);
				coord1 = i;
				if (alpha1 <= 0)
					break;
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] && !T2[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha2) {
				alpha2 = (singleton_bounds[i] - x[i]) / d[i];
				coord2 = i;
				if (alpha2 <= 0)
					break;
			}
		}
	}
	if (alpha1 <= alpha2) {
		T[coord1] = true;
		T_coord.push_back(coord1);
		alpha = alpha1;
		Atight.push_back(A[T_coord[t_size]]);
		t_size++;
		epsi += alpha;
		if (disp) {
			cout << "Coalition entering the tight set: " << coord1 + 1 << endl;
			cout << "Step size: " << alpha << endl;
			cout << "Epsilon: " << epsi << endl;
		}
	}
	else {
		T2[coord2] = true;
		Atight2.push_back(vector<bool>(n, false));
		Atight2[t2_size][coord2] = true;
		T2_coord.push_back(pow(2, coord2) - 1);
		t2_size++;
		alpha = alpha2;
		epsi += alpha;
		if (disp) {
			cout << "Coalition entering the tight set: " << pow(2, coord2) << endl;
			cout << "Step size: " << alpha << endl;
			cout << "Epsilon: " << epsi << endl;
		}
	}
	if (disp)
		cout << endl << "   ---===   STEP SIZE OBTAINED   ===---   " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New x point: " << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			excess[i] += alpha*Ad;
		}
	}
}

void DK_mem(bool &disp, unsigned short int &n, unsigned int &s, double &prec, vector<double> &singleton_bounds, vector<bool> &unsettled, vector<bool> &a, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x, vector<double> &excess, double &t1) {
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<vector<bool>> Asettled(n, vector<bool>(n, 0));
	Asettled[0] = vector<bool>(n, true);
	if (disp) {
		cout << "Starting point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	vector<double> d(n, 0);
	double epsi = 0;
	while (rank < n)
		iteration_mem(epsi, s, excess, prec, n, a, Arref, J, unsettled, rank, d, x, disp, Asettled, iter, piv, unsettled_p, singleton_bounds);
	t = cpuTime() - t1;
	cout << "DK finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_mem(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<bool>&a, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, vector<double> &singleton_bounds) {
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int l = 0;
	while (true) {
		if (unsettled[l] && abs(excess[l] - epsi) < prec) {
			T_coord.push_back(l);
			T[l] = true;
			break;
		}
		l++;
	}
	if (disp)
		cout << "Initial tight coalition: " << l + 1 << endl;
	bool pivo = true;
	unsigned int t_size = 1;
	vector<vector<bool>> Atight(1, vector<bool>(n, false));
	de2bi(l, Atight[0], n);
	unsigned short int t2_size = 0;
	vector<bool> T2(n, false);
	vector<vector<bool>> Atight2(0, vector<bool>(n, false));
	vector<unsigned int> T2_coord(0, 0);
	while (pivo)
		pivot_mem(epsi, s, excess, prec, n, a, Arref, J, unsettled, rank, d, x, disp, Asettled, pivo, t_size, Atight, T, T_coord, iter, piv, unsettled_p, singleton_bounds, Atight2, t2_size, T2, T2_coord);
}

void pivot_mem(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<bool>&a, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, bool &pivo, unsigned int &t_size, vector<vector<bool>> &Atight, vector<bool> &T, vector<unsigned int> &T_coord, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, vector<double> &singleton_bounds, vector<vector<bool>> &Atight2, unsigned short int &t2_size, vector<bool> &T2, vector<unsigned int> &T2_coord) {
	IloEnv env;
	IloModel model(env);
	IloNumVarArray lambda(env, t_size + t2_size + rank, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	IloExpr q(env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(env);
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				p += lambda[j];
			if (i == 0)
				q += lambda[j];
			if (j < rank && Asettled[j][i] == true)
				p += lambda[j + t_size];
		}
		if (rank > t_size) {
			for (unsigned short int j = t_size; j < rank; j++) {
				if (Asettled[j][i] == true)
					p += lambda[j + t_size];
			}
		}
		for (unsigned int j = 0; j < t2_size; j++) {
			if (i == 0)
				lambda[j + t_size + rank].setLB(0);
			if (Atight2[j][i])
				p += lambda[j + t_size + rank];
		}
		IloConstraint r = (p == 0);
		model.add(r);
	}
	IloConstraint r = (q == 1);
	model.add(r);
	IloObjective OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	IloCplex sr(model);
	sr.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sr.setOut(env.getNullStream());
	if (disp)
		cout << endl << "   ---===   CHECKING DUAL FEASIBILITY   ===---   " << endl << endl;
	bool feas = sr.solve();
	if (disp)
		cout << "dual feasibility: " << feas << endl;
	if (!feas) {
		piv++;
		env.end();
		if (disp)
			cout << endl << "   ---===   SOLVING IMPROVING DIRECTION LP   ===---   " << endl << endl;
		imprdir(d, n, t_size, Atight, rank, Asettled, disp, t2_size, Atight2);
		if (disp)
			cout << endl << "   ---===   IMPROVING DIRECTION OBTAINED   ===---   " << endl << endl;
		if (disp) {
			cout << "Improving direction:" << endl;
			for (unsigned short int i = 0; i < n; i++)
				cout << d[i] << "    ";
			cout << endl;
		}
		unsigned short int l = 0;
		for (unsigned short int l = 0; l < t2_size; l++) {
			if (T2[log2(T2_coord[l] + 1)] && d[log2(T2_coord[l] + 1)] > prec) {
				if (disp)
					cout << "COALITION LEAVING THE TIGHT SET: " << T2_coord[l] + 1 << endl;
				T2[log2(T2_coord[l] + 1)] = false;
				T2_coord.erase(T2_coord.begin() + l);
				Atight2.erase(Atight2.begin() + l);
				t2_size--;
				l--;
			}
		}
		if (disp)
			cout << endl << "   ---===   COMPUTING STEP SIZE   ===---   " << endl << endl;
		step_mem(T, unsettled, s, a, epsi, excess, d, n, x, disp, T_coord, prec, Atight, t_size, unsettled_p, singleton_bounds, Atight2, t2_size, T2, T2_coord);
	}
	else {
		vector<double> u(t_size, 0);
		for (unsigned short int i = 0; i < t_size; i++)
			u[i] = sr.getValue(lambda[i]);
		vector<double> u_impu(t2_size, 0);
		for (unsigned short int i = 0; i < t2_size; i++)
			u_impu[i] = sr.getValue(lambda[i + t_size + rank]);
		env.end();
		double m1;
		vec_min(m1, u);
		if (m1 < -prec) {
			if (disp)
				cout << endl << "   ---===   FEASIBLE SYSTEM WITH NEGATIVE DUAL VARIABLE   ===---   " << endl << endl;
			if (disp) {
				cout << "T_coord: ";
				for (unsigned short int l = 0; l < t_size; l++)
					cout << T_coord[l] << "  ";
				cout << endl << "u: ";
				for (unsigned short int l = 0; l < t_size; l++)
					cout << u[l] << "  ";
				cout << endl;
			}
			unsigned int l_coord = 0;
			unsigned int min_coord = s + 1;
			for (unsigned int l = 0; l < t_size; l++) {
				if (u[l] < -prec) {
					if (T_coord[l] < min_coord) {
						min_coord = T_coord[l];
						l_coord = l;
					}
				}
			}
			if (disp)
				cout << "Coalition leaving the tight set: " << min_coord + 1 << endl;
			T[min_coord] = false;
			T_coord.erase(T_coord.begin() + l_coord);
			Atight.erase(Atight.begin() + l_coord);
			t_size--;
			piv++;
			if (disp)
				cout << endl << "   ---===   SOLVING IMPROVING DIRECTION LP   ===---   " << endl << endl;
			imprdir(d, n, t_size, Atight, rank, Asettled, disp, t2_size, Atight2);
			if (disp)
				cout << endl << "   ---===   IMPROVING DIRECTION OBTAINED   ===---   " << endl << endl;
			if (disp) {
				cout << "Improving direction:" << endl;
				for (unsigned short int i = 0; i < n; i++) {
					cout << d[i] << "    ";
				}
				cout << endl;
			}
			if (disp)
				cout << endl << "   ---===   COMPUTING STEP SIZE   ===---   " << endl << endl;
			step_mem(T, unsettled, s, a, epsi, excess, d, n, x, disp, T_coord, prec, Atight, t_size, unsettled_p, singleton_bounds, Atight2, t2_size, T2, T2_coord);
		}
		else {
			pivo = false;
			iter++;
			if (disp)
				cout << endl << "   ---===   DUAL FEASIBLE POINT FOUND!   ===---   " << endl << endl;
			for (unsigned int i = 0; i < t_size; i++) {
				if (u[i] > prec) {
					if (binrank(Arref, J, Atight[i], n)) {
						rank++;
						if (disp)
							cout << "Rank increased to " << rank << " with " << T_coord[i] + 1 << "(and " << s - T_coord[i]  << ") getting settled." << endl;
						if (rank == n) {
							if (disp)
								cout << "Rank condition satisfied!" << endl;
							return;
						}
						Asettled[rank - 1] = Atight[i];
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						rowechform(Arref, J, Atight[i], n, rank);
					}
					else {
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
					}
				}
			}
			for (unsigned int i = 0; i < t2_size; i++) {
				if (u_impu[i] > prec) {
					if (binrank(Arref, J, Atight2[i], n)) {
						rank++;
						if (disp)
							cout << "Rank increased to " << rank << " with " << T2_coord[i] + 1 << "(and " << s - T2_coord[i] << ") getting settled." << endl;
						if (rank == n) {
							if (disp)
								cout << "Rank condition satisfied!" << endl;
							return;
						}
						Asettled[rank - 1] = Atight2[i];
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						rowechform(Arref, J, Atight2[i], n, rank);
					}
					else {
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						if (disp)
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
					}
				}
			}
			if (disp)
				cout << "Rank increased to: " << rank << endl;
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					de2bi(i, a, n);
					if (!(binrank(Arref, J, a, n))) {
						unsettled[i] = false;
						unsettled[s - 1 - i] = false;
					}
				}
			}
			for (unsigned short int i = 0; i < n; i++)
				if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
					unsettled_p[i] = false;
		}
	}
}

void step_mem(vector<bool> &T, vector<bool> &unsettled, unsigned int &s, vector<bool> &a, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, bool &disp, vector<unsigned int> &T_coord, double &prec, vector<vector<bool>> &Atight, unsigned int &t_size, vector<bool> &unsettled_p, vector<double> &singleton_bounds, vector<vector<bool>> &Atight2, unsigned short int &t2_size, vector<bool> &T2, vector<unsigned int> &T2_coord) {
	double alpha1 = DBL_MAX;
	double alpha2 = DBL_MAX;
	double alpha = DBL_MAX;
	double Ad;
	unsigned int coord1;
	unsigned int coord2;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha1) {
				alpha1 = (epsi - excess[i]) / (Ad - 1);
				coord1 = i;
				if (alpha1 <= 0)
					break;
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] && !T2[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha2) {
				alpha2 = (singleton_bounds[i] - x[i]) / d[i];
				coord2 = i;
				if (alpha2 <= 0)
					break;
			}
		}
	}
	if (alpha1 <= alpha2) {
		T[coord1] = true;
		T_coord.push_back(coord1);
		alpha = alpha1;
		de2bi(T_coord[t_size], a, n);
		Atight.push_back(a);
		t_size++;
		epsi += alpha;
		if (disp) {
			cout << "Coalition entering the tight set: " << coord1 + 1 << endl;
			cout << "Step size: " << alpha << endl;
			cout << "Epsilon: " << epsi << endl;
		}
	}
	else {
		T2[coord2] = true;
		Atight2.push_back(vector<bool>(n, false));
		Atight2[t2_size][coord2] = true;
		T2_coord.push_back(pow(2, coord2) - 1);
		t2_size++;
		alpha = alpha2;
		epsi += alpha;
		if (disp) {
			cout << "Coalition entering the tight set: " << pow(2, coord2) << endl;
			cout << "Step size: " << alpha << endl;
			cout << "Epsilon: " << epsi << endl;
		}
	}
	if (disp)
		cout << endl << "   ---===   STEP SIZE OBTAINED   ===---   " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New x point: " << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			de2bi(i, a, n);
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					Ad += d[j];
			}
			excess[i] += alpha*Ad;
		}
	}
}

void imprdir(vector<double>&d, unsigned short int &n, unsigned int &t_size, vector<vector<bool>> &Atight, unsigned short int &rank, vector<vector<bool>>&Asettled, bool &disp, unsigned short int &t2_size, vector<vector<bool>> &Atight2) {
	IloEnv dir_env;
	IloModel dir_model(dir_env);
	IloNumVarArray D(dir_env, n, -IloInfinity, IloInfinity);
	IloExpr dir_obj(dir_env);
	for (unsigned int i = 0; i < t_size; i++) {
		IloExpr eq(dir_env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Atight[i][j])
				eq += D[j];
		}
		IloConstraint r;
		r = (eq == 1);
		dir_model.add(r);
	}
	for (unsigned int i = 0; i < t2_size; i++) {
		for (unsigned short int j = 0; j < n; j++) {
			if (Atight2[i][j]) {
				IloConstraint q = (D[j] >= 0);
				dir_model.add(q);
			}
		}
	}
	for (unsigned short int i = 0; i < rank; i++) {
		IloExpr eq(dir_env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Asettled[i][j])
				eq += D[j];
		}
		IloConstraint r = (eq == 0);
		dir_model.add(r);
	}
	dir_model.add(IloMinimize(dir_env, dir_obj));
	IloCplex sr(dir_model);
	sr.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sr.setOut(dir_env.getNullStream());
	sr.solve();
	for (unsigned short int i = 0; i < n; i++)
		d[i] = sr.getValue(D[i]);
	dir_env.end();
}

void vec_min(double&m, vector<double>&x) {
	// finds the minimum value m of vector (double) x
	m = x[0];
	for (unsigned int i = 1; i != x.size(); i++) {
		if (x[i] < m)
			m = x[i];
	}
}