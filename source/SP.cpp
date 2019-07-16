/*
*    Nucleolus
*    SP.cpp
*    Purpose: finding the nucleolus of a cooperative game using the primal
*             sequence of Solymosi (1993) - On computing the nucleolus of
*             cooperative games
*             https://www.researchgate.net/publication/318017147
*
*    @author Marton Benedek
*    @version 1.1 16/07/2019
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

#include "SP.h"
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
	unsigned short int iter = 0;
	unsigned int piv = 0;
	unsigned int sr = 0;
	double t = 0;
	if (type == 3 || type == 5) {
		cout << "Generating game...";
		vector<bool> v(s + 1, 0);
		if (type == 3)
			type3(v, s, n);
		else if (type == 5)
			type5(v, s, n);
		cout << "done!" << endl;
		cout << "Running SP..." << endl;
		if (memo)
			SP_sg_mem(disp, n, v, iter, piv, sr, t, x, s, nlsu);
		else
			SP_sg(disp, n, v, iter, piv, sr, t, x, s, nlsu);
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
		cout << "Running SP..." << endl;
		if (memo)
			SP_mem(disp, n, v, iter, piv, sr, t, x, s, nlsu);
		else
			SP(disp, n, v, iter, piv, sr, t, x, s, nlsu);
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

void SP(bool &disp, unsigned short int &n, vector<double> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<vector<bool>> A(s + 1, vector<bool>(n, false));
	A_mx(A, n, s);
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	for (unsigned short int i = 0; i < n; i++){
		singleton_bounds[i] = v[pow(2, i) - 1];
		impu += singleton_bounds[i]
	}
	for (unsigned short int i = 0; i < n; i++)
		x[i] += (v[s] - impu) / n;
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	unsigned int sett = 1;
	vector<vector<bool>> Asettled(1, vector<bool>(n, true));
	vector<double> settled_values(1, v[s]);
	double epsi = 0;
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		IloExpr ineq(env);
		for (unsigned short int j = 0; j < n; j++) {
			if (A[i][j])
				ineq += X[j];
			if (i == 0) {
				eq += X[j];
				IloExpr Impu(env);
				Impu = X[j];
				IloRange impu_constr = (Impu >= singleton_bounds[j]);
				model.add(impu_constr);
			}
		}
		ineq += X[n];
		IloRange unsett_ineq = (ineq >= v[i]);
		model.add(unsett_ineq);
	}
	IloConstraint r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val,NULL,X,NULL,NULL,NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS;
	while (rank < n)
		iteration(unsettled, s, xS, n, A, x, v, epsi, prec, Arref, J, rank, disp, Asettled, sett, settled_values, iter, piv, sr, unsettled_p, singleton_bounds, nlsu);
	t = cpuTime() - t1;
	cout << "SP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	cout << "Subroutine solves needed: " << sr << endl;
}

void iteration(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds, bool &nlsu) {
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned short int t2_size = 0;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			xS = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					xS += x[j];
			}
			if (abs(v[i] - xS - epsi) < prec) {
				T[i] = true;
				T_coord.push_back(i);
				t_size++;
			}
		}
	}
	for (unsigned int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (abs(x[i] - singleton_bounds[i]) < prec) {
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
				t2_size++;
			}
		}
	}
	if (disp) {
		cout << "Tight coalitions:" << endl;
		for (unsigned short int i = 0; i < t_size; i++)
			cout << T_coord[i] + 1 << endl;
	}
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	unsigned int l = 0;
	for (unsigned int i = 0; i < t_size; i++)
		Atight[i] = A[T_coord[i]];
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		Atight2[i] = A[T2_coord[i]];
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sett, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord);
	if (disp) {
		cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
		cout << "MIN TIGHT SET FOUND!" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << endl;
	}
	if (rank == n)
		return;
	if (disp)
		cout << "Rank increased to: " << rank << endl;
	if (!nlsu) {
		for (unsigned int i = 0; i < s; i++) {
			if (unsettled[i]) {
				if (!(binrank(Arref, J, A[i], n))) {
					unsettled[i] = false;
					unsettled[s - 1 - i] = false;
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			unsettled_p[i] = false;
	}
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			IloExpr ineq(env);
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					ineq += X[j];
			}
			ineq += X[n];
			IloRange unsett_ineq = (ineq >= v[i]);
			model.add(unsett_ineq);
		}
	}
	for (unsigned int i = 1; i < sett; i++) {
		IloExpr settled(env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Asettled[i][j])
				settled += X[j];
			if (i == 1) {
				eq += X[j];
				if (unsettled_p[j]) {
					IloExpr Impu(env);
					Impu = X[j];
					IloRange impu_constr = (Impu >= singleton_bounds[j]);
					model.add(impu_constr);
				}
			}
		}
		IloRange sett_ineq = (settled >= settled_values[i]);
		model.add(sett_ineq);
	}
	IloRange r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "New solution point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void SP_sg(bool &disp, unsigned short int &n, vector<bool> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<vector<bool>> A(s + 1, vector<bool>(n, false));
	A_mx(A, n, s);
	vector<bool> singleton_bounds(n, 0);
	double impu = 0;
	for (unsigned short int i = 0; i < n; i++){
		singleton_bounds[i] = v[pow(2, i) - 1];
		impu += singleton_bounds[i];
	}
	for (unsigned short int i = 0; i < n; i++)
		x[i] += (v[s] - impu) / n;
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	unsigned int sett = 1;
	vector<vector<bool>> Asettled(1, vector<bool>(n, true));
	vector<double> settled_values(1, v[s]);
	double epsi = 0;
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		IloExpr ineq(env);
		for (unsigned short int j = 0; j < n; j++) {
			if (A[i][j])
				ineq += X[j];
			if (i == 0) {
				eq += X[j];
				IloExpr Impu(env);
				Impu = X[j];
				IloRange impu_constr = (Impu >= singleton_bounds[j]);
				model.add(impu_constr);
			}
		}
		ineq += X[n];
		IloRange unsett_ineq = (ineq >= v[i]);
		model.add(unsett_ineq);
	}
	IloConstraint r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS;
	while (rank < n)
		iteration_sg(unsettled, s, xS, n, A, x, v, epsi, prec, Arref, J, rank, disp, Asettled, sett, settled_values, iter, piv, sr, unsettled_p, singleton_bounds, nlsu);
	t = cpuTime() - t1;
	cout << "SP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	cout << "Subroutine solves needed: " << sr << endl;
}

void iteration_sg(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<bool> &singleton_bounds, bool &nlsu) {
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned short int t2_size = 0;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			xS = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					xS += x[j];
			}
			if (abs(v[i] - xS - epsi) < prec) {
				T[i] = true;
				T_coord.push_back(i);
				t_size++;
			}
		}
	}
	for (unsigned int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (abs(x[i] - singleton_bounds[i]) < prec) {
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
				t2_size++;
			}
		}
	}
	if (disp) {
		cout << "Tight coalitions:" << endl;
		for (unsigned short int i = 0; i < t_size; i++)
			cout << T_coord[i] + 1 << endl;
	}
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	unsigned int l = 0;
	for (unsigned int i = 0; i < t_size; i++)
		Atight[i] = A[T_coord[i]];
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		Atight2[i] = A[T2_coord[i]];
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	subroutine_sg(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sett, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord);
	if (disp) {
		cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
		cout << "MIN TIGHT SET FOUND!" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << endl;
	}
	if (rank == n)
		return;
	if (disp)
		cout << "Rank increased to: " << rank << endl;
	if (!nlsu) {
		for (unsigned int i = 0; i < s; i++) {
			if (unsettled[i]) {
				if (!(binrank(Arref, J, A[i], n))) {
					unsettled[i] = false;
					unsettled[s - 1 - i] = false;
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			unsettled_p[i] = false;
	}
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			IloExpr ineq(env);
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					ineq += X[j];
			}
			ineq += X[n];
			IloRange unsett_ineq = (ineq >= v[i]);
			model.add(unsett_ineq);
		}
	}
	for (unsigned int i = 1; i < sett; i++) {
		IloExpr settled(env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Asettled[i][j])
				settled += X[j];
			if (i == 1) {
				eq += X[j];
				if (unsettled_p[j]) {
					IloExpr Impu(env);
					Impu = X[j];
					IloRange impu_constr = (Impu >= singleton_bounds[j]);
					model.add(impu_constr);
				}
			}
		}
		IloRange sett_ineq = (settled >= settled_values[i]);
		model.add(sett_ineq);
	}
	IloRange r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "New solution point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void SP_mem(bool &disp, unsigned short int &n, vector<double> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<bool> a(n, false);
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	for (unsigned short int i = 0; i < n; i++){
		singleton_bounds[i] = v[pow(2, i) - 1];
		impu += singleton_bounds[i];
	}
	for (unsigned short int i = 0; i < n; i++)
		x[i] += (v[s] - impu) / n;
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	unsigned int sett = 1;
	vector<vector<bool>> Asettled(1, vector<bool>(n, true));
	vector<double> settled_values(1, v[s]);
	double epsi = 0;
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		IloExpr ineq(env);
		de2bi(i, a, n);
		for (unsigned short int j = 0; j < n; j++) {
			if (a[j])
				ineq += X[j];
			if (i == 0) {
				eq += X[j];
				IloExpr Impu(env);
				Impu = X[j];
				IloRange impu_constr = (Impu >= singleton_bounds[j]);
				model.add(impu_constr);
			}
		}
		ineq += X[n];
		IloRange unsett_ineq = (ineq >= v[i]);
		model.add(unsett_ineq);
	}
	IloConstraint r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS;
	while (rank < n)
		iteration_mem(unsettled, s, xS, n, a, x, v, epsi, prec, Arref, J, rank, disp, Asettled, sett, settled_values, iter, piv, sr, unsettled_p, singleton_bounds, nlsu);
	t = cpuTime() - t1;
	cout << "SP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	cout << "Subroutine solves needed: " << sr << endl;
}

void iteration_mem(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds, bool &nlsu) {
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned short int t2_size = 0;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			xS = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					xS += x[j];
			}
			if (abs(v[i] - xS - epsi) < prec) {
				T[i] = true;
				T_coord.push_back(i);
				t_size++;
			}
		}
	}
	for (unsigned int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (abs(x[i] - singleton_bounds[i]) < prec) {
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
				t2_size++;
			}
		}
	}
	if (disp) {
		cout << "Tight coalitions:" << endl;
		for (unsigned short int i = 0; i < t_size; i++)
			cout << T_coord[i] + 1 << endl;
	}
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	unsigned int l = 0;
	for (unsigned int i = 0; i < t_size; i++)
		de2bi(T_coord[i], Atight[i], n);
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		de2bi(T2_coord[i], Atight2[i], n);
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sett, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord);
	if (disp) {
		cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
		cout << "MIN TIGHT SET FOUND!" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << endl;
	}
	if (rank == n)
		return;
	if (disp)
		cout << "Rank increased to: " << rank << endl;
	if (!nlsu) {
		for (unsigned int i = 0; i < s; i++) {
			if (unsettled[i]) {
				de2bi(i, a, n);
				if (!(binrank(Arref, J, a, n))) {
					unsettled[i] = false;
					unsettled[s - 1 - i] = false;
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			unsettled_p[i] = false;
	}
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			IloExpr ineq(env);
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					ineq += X[j];
			}
			ineq += X[n];
			IloRange unsett_ineq = (ineq >= v[i]);
			model.add(unsett_ineq);
		}
	}
	for (unsigned int i = 1; i < sett; i++) {
		IloExpr settled(env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Asettled[i][j])
				settled += X[j];
			if (i == 1) {
				eq += X[j];
				if (unsettled_p[j]) {
					IloExpr Impu(env);
					Impu = X[j];
					IloRange impu_constr = (Impu >= singleton_bounds[j]);
					model.add(impu_constr);
				}
			}
		}
		IloRange sett_ineq = (settled >= settled_values[i]);
		model.add(sett_ineq);
	}
	IloRange r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++) {
		x[i] = lp.getValue(X[i]);
	}
	if (disp) {
		cout << "New solution point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void SP_sg_mem(bool &disp, unsigned short int &n, vector<bool> &v, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, unsigned int &s, bool &nlsu) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<bool> a(n, false);
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	for (unsigned short int i = 0; i < n; i++){
		singleton_bounds[i] = v[pow(2, i) - 1];
		impu += singleton_bounds[i];
	}
	for (unsigned short int i = 0; i < n; i++)
		x[i] += (v[s] - impu) / n;
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	unsigned int sett = 1;
	vector<vector<bool>> Asettled(1, vector<bool>(n, true));
	vector<double> settled_values(1, v[s]);
	double epsi = 0;
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		IloExpr ineq(env);
		de2bi(i, a, n);
		for (unsigned short int j = 0; j < n; j++) {
			if (a[j])
				ineq += X[j];
			if (i == 0) {
				eq += X[j];
				IloExpr Impu(env);
				Impu = X[j];
				IloRange impu_constr = (Impu >= singleton_bounds[j]);
				model.add(impu_constr);
			}
		}
		ineq += X[n];
		IloRange unsett_ineq = (ineq >= v[i]);
		model.add(unsett_ineq);
	}
	IloConstraint r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS;
	while (rank < n)
		iteration_sg_mem(unsettled, s, xS, n, a, x, v, epsi, prec, Arref, J, rank, disp, Asettled, sett, settled_values, iter, piv, sr, unsettled_p, singleton_bounds, nlsu);
	t = cpuTime() - t1;
	cout << "SP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	cout << "Subroutine solves needed: " << sr << endl;
}

void iteration_sg_mem(vector<bool> &unsettled, unsigned int &s, double &xS, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, vector<double> &settled_values, unsigned short int &iter, unsigned int &piv, unsigned int &sr, vector<bool> &unsettled_p, vector<double> &singleton_bounds, bool &nlsu) {
	IloEnv env;
	IloModel model(env);
	IloNumVarArray X(env, n + 1, -IloInfinity, IloInfinity);
	IloExpr obj(env);
	obj = X[n];
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned short int t2_size = 0;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			xS = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					xS += x[j];
			}
			if (abs(v[i] - xS - epsi) < prec) {
				T[i] = true;
				T_coord.push_back(i);
				t_size++;
			}
		}
	}
	for (unsigned int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (abs(x[i] - singleton_bounds[i]) < prec) {
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
				t2_size++;
			}
		}
	}
	if (disp) {
		cout << "Tight coalitions:" << endl;
		for (unsigned short int i = 0; i < t_size; i++)
			cout << T_coord[i] + 1 << endl;
	}
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	unsigned int l = 0;
	for (unsigned int i = 0; i < t_size; i++)
		de2bi(T_coord[i], Atight[i], n);
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		de2bi(T2_coord[i], Atight2[i], n);
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	subroutine_sg(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sett, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord);
	if (disp) {
		cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
		cout << "MIN TIGHT SET FOUND!" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << endl;
	}
	if (rank == n)
		return;
	if (disp)
		cout << "Rank increased to: " << rank << endl;
	if (!nlsu) {
		for (unsigned int i = 0; i < s; i++) {
			if (unsettled[i]) {
				de2bi(i, a, n);
				if (!(binrank(Arref, J, a, n))) {
					unsettled[i] = false;
					unsettled[s - 1 - i] = false;
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			unsettled_p[i] = false;
	}
	IloExpr eq(env);
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			IloExpr ineq(env);
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					ineq += X[j];
			}
			ineq += X[n];
			IloRange unsett_ineq = (ineq >= v[i]);
			model.add(unsett_ineq);
		}
	}
	for (unsigned int i = 1; i < sett; i++) {
		IloExpr settled(env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Asettled[i][j])
				settled += X[j];
			if (i == 1) {
				eq += X[j];
				if (unsettled_p[j]) {
					IloExpr Impu(env);
					Impu = X[j];
					IloRange impu_constr = (Impu >= singleton_bounds[j]);
					model.add(impu_constr);
				}
			}
		}
		IloRange sett_ineq = (settled >= settled_values[i]);
		model.add(sett_ineq);
	}
	IloRange r = (eq == v[s]);
	model.add(r);
	model.add(IloMinimize(env, obj));
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	IloNumArray x_val(env, n + 1);
	for (unsigned short int i = 0; i < n; i++)
		x_val[i] = x[i];
	x_val[n] = epsi;
	lp.setStart(x_val, NULL, X, NULL, NULL, NULL);
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getValue(X[i]);
	if (disp) {
		cout << "New solution point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	epsi = lp.getValue(X[n]);
	env.end();
	if (disp) {
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void subroutine(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, unsigned int &sr, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <double> &v, vector<unsigned int> &T2_coord) {
	unsigned int sumt = 0;
	vector<bool> t(t_size, false);
	unsigned short int sumt2 = 0;
	vector<bool> t2(t2_size, false);
	IloEnv sr_env;
	IloModel sr_model(sr_env);
	IloNumVarArray lambda(sr_env, t_size + t2_size + sett, 0, IloInfinity);
	lambda[t_size + t2_size].setLB(-IloInfinity);
	IloExpr sr_obj(sr_env);
	for (unsigned int i = 0; i < t_size + t2_size + sett; i++)
		sr_model.add(lambda[i]);
	IloExpr q(sr_env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(sr_env);
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				p += lambda[j];
			if (i == 0)
				q += lambda[j];
			if (j < sett && Asettled[j][i] == true)
				p += lambda[j + t_size + t2_size];
		}
		for (unsigned int j = 0; j < t2_size; j++) {
			if (Atight2[j][i] == true)
				p += lambda[j + t_size];
		}
		if (sett > t_size) {
			for (unsigned short int j = t_size; j < sett; j++) {
				if (Asettled[j][i] == true)
					p += lambda[j + t_size + t2_size];
			}
		}
		IloConstraint r = (p == 0);
		sr_model.add(r);
	}
	IloConstraint r = (q == 1);
	sr_model.add(r);
	IloObjective OBJ = IloMaximize(sr_env, sr_obj);
	sr_model.add(OBJ);
	IloCplex SR(sr_model);
	SR.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		SR.setOut(sr_env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING SUBROUTINE LP   ===---   " << endl << endl;
	bool feas = SR.solve();
	sr++;
	unsigned int i;
	while (feas) {
		subr_upd(Arref, J, i, q, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, t_size, t2_size, SR, lambda, rank, disp, sett, Asettled, settled_values, unsettled, T_coord, s, epsi, v, T2_coord);
		if (rank == n)
			return;
		else {
			i = 0;
			while (i < t_size) {
				if (t[i] == false) {
					if (!(binrank(Arref, J, Atight[i], n))) {
						U[i] = false;
						t[i] = true;
						q -= lambda[i];
						sumt++;
						sett++;
						Asettled.push_back(Atight[i]);
						settled_values.push_back(v[T_coord[i]] - epsi);
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
							cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
						if (disp)
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
					}
				}
				i++;
			}
			i = 0;
			while (i < t2_size) {
				if (t2[i] == false) {
					if (!(binrank(Arref, J, Atight2[i], n))) {
						U2[i] = false;
						t2[i] = true;
						sumt2++;
						sett++;
						Asettled.push_back(Atight2[i]);
						settled_values.push_back(v[T2_coord[i]]);
						if (disp)
							cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
						if (unsettled[T2_coord[i]]) {
							unsettled[T2_coord[i]] = false;
						}
						if (unsettled[s - 1 - T2_coord[i]]) {
							unsettled[s - 1 - T2_coord[i]] = false;
							if (disp)
								cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
						}
					}
				}
				i++;
			}
		}
		if (sumt == t_size)
			return;
		else {
			sr_model.remove(r);
			r = (q == 1);
			sr_model.add(r);
			if (disp)
				cout << endl << "   ---===   SOLVING SUBROUTINE LP AGAIN  ===---   " << endl << endl;
			feas = SR.solve();
			sr++;
		}
	}
	sr_env.end();
}

void subr_upd(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &q, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda, unsigned short int &rank, bool &disp, unsigned int &sett, vector<vector<bool>> &Asettled, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <double> &v, vector<unsigned int> &T2_coord) {
	i = 0;
	while (i < t_size && sumt < t_size) {
		if (t[i] == false && SR.getValue(lambda[i]) > prec) {
			U[i] = false;
			t[i] = true;
			q -= lambda[i];
			sumt++;
			if (binrank(Arref, J, Atight[i], n)) {
				if (disp)
					cout << "Rank increased to " << rank + 1 << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i] << ") getting settled." << endl;
				if (rank == n - 1) {
					rank++;
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				rowechform(Arref, J, Atight[i], n, rank);
				rank++;
				sett++;
				Asettled.push_back(Atight[i]);
				settled_values.push_back(v[T_coord[i]] - epsi);
				unsettled[T_coord[i]] = false;
				unsettled[s - 1 - T_coord[i]] = false;
				if (disp)
					cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
			}
			else {
				sett++;
				Asettled.push_back(Atight[i]);
				settled_values.push_back(v[T_coord[i]] - epsi);
				unsettled[T_coord[i]] = false;
				unsettled[s - 1 - T_coord[i]] = false;
				if (disp)
					cout << "SETTLED: " << T_coord[i] + 1 << " at " << v[T_coord[i]] - epsi << endl;
				if (disp)
					cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
	i = 0;
	while (i < t2_size && sumt2 < t2_size) {
		if (t2[i] == false && SR.getValue(lambda[i + t_size]) > prec) {
			U2[i] = false;
			t2[i] = true;
			sumt2++;
			if (binrank(Arref, J, Atight2[i], n)) {
				if (disp)
					cout << "Rank increased to " << rank + 1 << " with " << T2_coord[i] + 1 << " (and " << s - T2_coord[i] << ") getting settled." << endl;
				if (rank == n - 1) {
					rank++;
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				rowechform(Arref, J, Atight2[i], n, rank);
				rank++;
				sett++;
				Asettled.push_back(Atight2[i]);
				unsettled[T2_coord[i]] = false;
				settled_values.push_back(v[T2_coord[i]]);
				unsettled[s - 1 - T2_coord[i]] = false;
				if (disp)
					cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
			}
			else {
				sett++;
				Asettled.push_back(Atight2[i]);
				settled_values.push_back(v[T2_coord[i]]);
				unsettled[T2_coord[i]] = false;
				unsettled[s - 1 - T2_coord[i]] = false;
				if (disp)
					cout << "SETTLED: " << T2_coord[i] + 1 << " at " << v[T2_coord[i]] << endl;
				if (disp)
					cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
}

void subroutine_sg(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sett, unsigned int &sr, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <bool> &v, vector<unsigned int> &T2_coord) {
	unsigned int sumt = 0;
	vector<bool> t(t_size, false);
	unsigned short int sumt2 = 0;
	vector<bool> t2(t2_size, false);
	IloEnv sr_env;
	IloModel sr_model(sr_env);
	IloNumVarArray lambda(sr_env, t_size + t2_size + sett, 0, IloInfinity);
	lambda[t_size + t2_size].setLB(-IloInfinity);
	IloExpr sr_obj(sr_env);
	for (unsigned int i = 0; i < t_size + t2_size + sett; i++)
		sr_model.add(lambda[i]);
	IloExpr q(sr_env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(sr_env);
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				p += lambda[j];
			if (i == 0)
				q += lambda[j];
			if (j < sett && Asettled[j][i] == true)
				p += lambda[j + t_size + t2_size];
		}
		for (unsigned int j = 0; j < t2_size; j++) {
			if (Atight2[j][i] == true)
				p += lambda[j + t_size];
		}
		if (sett > t_size) {
			for (unsigned short int j = t_size; j < sett; j++) {
				if (Asettled[j][i] == true)
					p += lambda[j + t_size + t2_size];
			}
		}
		IloConstraint r = (p == 0);
		sr_model.add(r);
	}
	IloConstraint r = (q == 1);
	sr_model.add(r);
	IloObjective OBJ = IloMaximize(sr_env, sr_obj);
	sr_model.add(OBJ);
	IloCplex SR(sr_model);
	SR.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		SR.setOut(sr_env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING SUBROUTINE LP   ===---   " << endl << endl;
	bool feas = SR.solve();
	sr++;
	unsigned int i;
	while (feas) {
		subr_upd_sg(Arref, J, i, q, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, t_size, t2_size, SR, lambda, rank, disp, sett, Asettled, settled_values, unsettled, T_coord, s, epsi, v, T2_coord);
		if (rank == n)
			return;
		else {
			i = 0;
			while (i < t_size) {
				if (t[i] == false) {
					if (!(binrank(Arref, J, Atight[i], n))) {
						U[i] = false;
						t[i] = true;
						q -= lambda[i];
						sumt++;
						sett++;
						Asettled.push_back(Atight[i]);
						settled_values.push_back(v[T_coord[i]] - epsi);
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
					}
				}
				i++;
			}
			i = 0;
			while (i < t2_size) {
				if (t2[i] == false) {
					if (!(binrank(Arref, J, Atight2[i], n))) {
						U2[i] = false;
						t2[i] = true;
						sumt2++;
						sett++;
						Asettled.push_back(Atight2[i]);
						settled_values.push_back(v[T2_coord[i]]);
						if (unsettled[T2_coord[i]]) {
							unsettled[T2_coord[i]] = false;
						}
						if (unsettled[s - 1 - T2_coord[i]]) {
							unsettled[s - 1 - T2_coord[i]] = false;
							if (disp)
								cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
						}
					}
				}
				i++;
			}
		}
		if (sumt == t_size)
			return;
		else {
			sr_model.remove(r);
			r = (q == 1);
			sr_model.add(r);
			if (disp)
				cout << endl << "   ---===   SOLVING SUBROUTINE LP AGAIN  ===---   " << endl << endl;
			feas = SR.solve();
			sr++;
		}
	}
	sr_env.end();
}

void subr_upd_sg(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &q, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda, unsigned short int &rank, bool &disp, unsigned int &sett, vector<vector<bool>> &Asettled, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <bool> &v, vector<unsigned int> &T2_coord) {
	i = 0;
	while (i < t_size && sumt < t_size) {
		if (t[i] == false && SR.getValue(lambda[i]) > prec) {
			U[i] = false;
			t[i] = true;
			q -= lambda[i];
			sumt++;
			if (binrank(Arref, J, Atight[i], n)) {
				if (disp)
					cout << "Rank increased to " << rank + 1 << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i] << ") getting settled." << endl;
				if (rank == n - 1) {
					rank++;
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				rowechform(Arref, J, Atight[i], n, rank);
				rank++;
				sett++;
				Asettled.push_back(Atight[i]);
				settled_values.push_back(v[T_coord[i]] - epsi);
				unsettled[T_coord[i]] = false;
				unsettled[s - 1 - T_coord[i]] = false;
			}
			else {
				sett++;
				Asettled.push_back(Atight[i]);
				settled_values.push_back(v[T_coord[i]] - epsi);
				unsettled[T_coord[i]] = false;
				unsettled[s - 1 - T_coord[i]] = false;
				if (disp)
					cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
	i = 0;
	while (i < t2_size && sumt2 < t2_size) {
		if (t2[i] == false && SR.getValue(lambda[i + t_size]) > prec) {
			U2[i] = false;
			t2[i] = true;
			sumt2++;
			if (binrank(Arref, J, Atight2[i], n)) {
				if (disp)
					cout << "Rank increased to " << rank + 1 << " with " << T2_coord[i] + 1 << " (and " << s - T2_coord[i] << ") getting settled." << endl;
				if (rank == n - 1) {
					rank++;
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				rowechform(Arref, J, Atight2[i], n, rank);
				rank++;
				sett++;
				Asettled.push_back(Atight2[i]);
				unsettled[T2_coord[i]] = false;
				settled_values.push_back(v[T2_coord[i]]);
				unsettled[s - 1 - T2_coord[i]] = false;
			}
			else {
				sett++;
				Asettled.push_back(Atight2[i]);
				settled_values.push_back(v[T2_coord[i]]);
				unsettled[T2_coord[i]] = false;
				unsettled[s - 1 - T2_coord[i]] = false;
				if (disp)
					cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
}
