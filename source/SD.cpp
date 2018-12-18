/*
*    Nucleolus
*    SD.cpp
*    Purpose: finding the nucleolus of a cooperative game using the dual
*             sequence of Solymosi (1993) - On computing the nucleolus of
*             cooperative games
*             https://www.researchgate.net/publication/318017147
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

#include "SD.h"
#include "gen_game.h"

int main() {
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
		cout << "Running SD..." << endl;
		if (memo)
			SD_sg_mem(disp, n, s, v, iter, piv, t, x);
		else
			SD_sg(disp, n, s, v, iter, piv, t, x);
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
		cout << "Running SD..." << endl;
		if (memo)
			SD_mem(disp, n, s, v, iter, piv, t, x);
		else
			SD(disp, n, s, v, iter, piv, t, x);
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

void SD(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<vector<bool>> A(s + 1, vector<bool>(n, false));
	A_mx(A, n, s);
	vector<double> singleton_bounds(n, 0);
	for (unsigned short int i = 0; i < n; i++)
		singleton_bounds[i] = v[pow(2, i) - 1];
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<unsigned int> settled(n, 0);
	settled[0] = s;
	vector <double> settled_values(n, 0);
	settled_values[0] = v[s];
	vector<double> u(s + 1, 0);
	vector<double> u_impu(n, 0);
	IloEnv env;
	IloModel model(env);
	IloNumVarArray Lambda(env, s + 1, 0, IloInfinity);
	IloNumVarArray Lambda_impu(env, n, 0, IloInfinity);
	Lambda[s].setLB(-IloInfinity);
	IloRangeArray balanced(env, n);
	IloExprArray bal_eq(env, n);
	IloExpr pos_eq(env);
	IloNumVarArray impu(env, n);
	IloExpr obj(env);
	model.add(Lambda[s]);
	obj = v[s] * Lambda[s];
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr eq(env);
		model.add(Lambda_impu[i]);
		obj += v[pow(2, i) - 1] * Lambda_impu[i];
		eq += Lambda_impu[i];
		for (unsigned int j = 0; j < s; j++) {
			if (i == 0) {
				model.add(Lambda[j]);
				pos_eq += Lambda[j];
				obj += v[j] * Lambda[j];
			}
			if (A[j][i])
				eq += Lambda[j];
		}
		eq += Lambda[s];
		bal_eq[i] = eq;
		IloRange bal = (bal_eq[i] == 0);
		balanced[i] = bal;
		model.add(balanced[i]);
	}
	IloConstraint pos = (pos_eq == 1);
	model.add(pos);
	IloObjective OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	double epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	while (rank < n)
		iteration(unsettled, settled, s, n, A, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values);
	env.end();
	if (disp) {
		cout << "Settled coalitions:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << settled[i] + 1 << " at " << settled_values[i] << endl;
	}
	IloEnv sol_env;
	IloModel sol_model(sol_env);
	IloNumVarArray X(sol_env, n, -IloInfinity, IloInfinity);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(sol_env);
		for (unsigned short int j = 0; j < n; j++) {
			if (A[settled[i]][j])
				p += X[j];
		}
		IloConstraint q = (p == settled_values[i]);
		sol_model.add(q);
	}
	IloExpr sol_obj(sol_env);
	sol_model.add(IloMaximize(sol_env, sol_obj));
	IloCplex sol_lp(sol_model);
	sol_lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sol_lp.setOut(sol_env.getNullStream());
	sol_lp.solve();
	for (unsigned int i = 0; i < n; i++)
		x[i] = sol_lp.getValue(X[i]);
	sol_env.end();
	t = cpuTime() - t1;
	cout << "SD finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (u[i] > prec) {
				if (binrank(Arref, J, A[i], n)) {
					settled[rank] = i;
					settled_values[rank] = v[i] - epsi;
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << i + 1 << " getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					unsettled[i] = false;
					obj -= epsi*Lambda[i];
					pos_eq -= Lambda[i];
					rowechform(Arref, J, A[i], n, rank);
				}
				else {
					unsettled[i] = false;
					pos_eq -= Lambda[i];
					obj -= epsi * Lambda[i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (u_impu[i] > prec) {
				if (binrank(Arref, J, A[pow(2, i) - 1], n)) {
					unsettled_p[i] = false;
					settled[rank] = pow(2, i) - 1;
					settled_values[rank] = v[pow(2, i) - 1];
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << pow(2, i) << " (and " << s - pow(2, i) + 1 << ") getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					if (unsettled[pow(2, i) - 1]) {
						unsettled[pow(2, i) - 1] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (A[pow(2, i) - 1][j])
								bal_eq[j] -= Lambda[pow(2, i) - 1];
						}
						pos_eq -= Lambda[pow(2, i) - 1];
						obj -= v[pow(2, i) - 1] * Lambda[pow(2, i) - 1];
					}
					if (unsettled[s - pow(2, i)]) {
						unsettled[s - pow(2, i)] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (A[s - pow(2, i)][j])
								bal_eq[j] -= Lambda[s - pow(2, i)];
						}
						pos_eq -= Lambda[s - pow(2, i)];
						obj -= v[s - pow(2, i)] * Lambda[s - pow(2, i)];
					}
					rowechform(Arref, J, A[pow(2, i) - 1], n, rank);
				}
				else {
					unsettled_p[i] = false;
					if (disp)
						cout << pow(2, i) << " and " << s - pow(2, i) + 1 << " got settled without rank increase." << endl;
				}
			}
		}
	}
	if (disp)
		cout << "Rank increased to: " << rank << endl;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (!(binrank(Arref, J, A[i], n))) {
				unsettled[i] = false;
				for (unsigned short int j = 0; j < n; j++) {
					if (A[i][j])
						bal_eq[j] -= Lambda[i];
				}
				pos_eq -= Lambda[i];
				obj -= v[i] * Lambda[i];
				if (unsettled[s - 1 - i]) {
					unsettled[s - 1 - i] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (A[s - 1 - i][j])
							bal_eq[j] -= Lambda[s - 1 - i];
					}
					pos_eq -= Lambda[s - 1 - i];
					obj -= v[s - 1 - i] * Lambda[s - 1 - i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false) {
			bal_eq[i] -= Lambda_impu[i];
			obj -= singleton_bounds[i] * Lambda_impu[i];
			unsettled_p[i] = false;
		}
		model.remove(balanced[i]);
		balanced[i] = (bal_eq[i] == 0);
		model.add(balanced[i]);
	}
	model.remove(pos);
	pos = (pos_eq == 1);
	model.add(pos);
	model.remove(OBJ);
	OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	if (disp)
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	bool feas = lp.solve();
	if (disp)
		cout << "LP feasibility: " << feas << endl;
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights:" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void SD_sg(bool &disp, unsigned short int &n, unsigned int &s, vector<bool> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<vector<bool>> A(s + 1, vector<bool>(n, false));
	A_mx(A, n, s);
	vector<double> singleton_bounds(n, 0);
	for (unsigned short int i = 0; i < n; i++)
		singleton_bounds[i] = v[pow(2, i) - 1];
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<unsigned int> settled(n, 0);
	settled[0] = s;
	vector <double> settled_values(n, 0);
	settled_values[0] = v[s];
	vector<double> u(s + 1, 0);
	vector<double> u_impu(n, 0);
	IloEnv env;
	IloModel model(env);
	IloNumVarArray Lambda(env, s + 1, 0, IloInfinity);
	IloNumVarArray Lambda_impu(env, n, 0, IloInfinity);
	Lambda[s].setLB(-IloInfinity);
	IloRangeArray balanced(env, n);
	IloExprArray bal_eq(env, n);
	IloExpr pos_eq(env);
	IloNumVarArray impu(env, n);
	IloExpr obj(env);
	model.add(Lambda[s]);
	obj = v[s] * Lambda[s];
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr eq(env);
		model.add(Lambda_impu[i]);
		obj += v[pow(2, i) - 1] * Lambda_impu[i];
		eq += Lambda_impu[i];
		for (unsigned int j = 0; j < s; j++) {
			if (i == 0) {
				model.add(Lambda[j]);
				pos_eq += Lambda[j];
				obj += v[j] * Lambda[j];
			}
			if (A[j][i])
				eq += Lambda[j];
		}
		eq += Lambda[s];
		bal_eq[i] = eq;
		IloRange bal = (bal_eq[i] == 0);
		balanced[i] = bal;
		model.add(balanced[i]);
	}
	IloConstraint pos = (pos_eq == 1);
	model.add(pos);
	IloObjective OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	double epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	while (rank < n)
		iteration_sg(unsettled, settled, s, n, A, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values);
	env.end();
	if (disp) {
		cout << "Settled coalitions:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << settled[i] + 1 << " at " << settled_values[i] << endl;
	}
	IloEnv sol_env;
	IloModel sol_model(sol_env);
	IloNumVarArray X(sol_env, n, -IloInfinity, IloInfinity);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(sol_env);
		for (unsigned short int j = 0; j < n; j++) {
			if (A[settled[i]][j])
				p += X[j];
		}
		IloConstraint q = (p == settled_values[i]);
		sol_model.add(q);
	}
	IloExpr sol_obj(sol_env);
	sol_model.add(IloMaximize(sol_env, sol_obj));
	IloCplex sol_lp(sol_model);
	sol_lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sol_lp.setOut(sol_env.getNullStream());
	sol_lp.solve();
	for (unsigned int i = 0; i < n; i++)
		x[i] = sol_lp.getValue(X[i]);
	sol_env.end();
	t = cpuTime() - t1;
	cout << "SD finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_sg(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (u[i] > prec) {
				if (binrank(Arref, J, A[i], n)) {
					settled[rank] = i;
					settled_values[rank] = v[i] - epsi;
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << i + 1 << " getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					unsettled[i] = false;
					obj -= epsi*Lambda[i];
					pos_eq -= Lambda[i];
					rowechform(Arref, J, A[i], n, rank);
				}
				else {
					unsettled[i] = false;
					pos_eq -= Lambda[i];
					obj -= epsi * Lambda[i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (u_impu[i] > prec) {
				if (binrank(Arref, J, A[pow(2, i) - 1], n)) {
					unsettled_p[i] = false;
					settled[rank] = pow(2, i) - 1;
					settled_values[rank] = v[pow(2, i) - 1];
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << pow(2, i) << " (and " << s - pow(2, i) + 1 << ") getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					if (unsettled[pow(2, i) - 1]) {
						unsettled[pow(2, i) - 1] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (A[pow(2, i) - 1][j])
								bal_eq[j] -= Lambda[pow(2, i) - 1];
						}
						pos_eq -= Lambda[pow(2, i) - 1];
						obj -= v[pow(2, i) - 1] * Lambda[pow(2, i) - 1];
					}
					if (unsettled[s - pow(2, i)]) {
						unsettled[s - pow(2, i)] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (A[s - pow(2, i)][j])
								bal_eq[j] -= Lambda[s - pow(2, i)];
						}
						pos_eq -= Lambda[s - pow(2, i)];
						obj -= v[s - pow(2, i)] * Lambda[s - pow(2, i)];
					}
					rowechform(Arref, J, A[pow(2, i) - 1], n, rank);
				}
				else {
					unsettled_p[i] = false;
					if (disp)
						cout << pow(2, i) << " and " << s - pow(2, i) + 1 << " got settled without rank increase." << endl;
				}
			}
		}
	}
	if (disp)
		cout << "Rank increased to: " << rank << endl;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (!(binrank(Arref, J, A[i], n))) {
				unsettled[i] = false;
				for (unsigned short int j = 0; j < n; j++) {
					if (A[i][j])
						bal_eq[j] -= Lambda[i];
				}
				pos_eq -= Lambda[i];
				obj -= v[i] * Lambda[i];
				if (unsettled[s - 1 - i]) {
					unsettled[s - 1 - i] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (A[s - 1 - i][j])
							bal_eq[j] -= Lambda[s - 1 - i];
					}
					pos_eq -= Lambda[s - 1 - i];
					obj -= v[s - 1 - i] * Lambda[s - 1 - i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false) {
			bal_eq[i] -= Lambda_impu[i];
			obj -= singleton_bounds[i] * Lambda_impu[i];
			unsettled_p[i] = false;
		}
		model.remove(balanced[i]);
		balanced[i] = (bal_eq[i] == 0);
		model.add(balanced[i]);
	}
	model.remove(pos);
	pos = (pos_eq == 1);
	model.add(pos);
	model.remove(OBJ);
	OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	if (disp)
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void SD_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<bool> a(n, false);
	vector<double> singleton_bounds(n, 0);
	for (unsigned short int i = 0; i < n; i++)
		singleton_bounds[i] = v[pow(2, i) - 1];
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<unsigned int> settled(n, 0);
	settled[0] = s;
	vector <double> settled_values(n, 0);
	settled_values[0] = v[s];
	vector<double> u(s + 1, 0);
	vector<double> u_impu(n, 0);
	IloEnv env;
	IloModel model(env);
	IloNumVarArray Lambda(env, s + 1, 0, IloInfinity);
	IloNumVarArray Lambda_impu(env, n, 0, IloInfinity);
	Lambda[s].setLB(-IloInfinity);
	IloRangeArray balanced(env, n);
	IloExprArray bal_eq(env, n);
	IloExpr pos_eq(env);
	IloNumVarArray impu(env, n);
	IloExpr obj(env);
	model.add(Lambda[s]);
	obj = v[s] * Lambda[s];
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr eq(env);
		model.add(Lambda_impu[i]);
		obj += v[pow(2, i) - 1] * Lambda_impu[i];
		eq += Lambda_impu[i];		
		for (unsigned int j = 0; j < s; j++) {
			de2bi(j, a, n);
			if (i == 0) {
				model.add(Lambda[j]);
				pos_eq += Lambda[j];
				obj += v[j] * Lambda[j];
			}
			if (a[i])
				eq += Lambda[j];
		}
		eq += Lambda[s];
		bal_eq[i] = eq;
		IloRange bal = (bal_eq[i] == 0);
		balanced[i] = bal;
		model.add(balanced[i]);
	}
	IloConstraint pos = (pos_eq == 1);
	model.add(pos);
	IloObjective OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	double epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	while (rank < n)
		iteration_mem(unsettled, settled, s, n, a, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values);
	env.end();
	if (disp) {
		cout << "Settled coalitions:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << settled[i] + 1 << " at " << settled_values[i] << endl;
	}
	IloEnv sol_env;
	IloModel sol_model(sol_env);
	IloNumVarArray X(sol_env, n, -IloInfinity, IloInfinity);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(sol_env);
		de2bi(settled[i], a, n);
		for (unsigned short int j = 0; j < n; j++) {
			if (a[j])
				p += X[j];
		}
		IloConstraint q = (p == settled_values[i]);
		sol_model.add(q);
	}
	IloExpr sol_obj(sol_env);
	sol_model.add(IloMaximize(sol_env, sol_obj));
	IloCplex sol_lp(sol_model);
	sol_lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sol_lp.setOut(sol_env.getNullStream());
	sol_lp.solve();
	for (unsigned int i = 0; i < n; i++)
		x[i] = sol_lp.getValue(X[i]);
	sol_env.end();
	t = cpuTime() - t1;
	cout << "SD finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_mem(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (u[i] > prec) {
				de2bi(i, a, n);
				if (binrank(Arref, J, a, n)) {
					settled[rank] = i;
					settled_values[rank] = v[i] - epsi;
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << i + 1 << " getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					unsettled[i] = false;
					obj -= epsi*Lambda[i];
					pos_eq -= Lambda[i];
					rowechform(Arref, J, a, n, rank);
				}
				else {
					unsettled[i] = false;
					pos_eq -= Lambda[i];
					obj -= epsi * Lambda[i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (u_impu[i] > prec) {
				unsigned int idx = pow(2, i) - 1;
				de2bi(idx, a, n);
				if (binrank(Arref, J, a, n)) {
					unsettled_p[i] = false;
					settled[rank] = pow(2, i) - 1;
					settled_values[rank] = v[pow(2, i) - 1];
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << pow(2, i) << " (and " << s - pow(2, i) + 1 << ") getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					if (unsettled[pow(2, i) - 1]) {
						unsettled[pow(2, i) - 1] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (a[j])
								bal_eq[j] -= Lambda[pow(2, i) - 1];
						}
						pos_eq -= Lambda[pow(2, i) - 1];
						obj -= v[pow(2, i) - 1] * Lambda[pow(2, i) - 1];
					}
					if (unsettled[s - pow(2, i)]) {
						unsettled[s - pow(2, i)] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (a[j])
								bal_eq[j] -= Lambda[s - pow(2, i)];
						}
						pos_eq -= Lambda[s - pow(2, i)];
						obj -= v[s - pow(2, i)] * Lambda[s - pow(2, i)];
					}
					rowechform(Arref, J, a, n, rank);
				}
				else {
					unsettled_p[i] = false;
					if (disp)
						cout << pow(2, i) << " and " << s - pow(2, i) + 1 << " got settled without rank increase." << endl;
				}
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
				for (unsigned short int j = 0; j < n; j++) {
					if (a[j])
						bal_eq[j] -= Lambda[i];
				}
				pos_eq -= Lambda[i];
				obj -= v[i] * Lambda[i];
				if (unsettled[s - 1 - i]) {
					unsigned int idx = s - 1 - i;
					de2bi(idx, a, n);
					unsettled[s - 1 - i] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (a[j])
							bal_eq[j] -= Lambda[s - 1 - i];
					}
					pos_eq -= Lambda[s - 1 - i];
					obj -= v[s - 1 - i] * Lambda[s - 1 - i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false) {
			bal_eq[i] -= Lambda_impu[i];
			obj -= singleton_bounds[i] * Lambda_impu[i];
			unsettled_p[i] = false;
		}
		model.remove(balanced[i]);
		balanced[i] = (bal_eq[i] == 0);
		model.add(balanced[i]);
	}
	model.remove(pos);
	pos = (pos_eq == 1);
	model.add(pos);
	model.remove(OBJ);
	OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	if (disp)
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void SD_sg_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<bool> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x) {
	double prec = pow(10, -6);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	double t1 = cpuTime();
	vector<bool> a(n, false);
	vector<double> singleton_bounds(n, 0);
	for (unsigned short int i = 0; i < n; i++)
		singleton_bounds[i] = v[pow(2, i) - 1];
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<unsigned int> settled(n, 0);
	settled[0] = s;
	vector <double> settled_values(n, 0);
	settled_values[0] = v[s];
	vector<double> u(s + 1, 0);
	vector<double> u_impu(n, 0);
	IloEnv env;
	IloModel model(env);
	IloNumVarArray Lambda(env, s + 1, 0, IloInfinity);
	IloNumVarArray Lambda_impu(env, n, 0, IloInfinity);
	Lambda[s].setLB(-IloInfinity);
	IloRangeArray balanced(env, n);
	IloExprArray bal_eq(env, n);
	IloExpr pos_eq(env);
	IloNumVarArray impu(env, n);
	IloExpr obj(env);
	model.add(Lambda[s]);
	obj = v[s] * Lambda[s];
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr eq(env);
		model.add(Lambda_impu[i]);
		obj += v[pow(2, i) - 1] * Lambda_impu[i];
		eq += Lambda_impu[i];
		for (unsigned int j = 0; j < s; j++) {
			de2bi(j, a, n);
			if (i == 0) {
				model.add(Lambda[j]);
				pos_eq += Lambda[j];
				obj += v[j] * Lambda[j];
			}
			if (a[i])
				eq += Lambda[j];
		}
		eq += Lambda[s];
		bal_eq[i] = eq;
		IloRange bal = (bal_eq[i] == 0);
		balanced[i] = bal;
		model.add(balanced[i]);
	}
	IloConstraint pos = (pos_eq == 1);
	model.add(pos);
	IloObjective OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	IloCplex lp(model);
	lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	lp.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
	if (!disp)
		lp.setOut(env.getNullStream());
	else
		cout << endl << "   ---===   SOLVING THE FIRST LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	double epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	while (rank < n)
		iteration_sg_mem(unsettled, settled, s, n, a, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values);
	env.end();
	if (disp) {
		cout << "Settled coalitions:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << settled[i] + 1 << " at " << settled_values[i] << endl;
	}
	IloEnv sol_env;
	IloModel sol_model(sol_env);
	IloNumVarArray X(sol_env, n, -IloInfinity, IloInfinity);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(sol_env);
		de2bi(settled[i], a, n);
		for (unsigned short int j = 0; j < n; j++) {
			if (a[j])
				p += X[j];
		}
		IloConstraint q = (p == settled_values[i]);
		sol_model.add(q);
	}
	IloExpr sol_obj(sol_env);
	sol_model.add(IloMaximize(sol_env, sol_obj));
	IloCplex sol_lp(sol_model);
	sol_lp.setParam(IloCplex::Param::RootAlgorithm, 1);
	if (!disp)
		sol_lp.setOut(sol_env.getNullStream());
	sol_lp.solve();
	for (unsigned int i = 0; i < n; i++)
		x[i] = sol_lp.getValue(X[i]);
	sol_env.end();
	t = cpuTime() - t1;
	cout << "SD finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_sg_mem(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (u[i] > prec) {
				de2bi(i, a, n);
				if (binrank(Arref, J, a, n)) {
					settled[rank] = i;
					settled_values[rank] = v[i] - epsi;
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << i + 1 << " getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					unsettled[i] = false;
					obj -= epsi*Lambda[i];
					pos_eq -= Lambda[i];
					rowechform(Arref, J, a, n, rank);
				}
				else {
					unsettled[i] = false;
					pos_eq -= Lambda[i];
					obj -= epsi * Lambda[i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (u_impu[i] > prec) {
				unsigned int idx = pow(2, i) - 1;
				de2bi(idx, a, n);
				if (binrank(Arref, J, a, n)) {
					unsettled_p[i] = false;
					settled[rank] = pow(2, i) - 1;
					settled_values[rank] = v[pow(2, i) - 1];
					rank++;
					if (disp)
						cout << "Rank increased to " << rank << " with " << pow(2, i) << " (and " << s - pow(2, i) + 1 << ") getting settled." << endl;
					if (rank == n) {
						if (disp)
							cout << "Rank condition satisfied!" << endl;
						return;
					}
					if (unsettled[pow(2, i) - 1]) {
						unsettled[pow(2, i) - 1] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (a[j])
								bal_eq[j] -= Lambda[pow(2, i) - 1];
						}
						pos_eq -= Lambda[pow(2, i) - 1];
						obj -= v[pow(2, i) - 1] * Lambda[pow(2, i) - 1];
					}
					if (unsettled[s - pow(2, i)]) {
						unsettled[s - pow(2, i)] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (a[j])
								bal_eq[j] -= Lambda[s - pow(2, i)];
						}
						pos_eq -= Lambda[s - pow(2, i)];
						obj -= v[s - pow(2, i)] * Lambda[s - pow(2, i)];
					}
					rowechform(Arref, J, a, n, rank);
				}
				else {
					unsettled_p[i] = false;
					if (disp)
						cout << pow(2, i) << " and " << s - pow(2, i) + 1 << " got settled without rank increase." << endl;
				}
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
				for (unsigned short int j = 0; j < n; j++) {
					if (a[j])
						bal_eq[j] -= Lambda[i];
				}
				pos_eq -= Lambda[i];
				obj -= v[i] * Lambda[i];
				if (unsettled[s - 1 - i]) {
					unsigned int idx = s - 1 - i;
					de2bi(idx,a,n),
					unsettled[s - 1 - i] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (a[j])
							bal_eq[j] -= Lambda[s - 1 - i];
					}
					pos_eq -= Lambda[s - 1 - i];
					obj -= v[s - 1 - i] * Lambda[s - 1 - i];
				}
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false) {
			bal_eq[i] -= Lambda_impu[i];
			obj -= singleton_bounds[i] * Lambda_impu[i];
			unsettled_p[i] = false;
		}
		model.remove(balanced[i]);
		balanced[i] = (bal_eq[i] == 0);
		model.add(balanced[i]);
	}
	model.remove(pos);
	pos = (pos_eq == 1);
	model.add(pos);
	model.remove(OBJ);
	OBJ = IloMaximize(env, obj);
	model.add(OBJ);
	if (disp)
		cout << endl << "   ---===   SOLVING THE " << iter + 1 << "-TH LP   ===---   " << endl << endl;
	lp.solve();
	piv += lp.getNiterations();
	iter++;
	for (unsigned int i = 0; i < s + 1; i++)
		u[i] = lp.getValue(Lambda[i]);
	for (unsigned int i = 0; i < n; i++)
		u_impu[i] = lp.getValue(Lambda_impu[i]);
	epsi = lp.getObjValue();
	if (disp) {
		cout << "Non-zero balancing weights (least core):" << endl;
		for (unsigned short int i = 0; i < s + 1; i++) {
			if (u[i] != 0)
				cout << "u[" << i + 1 << "]=" << u[i] << endl;
		}
		cout << "With u_impu:" << endl;
		for (unsigned short int i = 0; i < n; i++) {
			if (u_impu[i] != 0)
				cout << "u_impu[" << i + 1 << "]=" << u_impu[i] << endl;
		}
		for (unsigned short int i = 0; i < n; i++)
			x[i] = lp.getDual(balanced[i]);
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}