/*
*    Nucleolus
*    DP.cpp
*    Purpose: finding the nucleolus of a cooperative game using the
*             dual-primal sequence of Benedek (2019) - Computing the
*             nucleolus of cooperative games
*
*    @author Marton Benedek
*    @version 1.0 16/07/2019
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

#include "DP.h"
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
		cout << "Running DP..." << endl;
		if (memo)
			DP_sg_mem(disp, n, s, v, iter, piv, t, x, sr, nlsu);
		else
			DP_sg(disp, n, s, v, iter, piv, t, x, sr, nlsu);
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
		cout << "Running DP..." << endl;
		if (memo)
			DP_mem(disp, n, s, v, iter, piv, t, x, sr, nlsu);
		else
			DP(disp, n, s, v, iter, piv, t, x, sr, nlsu);
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

void DP(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x, unsigned int &sr, bool &nlsu) {
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
		balanced[i] = (bal_eq[i] == 0);
		
	}
	model.add(balanced);
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS = 0;
	while (rank < n)
		iteration(unsettled, settled, s, n, A, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values, xS, sr, nlsu);
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
	cout << "DP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values, double &xS, unsigned int &sr, bool &nlsu) {
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
	
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	vector<unsigned short int> T2_coord_p(0, 0);
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
				T2_coord_p.push_back(i);
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
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, settled, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, Lambda, obj, pos_eq, T2_coord_p, bal_eq, unsettled_p);
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
	if (!nlsu){
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void subroutine(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp, vector<unsigned int> &settled, unsigned int &sr, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <double> &v, vector<unsigned int> &T2_coord, IloNumVarArray &Lambda, IloExpr &obj, IloExpr &pos_eq, vector<unsigned short int> &T2_coord_p, IloExprArray &bal_eq, vector<bool> &unsettled_p) {
	unsigned int sumt = 0;
	vector<bool> t(t_size, false);
	unsigned short int sumt2 = 0;
	vector<bool> t2(t2_size, false);
	IloEnv sr_env;
	IloModel sr_model(sr_env);
	IloNumVarArray lambda(sr_env, t_size + t2_size + rank, 0, IloInfinity);
	vector<double> lambdi(t_size + t2_size, 0);
	IloExpr sr_obj(sr_env);
	IloRangeArray bal(sr_env, n + 1);
	IloExprArray sr_bal_eq(sr_env, n + 1);
	//for (unsigned int i = 0; i < t_size + t2_size + rank; i++)
	//	sr_model.add(lambda[i]);
	for (unsigned short int i = 0; i < rank; i++)
		lambda[t_size + t2_size + i].setLB(-IloInfinity);
	IloExpr sr_pos_eq(sr_env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr eq(sr_env);
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				eq += lambda[j];
			if (i == 0){
				sr_pos_eq += lambda[j];
				sr_obj += lambda[j];
			}
			if (j < rank && A[settled[j]][i] == true)
				p += lambda[j + t_size + t2_size];
		}
		for (unsigned int j = 0; j < t2_size; j++) {
			if (Atight2[j][i] == true)
				eq += lambda[j + t_size];
			if (i == 0)
				sr_obj += lambda[j + t_size];
		}
		if (rank > t_size) {
			for (unsigned short int j = t_size; j < rank; j++) {
				if (A[settled[j]][i] == true)
					eq += lambda[j + t_size + t2_size];
			}
		}
		sr_bal_eq[i] = eq;
		IloRange r = (eq == 0);
		bal[i] = r;
		sr_model.add(bal[i]);
	}
	sr_bal_eq[n] = sr_pos_eq;
	IloRange r = (sr_pos_eq == 1);
	bal[n] = r;
	sr_model.add(bal[n]);
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
	if (feas) {
		for (unsigned short int j = 0; j < t_size + t2_size; j++)
			lambdi[j] = SR.getValue(lambda[j]);
	}
	unsigned int i;
	while (feas) {
		subr_upd(Arref, J, i, sr_pos_eq, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, t_size, t2_size, SR, lambda, rank, disp, settled, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, Lambda, obj, pos_eq, T2_coord_p, unsettled_p, sr_bal_eq, lambdi, bal_eq);
		if (rank == n)
			return;
		else {
			i = 0;
			while (i < t_size) {
				if (t[i] == false) {
					if (!(binrank(Arref, J, Atight[i], n))) {
						U[i] = false;
						t[i] = true;
						sr_pos_eq -= lambda[i];
						sr_obj -= lambda[i];
						for (unsigned short int j = 0; j < n; j++){
							if (Atight[i][j])
								sr_bal_eq[j] -= lambda[i];
						}
						sumt++;
						unsettled[T_coord[i]] = false;
						obj -= epsi*Lambda[T_coord[i]];
						pos_eq -= Lambda[T_coord[i]];
						//unsettled[s - 1 - T_coord[i]] = false;
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
						sr_obj -= lambda[i + t_size];
						for (unsigned short int j = 0; j < n; j++) {
							if (Atight2[i][j])
								sr_bal_eq[j] -= lambda[i + t_size];
						}
						sumt2++;
						unsettled_p[T2_coord_p[i]] = false;
						if (disp)
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
					}
				}
				i++;
			}
		}
		if (sumt == t_size)
			return;
		else {
			sr_model.remove(bal[n]);
			sr_bal_eq[n] = sr_pos_eq;
			r = (sr_pos_eq == 1);
			bal[n] = r;
			sr_model.add(bal[n]);
			sr_model.remove(OBJ);
			OBJ = IloMinimize(sr_env, sr_obj);
			sr_model.add(OBJ);
			if (disp)
				cout << endl << "   ---===   SOLVING SUBROUTINE LP AGAIN  ===---   " << endl << endl;
			feas = SR.solve();
			if (disp)
				cout << "subroutine feasibility: " << feas << endl;
			sr++;
			if (feas) {
				for (unsigned short int j = 0; j < t_size + t2_size; j++)
					lambdi[j] = SR.getValue(lambda[j]);
			}
		}
	}
	sr_env.end();
}

void subr_upd(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &sr_pos_eq, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda, unsigned short int &rank, bool &disp, vector<unsigned int> &settled, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <double> &v, vector<unsigned int> &T2_coord, IloNumVarArray &Lambda, IloExpr &obj, IloExpr &pos_eq, vector<unsigned short int> &T2_coord_p, vector<bool> &unsettled_p, IloExprArray &sr_bal_eq, vector<double> &lambdi, IloExprArray &bal_eq) {
	i = 0;
	while (i < t_size && sumt < t_size) {
		if (t[i] == false && lambdi[i] > prec) {
			U[i] = false;
			t[i] = true;
			sr_pos_eq -= lambda[i];
			sr_obj -= lambda[i];
			sumt++;
			if (binrank(Arref, J, Atight[i], n)) {
				settled[rank] = T_coord[i];
				settled_values[rank] = v[T_coord[i]] - epsi;
				if (disp)
					cout << "Rank increased to " << rank + 1 << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i] << ") getting settled." << endl;
				if (rank == n - 1) {
					rank++;
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				lambda[i].setLB(-IloInfinity);
				rowechform(Arref, J, Atight[i], n, rank);
				rank++;
				unsettled[T_coord[i]] = false;
				//unsettled[s - 1 - T_coord[i]] = false;
				obj -= epsi*Lambda[T_coord[i]];
				pos_eq -= Lambda[T_coord[i]];
			}
			else {
				for (unsigned short int j = 0; j < n; j++){
					if (Atight[i][j])
						sr_bal_eq[j] -= lambda[i];
				}
				unsettled[T_coord[i]] = false;
				//unsettled[s - 1 - T_coord[i]] = false;
				obj -= epsi*Lambda[T_coord[i]];
				pos_eq -= Lambda[T_coord[i]];
				if (disp)
					cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
	i = 0;
	while (i < t2_size && sumt2 < t2_size) {
		if (t2[i] == false && lambdi[i + t_size] > prec) {
			U2[i] = false;
			t2[i] = true;
			sumt2++;
			sr_obj -= lambda[i + t_size];
			if (binrank(Arref, J, Atight2[i], n)) {
				settled[rank] = T2_coord[i];
				settled_values[rank] = v[T2_coord[i]];
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
				lambda[t_size + i].setLB(-IloInfinity);
				if (unsettled[T2_coord[i]]) {
					unsettled[T2_coord[i]] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (Atight2[i][j])
							bal_eq[j] -= Lambda[T2_coord[i]];
					}
					pos_eq -= Lambda[T2_coord[i]];
					obj -= v[T2_coord[i]] * Lambda[T2_coord[i]];
				}
				if (unsettled[s - 1 - T2_coord[i]]) {
					unsettled[s - 1 - T2_coord[i]] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (!(Atight2[i][j]))
							bal_eq[j] -= Lambda[s - 1 - T2_coord[i]];
					}
					pos_eq -= Lambda[s - 1 - T2_coord[i]];
					obj -= v[s - 1 - T2_coord[i]] * Lambda[s - 1 - T2_coord[i]];
				}
			}
			else {
				for (unsigned short int j = 0; j < n; j++) {
					if (Atight2[i][j])
						sr_bal_eq[j] -= lambda[i + t_size];
				}
				unsettled_p[T2_coord_p[i]] = false;
				if (disp)
					cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
}

void DP_sg(bool &disp, unsigned short int &n, unsigned int &s, vector<bool> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x, unsigned int &sr, bool &nlsu) {
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
		balanced[i] = (bal_eq[i] == 0);
		
	}
	model.add(balanced);
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS = 0;
	while (rank < n)
		iteration_sg(unsettled, settled, s, n, A, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values, xS, sr, nlsu);
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
	cout << "DP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_sg(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<vector<bool>> &A, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values, double &xS, unsigned int &sr, bool &nlsu) {
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
	
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	vector<unsigned short int> T2_coord_p(0, 0);
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
				T2_coord_p.push_back(i);
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
	subroutine_sg(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, settled, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, Lambda, obj, pos_eq, T2_coord_p, bal_eq, unsettled_p);
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
	if (!nlsu){
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void DP_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x, unsigned int &sr, bool &nlsu) {
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
			if (i == 0) {
				model.add(Lambda[j]);
				pos_eq += Lambda[j];
				obj += v[j] * Lambda[j];
			}
			de2bi(j, a, n);
			if (a[i])
				eq += Lambda[j];
		}
		eq += Lambda[s];
		bal_eq[i] = eq;
		balanced[i] = (bal_eq[i] == 0);
		
	}
	model.add(balanced);
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS = 0;
	while (rank < n)
		iteration_mem(unsettled, settled, s, n, a, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values, xS, sr, nlsu);
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
			de2bi([settled[i]], a, n);
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
	cout << "DP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_mem(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<double> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values, double &xS, unsigned int &sr, bool &nlsu) {
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
	for (unsigned short int i = 0; i < n; i++)
		a[i] = false;
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (u_impu[i] > prec) {
				a[i] = true;
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
							if (!(a[j]))
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
				a[i] = false;
			}
		}
	}
	
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	vector<unsigned short int> T2_coord_p(0, 0);
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
				T2_coord_p.push_back(i);
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
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, settled, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, Lambda, obj, pos_eq, T2_coord_p, bal_eq, unsettled_p);
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
	if (!nlsu){
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
						unsettled[s - 1 - i] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (!(a[j]))
								bal_eq[j] -= Lambda[s - 1 - i];
						}
						pos_eq -= Lambda[s - 1 - i];
						obj -= v[s - 1 - i] * Lambda[s - 1 - i];
					}
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}

void DP_sg_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<bool> &v, unsigned short int &iter, unsigned int &piv, double &t, vector<double> &x, unsigned int &sr, bool &nlsu) {
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
			if (i == 0) {
				model.add(Lambda[j]);
				pos_eq += Lambda[j];
				obj += v[j] * Lambda[j];
			}
			de2bi(j, a, n);
			if (a[i])
				eq += Lambda[j];
		}
		eq += Lambda[s];
		bal_eq[i] = eq;
		balanced[i] = (bal_eq[i] == 0);
		
	}
	model.add(balanced);
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "Least core solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Least core value: " << epsi << endl;
		cout << endl << "   ---===   FIRST LP SOLVED   ===---   " << endl << endl;
	}
	double xS = 0;
	while (rank < n)
		iteration_sg_mem(unsettled, settled, s, n, a, x, v, epsi, prec, Arref, J, rank, disp, model, Lambda, obj, lp, env, u, OBJ, pos_eq, pos, balanced, bal_eq, iter, piv, unsettled_p, Lambda_impu, u_impu, singleton_bounds, settled_values, xS, sr, nlsu);
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
			de2bi([settled[i]], a, n);
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
	cout << "DP finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
}

void iteration_sg_mem(vector<bool> &unsettled, vector<unsigned int> &settled, unsigned int &s, unsigned short int &n, vector<bool> &a, vector<double> &x, vector<bool> &v, double &epsi, double &prec, vector<vector<double>> &Arref, vector<bool> &J, unsigned short int &rank, bool &disp, IloModel &model, IloNumVarArray &Lambda, IloExpr &obj, IloCplex &lp, IloEnv &env, vector<double>&u, IloObjective &OBJ, IloExpr &pos_eq, IloConstraint &pos, IloRangeArray &balanced, IloExprArray &bal_eq, unsigned short int &iter, unsigned int &piv, vector<bool> &unsettled_p, IloNumVarArray &Lambda_impu, vector<double> &u_impu, vector<double> &singleton_bounds, vector<double> &settled_values, double &xS, unsigned int &sr, bool &nlsu) {
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
	for (unsigned short int i = 0; i < n; i++)
		a[i] = false;
	for (unsigned short int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (u_impu[i] > prec) {
				a[i] = true;
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
							if (!(a[j]))
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
				a[i] = false;
			}
		}
	}
	
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	vector<unsigned short int> T2_coord_p(0, 0);
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
				T2_coord_p.push_back(i);
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
	subroutine_sg(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, settled, sr, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, Lambda, obj, pos_eq, T2_coord_p, bal_eq, unsettled_p);
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
	if (!nlsu){
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
						unsettled[s - 1 - i] = false;
						for (unsigned short int j = 0; j < n; j++) {
							if (!(a[j]))
								bal_eq[j] -= Lambda[s - 1 - i];
						}
						pos_eq -= Lambda[s - 1 - i];
						obj -= v[s - 1 - i] * Lambda[s - 1 - i];
					}
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
	for (unsigned short int i = 0; i < n; i++)
		x[i] = lp.getDual(balanced[i]);
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
		cout << "New solution:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
		cout << "Epsilon: " << epsi << endl;
		cout << endl << "   ---===   " << iter << "-TH LP SOLVED   ===---   " << endl << endl;
	}
}
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
	if (!nlsu){
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

void subroutine_sg(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, unsigned short int &rank, bool &disp, vector<unsigned int> &settled, unsigned int &sr, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <bool> &v, vector<unsigned int> &T2_coord, IloNumVarArray &Lambda, IloExpr &obj, IloExpr &pos_eq, vector<unsigned short int> &T2_coord_p, IloExprArray &bal_eq, vector<bool> &unsettled_p) {
	unsigned int sumt = 0;
	vector<bool> t(t_size, false);
	unsigned short int sumt2 = 0;
	vector<bool> t2(t2_size, false);
	IloEnv sr_env;
	IloModel sr_model(sr_env);
	IloNumVarArray lambda(sr_env, t_size + t2_size + rank, 0, IloInfinity);
	vector<double> lambdi(t_size + t2_size, 0);
	IloExpr sr_obj(sr_env);
	IloRangeArray bal(sr_env, n + 1);
	IloExprArray sr_bal_eq(sr_env, n + 1);
	//for (unsigned int i = 0; i < t_size + t2_size + rank; i++)
	//	sr_model.add(lambda[i]);
	for (unsigned short int i = 0; i < rank; i++)
		lambda[t_size + t2_size + i].setLB(-IloInfinity);
	IloExpr sr_pos_eq(sr_env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr eq(sr_env);
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				eq += lambda[j];
			if (i == 0){
				sr_pos_eq += lambda[j];
				sr_obj += lambda[j];
			}
			if (j < rank && A[settled[j]][i] == true)
				p += lambda[j + t_size + t2_size];
		}
		for (unsigned int j = 0; j < t2_size; j++) {
			if (Atight2[j][i] == true)
				eq += lambda[j + t_size];
			if (i == 0)
				sr_obj += lambda[j + t_size];
		}
		if (rank > t_size) {
			for (unsigned short int j = t_size; j < rank; j++) {
				if (A[settled[j]][i] == true)
					eq += lambda[j + t_size + t2_size];
			}
		}
		sr_bal_eq[i] = eq;
		IloRange r = (eq == 0);
		bal[i] = r;
		sr_model.add(bal[i]);
	}
	sr_bal_eq[n] = sr_pos_eq;
	IloRange r = (sr_pos_eq == 1);
	bal[n] = r;
	sr_model.add(bal[n]);
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
	if (feas) {
		for (unsigned short int j = 0; j < t_size + t2_size; j++)
			lambdi[j] = SR.getValue(lambda[j]);
	}
	unsigned int i;
	while (feas) {
		subr_upd_sg(Arref, J, i, sr_pos_eq, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, t_size, t2_size, SR, lambda, rank, disp, settled, settled_values, unsettled, T_coord, s, epsi, v, T2_coord, Lambda, obj, pos_eq, T2_coord_p, unsettled_p, sr_bal_eq, lambdi, bal_eq);
		if (rank == n)
			return;
		else {
			i = 0;
			while (i < t_size) {
				if (t[i] == false) {
					if (!(binrank(Arref, J, Atight[i], n))) {
						U[i] = false;
						t[i] = true;
						sr_pos_eq -= lambda[i];
						sr_obj -= lambda[i];
						for (unsigned short int j = 0; j < n; j++){
							if (Atight[i][j])
								sr_bal_eq[j] -= lambda[i];
						}
						sumt++;
						unsettled[T_coord[i]] = false;
						obj -= epsi*Lambda[T_coord[i]];
						pos_eq -= Lambda[T_coord[i]];
						//unsettled[s - 1 - T_coord[i]] = false;
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
						sr_obj -= lambda[i + t_size];
						for (unsigned short int j = 0; j < n; j++) {
							if (Atight2[i][j])
								sr_bal_eq[j] -= lambda[i + t_size];
						}
						sumt2++;
						unsettled_p[T2_coord_p[i]] = false;
						if (disp)
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
					}
				}
				i++;
			}
		}
		if (sumt == t_size)
			return;
		else {
			sr_model.remove(bal[n]);
			sr_bal_eq[n] = sr_pos_eq;
			r = (sr_pos_eq == 1);
			bal[n] = r;
			sr_model.add(bal[n]);
			sr_model.remove(OBJ);
			OBJ = IloMinimize(sr_env, sr_obj);
			sr_model.add(OBJ);
			if (disp)
				cout << endl << "   ---===   SOLVING SUBROUTINE LP AGAIN  ===---   " << endl << endl;
			feas = SR.solve();
			if (disp)
				cout << "subroutine feasibility: " << feas << endl;
			sr++;
			if (feas) {
				for (unsigned short int j = 0; j < t_size + t2_size; j++)
					lambdi[j] = SR.getValue(lambda[j]);
			}
		}
	}
	sr_env.end();
}

void subr_upd_sg(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &sr_pos_eq, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned short int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &t_size, unsigned short int &t2_size, IloCplex &SR, IloNumVarArray &lambda, unsigned short int &rank, bool &disp, vector<unsigned int> &settled, vector <double> &settled_values, vector<bool> &unsettled, vector<unsigned int> &T_coord, unsigned int &s, double &epsi, vector <bool> &v, vector<unsigned int> &T2_coord, IloNumVarArray &Lambda, IloExpr &obj, IloExpr &pos_eq, vector<unsigned short int> &T2_coord_p, vector<bool> &unsettled_p, IloExprArray &sr_bal_eq, vector<double> &lambdi, IloExprArray &bal_eq) {
	i = 0;
	while (i < t_size && sumt < t_size) {
		if (t[i] == false && lambdi[i] > prec) {
			U[i] = false;
			t[i] = true;
			sr_pos_eq -= lambda[i];
			sr_obj -= lambda[i];
			sumt++;
			if (binrank(Arref, J, Atight[i], n)) {
				settled[rank] = T_coord[i];
				settled_values[rank] = v[T_coord[i]] - epsi;
				if (disp)
					cout << "Rank increased to " << rank + 1 << " with " << T_coord[i] + 1 << " (and " << s - T_coord[i] << ") getting settled." << endl;
				if (rank == n - 1) {
					rank++;
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				lambda[i].setLB(-IloInfinity);
				rowechform(Arref, J, Atight[i], n, rank);
				rank++;
				unsettled[T_coord[i]] = false;
				//unsettled[s - 1 - T_coord[i]] = false;
				obj -= epsi*Lambda[T_coord[i]];
				pos_eq -= Lambda[T_coord[i]];
			}
			else {
				for (unsigned short int j = 0; j < n; j++){
					if (Atight[i][j])
						sr_bal_eq[j] -= lambda[i];
				}
				unsettled[T_coord[i]] = false;
				//unsettled[s - 1 - T_coord[i]] = false;
				obj -= epsi*Lambda[T_coord[i]];
				pos_eq -= Lambda[T_coord[i]];
				if (disp)
					cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
	i = 0;
	while (i < t2_size && sumt2 < t2_size) {
		if (t2[i] == false && lambdi[i + t_size] > prec) {
			U2[i] = false;
			t2[i] = true;
			sumt2++;
			sr_obj -= lambda[i + t_size];
			if (binrank(Arref, J, Atight2[i], n)) {
				settled[rank] = T2_coord[i];
				settled_values[rank] = v[T2_coord[i]];
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
				lambda[t_size + i].setLB(-IloInfinity);
				if (unsettled[T2_coord[i]]) {
					unsettled[T2_coord[i]] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (Atight2[i][j])
							bal_eq[j] -= Lambda[T2_coord[i]];
					}
					pos_eq -= Lambda[T2_coord[i]];
					obj -= v[T2_coord[i]] * Lambda[T2_coord[i]];
				}
				if (unsettled[s - 1 - T2_coord[i]]) {
					unsettled[s - 1 - T2_coord[i]] = false;
					for (unsigned short int j = 0; j < n; j++) {
						if (!(Atight2[i][j]))
							bal_eq[j] -= Lambda[s - 1 - T2_coord[i]];
					}
					pos_eq -= Lambda[s - 1 - T2_coord[i]];
					obj -= v[s - 1 - T2_coord[i]] * Lambda[s - 1 - T2_coord[i]];
				}
			}
			else {
				for (unsigned short int j = 0; j < n; j++) {
					if (Atight2[i][j])
						sr_bal_eq[j] -= lambda[i + t_size];
				}
				unsettled_p[T2_coord_p[i]] = false;
				if (disp)
					cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
			}
		}
		i++;
	}
}