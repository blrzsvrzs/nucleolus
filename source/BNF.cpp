/*
*    Nucleolus
*    BNF.cpp
*    Purpose: finding the nucleolus of a cooperative game using
*             Benedek et al. (2018) - Finding and verifying the nucleolus
*             of cooperative games
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

#include "BNF.h"
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
	vector<double> singleton_bounds(n, 0);
	double impu = 0;
	double prec = pow(10, -6);
	vector<double> excess(s, 0);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
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
		cout << "Running BNF..." << endl;
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
			BNF_mem(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, a, t1, singleton_bounds, nlsu);
		}
		else {
			vector<vector<bool>> A(s + 1, vector<bool>(n, false));
			A_mx(A, n, s);
			excess_init_sg(excess, unsettled, A, x, v, s, n);
			BNF(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);
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
		cout << "Running BNF..." << endl;
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
			BNF_mem(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, a, t1, singleton_bounds, nlsu);
		}
		else {
			vector<vector<bool>> A(s + 1, vector<bool>(n, false));
			A_mx(A, n, s);
			excess_init(excess, unsettled, A, x, v, s, n);
			BNF(disp, n, s, excess, prec, unsettled, iter, piv, sr, t, x, A, t1, singleton_bounds, nlsu);
		}
	}
	ofstream res;
	res.open("results.txt", ofstream::out | ofstream::trunc);
	res << seed << endl << t << endl << iter << endl << piv << endl << sr << endl;
	for (unsigned int i = 0; i < n; i++)
		res << fixed << setprecision(17) << x[i] << endl;
	res.close();
	cout << "Press 0 then Enter to quit: ";
	double quit;
	cin >> quit;
	cin.get();
	return 0;
}

void BNF(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec, vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, vector<vector<bool>> &A, double &t1, vector<double> &singleton_bounds, bool &nlsu) {
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool>J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<vector<bool>> Asettled(n, vector<bool>(n, 0));
	Asettled[0] = vector<bool>(n, true);
	if (disp) {
		cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
		cout << "Starting point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
			cout << x[i] << endl;
	}
	vector<double> d(n, 0);
	double epsi = 0;
	double epsi_old = -DBL_MAX;
	while (rank < n)
		pivot(epsi, s, excess, prec, n, A, Arref, J, unsettled, rank, d, x, disp, Asettled, piv, sr, iter, unsettled_p, singleton_bounds, epsi_old, nlsu);
	t = cpuTime() - t1;
	cout << "BNF finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	cout << "Subroutine solves needed: " << sr << endl;
}

void pivot(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<vector<bool>>&A, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv, unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds, double &epsi_old, bool &nlsu) {
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned int t_size = 0;
	tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size);
	unsigned short int t2_size = 0;
	tight_coal2(T2, x, singleton_bounds, prec, n, T2_coord, unsettled_p, t2_size);
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t_size; i++)
		Atight[i] = A[T_coord[i]];
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		Atight2[i]=A[T2_coord[i]];
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	bool u = true;
	bool settled = false;
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sr_count, u, s, T_coord, T2_coord, unsettled, epsi_old, epsi, unsettled_p, settled, nlsu);
	if (disp)
		cout << endl << "  --==  subroutine finished  ==--  " << endl << endl;
	if (settled)
		iter++;
	if (disp) {
		cout << "T:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "U:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "T0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (!U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
		cout << "U0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
	}
	if (u) {
		piv++;
		if (disp)
			cout << endl << "  --==  solving improving direction LP  ==--  " << endl << endl;
		imprdir(d, n, t_size, t2_size, Atight, Atight2, U, U2, rank, Asettled, disp);
		if (disp)
			cout << endl << "  --==  improving direction obtained  ==--  " << endl << endl;
		if (disp) {
			cout << "Improving direction:" << endl;
			for (unsigned short int i = 0; i < n; i++) {
				cout << d[i] << "    ";
			}
			cout << endl;
		}
		if (disp)
			cout << endl << "  --==  computing step size  ==--  " << endl << endl;
		step(T, T2, unsettled, unsettled_p, s, A, epsi, excess, d, n, x, singleton_bounds, disp, prec);
	}
	else {
		if (disp)
			cout << endl << "  --==  minimal tight set found!  ==--  " << endl << endl;
		if (rank == n)
			return;
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
		for (unsigned short int i = 0; i < n; i++)
			if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
				unsettled_p[i] = false;
	}
	if (disp && settled)
		cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
}

void subroutine(vector<bool>&U, vector<bool>&U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &tight_size, unsigned short int &tight2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sr_count, bool &u, unsigned int &s, vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, vector<bool> &unsettled, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled, bool &nlsu) {
	unsigned int sumt = 0;
	vector<bool> t(tight_size, false);
	unsigned int sumt2 = 0;
	vector<bool> t2(tight2_size, false);
	IloEnv env;
	IloModel model(env);
	IloNumVarArray lambda(env, tight_size + tight2_size + rank, 0, IloInfinity);
	IloExpr obj(env);
	for (unsigned int i = 0; i < tight_size + tight2_size + rank; i++) {
		if (i < tight_size + tight2_size)
			obj += lambda[i];
		else
			lambda[i].setLB(-IloInfinity);
	}
	IloExpr q(env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(env);
		for (unsigned int j = 0; j < tight_size; j++) {
			if (Atight[j][i] == true)
				p += lambda[j];
			if (i == 0)
				q += lambda[j];
			if (j < rank && Asettled[j][i] == true)
				p += lambda[j + tight_size + tight2_size];
		}
		for (unsigned int j = 0; j < tight2_size; j++) {
			if (Atight2[j][i] == true)
				p += lambda[j + tight_size];
		}
		if (rank > tight_size) {
			for (unsigned short int j = tight_size; j < rank; j++) {
				if (Asettled[j][i] == true)
					p += lambda[j + tight_size + tight2_size];
			}
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
		cout << endl << "  --==  solving subroutine LP  ==--  " << endl << endl;
	bool feas = sr.solve();
	if (disp)
		cout << "subroutine feasibility: " << feas << endl;
	if (feas && nlsu)
		settled = true;
	sr_count++;
	unsigned int i;
	while (feas) {
		subr_upd(Arref, J, i, q, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, tight_size, tight2_size, sr, lambda, obj, rank, unsettled, Asettled, disp, s, T_coord, T2_coord, epsi_old, epsi, unsettled_p, settled);
		if (rank == n) {
			u = false;
			return;
		}
		if (sumt < tight_size) {
			i = 0;
			while (i < tight_size) {
				if (t[i] == false) {
					if (!(binrank(Arref, J, Atight[i], n))) {
						U[i] = false;
						t[i] = true;
						q -= lambda[i];
						obj -= lambda[i];
						sumt++;
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase." << endl;
						if (sumt == tight_size && sumt2 == tight2_size) {
							u = false;
							return;
						}
					}
				}
				i++;
			}
			i = 0;
			while (i < tight2_size) {
				if (t2[i] == false) {
					if (!(binrank(Arref, J, Atight2[i], n))) {
						U2[i] = false;
						t2[i] = true;
						obj -= lambda[i + tight_size];
						sumt2++;
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						if (disp)
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i] << " got settled without rank increase." << endl;
						if (sumt == tight_size && sumt2 == tight2_size) {
							u = false;
							return;
						}
					}
				}
				i++;
			}
			for (unsigned short int i = 0; i < n; i++)
				if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
					unsettled_p[i] = false;
			model.remove(r);
			r = (q == 1);
			model.add(r);
			model.remove(OBJ);
			OBJ = IloMinimize(env, obj);
			model.add(OBJ);
			if (disp)
				cout << endl << "  --==  solving subroutine LP again  ==--  " << endl << endl;
			feas = sr.solve();
			if (disp)
				cout << "subroutine feasibility: " << feas << endl;
			sr_count++;
		}
		else {
			u = false;
			return;
		}
	}
	env.end();
	return;
}

void subr_upd(vector<vector<double>>&Arref, vector<bool>&J, unsigned int &i, IloExpr &q, unsigned short int &n, double &prec, vector<bool>&U, vector<bool>&U2, unsigned int &sumt, unsigned int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &tight_size, unsigned short int &tight2_size, IloCplex &sr, IloNumVarArray &lambda, IloExpr &obj, unsigned short int &rank, vector<bool> &unsettled, vector<vector<bool>> &Asettled, bool &disp, unsigned int &s, vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled) {
	i = 0;
	vector<double> lambdi(tight_size+tight2_size, 0);
	for (unsigned int j = 0; j < tight_size; j++){
		if (t[j] == false)
			lambdi[j] = sr.getValue(lambda[j]);
	}
	for (unsigned short int j = 0; j < tight2_size; j++) {
		if (t2[j] == false)
			lambdi[j+tight_size] = sr.getValue(lambda[j+tight_size]);
	}
	while (i < tight_size && sumt < tight_size) {
		if (lambdi[i] > prec) {
			U[i] = false;
			t[i] = true;
			q -= lambda[i];
			obj -= lambda[i];
			sumt++;
			unsettled[T_coord[i]] = false;
			unsettled[s - 1 - T_coord[i]] = false;
			if (binrank(Arref, J, Atight[i], n)) {
				rank++;
				if (epsi > epsi_old) {
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T_coord[i] << " settled as well)" << endl;
				rowechform(Arref, J, Atight[i], n, rank);
				Asettled[rank-1]=Atight[i];
				if (rank == n) {
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				lambda[i].setLB(-IloInfinity);
			}
			else {
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, got settled (with " << s - T_coord[i] << ") without rank increase" << endl;
			}
		}
		i++;
	}
	i = 0;
	while (i < tight2_size && sumt2 < tight2_size) {
		if (lambdi[i + tight_size] > prec) {
			U2[i] = false;
			t2[i] = true;
			sumt2++;
			obj -= lambda[i + tight_size];
			unsettled[T2_coord[i]] = false;
			unsettled[s - 1 - T2_coord[i]] = false;
			if (binrank(Arref, J, Atight2[i], n)) {
				rank++;
				if (epsi > epsi_old) {
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
					cout << "lambda_" << T2_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T2_coord[i] << " settled as well)" << endl;
				rowechform(Arref, J, Atight2[i], n, rank);
				Asettled[rank-1]=Atight2[i];
				if (rank == n) {
					if (disp)
						cout << "Rank condition satisfied!" << endl;
					return;
				}
				lambda[i + tight_size].setLB(-IloInfinity);
			}
			else {
				if (disp)
					cout << "lambda_" << T_coord[i] + 1 << " > 0, got settled (with " << s - T_coord[i] << ") without rank increase" << endl;
			}
		}
		i++;
	}
	for (unsigned short int i = 0; i < n; i++)
		if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			unsettled_p[i] = false;
}

void imprdir(vector<double>&d, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<bool> &U, vector<bool> &U2, unsigned short int &rank, vector<vector<bool>>&Asettled, bool &disp) {
	IloEnv dir_env;
	IloModel dir_model(dir_env);
	IloNumVarArray D(dir_env, n, -IloInfinity, IloInfinity);
	IloExpr dir_obj(dir_env);
	for (unsigned int i = 0; i < t_size; i++) {
		IloExpr ineq(dir_env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Atight[i][j]) {
				dir_obj += D[j];
				ineq += D[j];
			}
		}
		IloConstraint r;
		if (U[i])
			r = (ineq >= 1);
		else
			r = (ineq == 0);
		dir_model.add(r);
	}
	for (unsigned int i = 0; i < t2_size; i++) {
		IloExpr ineq(dir_env);
		for (unsigned short int j = 0; j < n; j++) {
			if (Atight2[i][j]) {
				ineq += D[j];
			}
		}
		IloConstraint r;
		r = (ineq >= 0);
		dir_model.add(r);
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
	for (unsigned short int i = 0; i < n; i++) {
		d[i] = sr.getValue(D[i]);
	}
	dir_env.end();
}

void step(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s, vector<vector<bool>> &A, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec) {
	double alpha = DBL_MAX;
	double Ad;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha)
				alpha = (epsi - excess[i]) / (Ad - 1);
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (!T2[i] && unsettled_p[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha)
				alpha = (singleton_bounds[i] - x[i]) / d[i];
		}
	}
	if (disp)
		cout << "Step size: " << alpha << endl;
	if (disp)
		cout << endl << "  --==  step size obtained  ==--  " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New point: " << endl;
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

void BNF_mem(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec, vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, vector<bool> &a, double &t1, vector<double> &singleton_bounds, bool &nlsu) {
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
	double epsi_old = -DBL_MAX;
	while (rank < n)
		pivot_mem(epsi, s, excess, prec, n, a, Arref, J, unsettled, rank, d, x, disp, Asettled, piv, sr, iter, unsettled_p, singleton_bounds, epsi_old, nlsu);
	t = cpuTime() - t1;
	cout << "BNF finished!" << endl;
	cout << "The nucleolus solution:" << endl;
	for (unsigned short int i = 0; i < n; i++)
		cout << x[i] << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Pivots needed: " << piv << endl;
	cout << "Subroutine solves needed: " << sr << endl;
}

void pivot_mem(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<bool>&a, vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv, unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds, double &epsi_old, bool &nlsu) {
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned int t_size = 0;
	tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size);
	unsigned short int t2_size = 0;
	tight_coal2(T2, x, singleton_bounds, prec, n, T2_coord, unsettled_p, t2_size);
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t_size; i++)
		de2bi(T_coord[i], Atight[i], n);
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
		de2bi(T2_coord[i], Atight2[i], n);
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	bool u = true;
	bool settled = false;
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sr_count, u, s, T_coord, T2_coord, unsettled, epsi_old, epsi, unsettled_p, settled, nlsu);
	if (disp)
		cout << endl << "   ---===   SUBROUTINE FINISHED   ===---   " << endl << endl;
	if (settled)
		iter++;
	if (disp) {
		cout << "T:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (!U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "U:" << endl;
		for (unsigned int i = 0; i < t_size; i++) {
			if (U[i])
				cout << T_coord[i] + 1 << endl;
		}
		cout << "T0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (!U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
		cout << "U0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++) {
			if (U2[i])
				cout << T2_coord[i] + 1 << endl;
		}
	}
	if (u) {
		piv++;
		if (disp)
			cout << endl << "   ---===   SOLVING IMPROVING DIRECTION LP   ===---   " << endl << endl;
		imprdir(d, n, t_size, t2_size, Atight, Atight2, U, U2, rank, Asettled, disp);
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
		step_mem(T, T2, unsettled, unsettled_p, s, a, epsi, excess, d, n, x, singleton_bounds, disp, prec);
	}
	else {
		if (disp)
			cout << "Min tight set found! Rank increased to: " << rank << endl;
		if (rank == n)
			return;
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
		for (unsigned short int i = 0; i < n; i++)
			if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
				unsettled_p[i] = false;
	}
}

void step_mem(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s, vector<bool> &a, double &epsi, vector<double>&excess, vector<double> &d, unsigned short int &n, vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec) {
	double alpha = DBL_MAX;
	double Ad;
	for (unsigned int i = 0; i < s; i++) {
		if (!T[i] && unsettled[i]) {
			Ad = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha)
				alpha = (epsi - excess[i]) / (Ad - 1);
		}
	}
	for (unsigned short int i = 0; i < n; i++) {
		if (!T2[i] && unsettled_p[i]) {
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha)
				alpha = (singleton_bounds[i] - x[i]) / d[i];
		}
	}
	if (disp)
		cout << "Step size: " << alpha << endl;
	if (disp)
		cout << endl << "   ---===   STEP SIZE OBTAINED   ===---   " << endl << endl;
	for (unsigned short int i = 0; i < n; i++)
		x[i] += alpha*d[i];
	if (disp) {
		cout << "New x point: " << endl;
		for (unsigned short int i = 0; i < n; i++) {
			cout << x[i] << endl;
		}
	}
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			Ad = 0;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					Ad += d[j];
			}
			excess[i] += alpha*Ad;
		}
	}
}

void tight_coal(vector<bool>&T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (abs(excess[i] - epsi) < prec) {
				t_size++;
				T[i] = true;
				T_coord.push_back(i);
			}
		}
	}
}

void tight_coal2(vector<bool>&T2, vector<double> &x, vector<double> &singleton_bounds, double &prec, unsigned short int &n, vector<unsigned int> &T2_coord, vector<bool> &unsettled_p, unsigned short int &t2_size) {
	for (unsigned int i = 0; i < n; i++) {
		if (unsettled_p[i]) {
			if (abs(x[i] - singleton_bounds[i]) < prec) {
				t2_size++;
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
			}
		}
	}
}
