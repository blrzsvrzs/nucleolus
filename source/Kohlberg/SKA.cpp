/*
*    Nucleolus
*    SKA.cpp
*    Purpose: Verifying whether a given solution is the nucleolus of a
*             cooperative game or not using the Simplified Kohlberg
*             Algorithm
*
*    @author Marton Benedek
*    @version 1.0 19/12/2018
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

#include "SKA.h"
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
	vector<double> x(n, 0.0);
	inp.open("sol3.txt");
	for (unsigned short int i = 0; i < n; i++)
		inp >> x[i];
	inp.close();
	cout << "done!" << endl;
	if (seed == 0)
		seed = GetTickCount();
	srand(seed);
	unsigned int s = pow(2, n) - 2;
	double prec = pow(10, -6);
	vector<double> excess(s, 0);
	vector<bool> unsettled(s + 1, true);
	unsettled[s] = false;
	bool z = false;
	unsigned short int iter = 0;
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
		cout << "Running SKA..." << endl;
		double t1 = cpuTime();
		if (memo) {
			vector<bool> a(n, false);
			vector<bool> T2(n, false);
			vector<unsigned int> T2_coord(0, 0);
			unsigned short int t2_size = 0;
			vector<vector<bool>> Atight2(0, vector<bool>(n, false));
			for (unsigned int i = 0; i < n; i++) {
				if (abs(x[i] - v[pow(2, i) - 1]) < prec) {
					T2[i] = true;
					T2_coord.push_back(pow(2, i) - 1);
					a[i] = true;
					Atight2.push_back(a);
					t2_size++;
				}
			}
			excess_init_sg_mem(excess, unsettled, a, x, v, s, n);
			SKA_mem(excess, prec, s, unsettled, a, n, disp, sr, t, z, iter, t2_size, T2_coord, Atight2, t1);
		}
		else {
			vector<vector<bool>> A(s + 1, vector<bool>(n, false));
			A_mx(A, n, s);
			vector<bool> T2(n, false);
			vector<unsigned int> T2_coord(0, 0);
			unsigned short int t2_size = 0;
			vector<vector<bool>> Atight2(0, vector<bool>(n, false));
			for (unsigned int i = 0; i < n; i++) {
				if (abs(x[i] - v[pow(2, i) - 1]) < prec) {
					T2[i] = true;
					T2_coord.push_back(pow(2, i) - 1);
					Atight2.push_back(A[pow(2, i) - 1]);
					t2_size++;
				}
			}
			excess_init_sg(excess, unsettled, A, x, v, s, n);
			SKA(excess, prec, s, unsettled, A, n, disp, sr, t, z, iter, t2_size, T2_coord, Atight2, t1);
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
		cout << "Running SKA..." << endl;
		double t1 = cpuTime();
		if (memo) {
			vector<bool> a(n, false);
			vector<bool> T2(n, false);
			vector<unsigned int> T2_coord(0, 0);
			unsigned short int t2_size = 0;
			vector<vector<bool>> Atight2(0, vector<bool>(n, false));
			for (unsigned int i = 0; i < n; i++) {
				if (abs(x[i] - v[pow(2, i) - 1]) < prec) {
					T2[i] = true;
					T2_coord.push_back(pow(2, i) - 1);
					a[i] = true;
					Atight2.push_back(a);
					t2_size++;
				}
			}
			excess_init_mem(excess, unsettled, a, x, v, s, n);
			SKA_mem(excess, prec, s, unsettled, a, n, disp, sr, t, z, iter, t2_size, T2_coord, Atight2, t1);
		}
		else {
			vector<vector<bool>> A(s + 1, vector<bool>(n, false));
			A_mx(A, n, s);
			vector<bool> T2(n, false);
			vector<unsigned int> T2_coord(0, 0);
			unsigned short int t2_size = 0;
			vector<vector<bool>> Atight2(0, vector<bool>(n, false));
			for (unsigned int i = 0; i < n; i++) {
				if (abs(x[i] - v[pow(2, i) - 1]) < prec) {
					T2[i] = true;
					T2_coord.push_back(pow(2, i) - 1);
					Atight2.push_back(A[pow(2, i) - 1]);
					t2_size++;
				}
			}
			excess_init(excess, unsettled, A, x, v, s, n);
			SKA(excess, prec, s, unsettled, A, n, disp, sr, t, z, iter, t2_size, T2_coord, Atight2, t1);
		}
	}
	ofstream res;
	res.open("results.txt", ofstream::out | ofstream::trunc);
	res << seed << endl << z << endl << t << endl << iter << endl << sr << endl;
	for (unsigned int i = 0; i < n; i++)
		res << fixed << setprecision(17) << x[i] << endl;
	res.close();
	cout << "Press 0 then Enter to quit: ";
	double quit;
	cin >> quit;
	cin.get();
	return 0;
}

void SKA(vector<double> &excess, double &prec, unsigned int &s, vector<bool> &unsettled, vector<vector<bool>> &A, unsigned short int &n, bool &disp, unsigned int &sr, double &t, bool &z, unsigned short int &iter, unsigned short int &t2_size, vector<unsigned int> &T2_coord, vector<vector<bool>> &Atight2, double &t1) {
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	vector<bool>J(n, true);
	double epsi = 0;
	if (disp) {
		cout << "On the imputation boundary: ";
		for (unsigned short int i = 0; i < t2_size; i++) {
			cout << log2(T2_coord[i] + 1) + 1;
			if (i < t2_size - 1)
				cout << ", ";
		}
		cout << endl;
	}
	if (disp)
		cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	unsigned int t_size_current = 0;
	vector<vector<bool>> Atight(0, vector<bool>(n, false));
	Arref[0] = vector<double>(n, 1);
	J[0] = false;
	bool jmax = true;
	tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size, Atight, A, n);
	if (disp) {
		cout << "Number of tight coalitions: " << t_size << " (out of " << s << ")" << endl;
		if (t_size < 100) {
			for (unsigned short int i = 0; i < t_size; i++)
				cout << T_coord[i] << endl;
		}
	}
	while (jmax) {
		vector<bool> U(t_size, false);
		balanced(t_size, t2_size, Atight, Atight2, sr, Arref, J, U, n, prec, disp);
		if (disp)
			cout << endl << "  --=  subroutine finished  ==--  " << endl << endl;
		jmax = false;
		vec_maxb(jmax, J);
		iter++;
		bool u = true;
		vec_minb(u, U);
		if (!jmax && u)
			break;
		if (u) {
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					if (T[i])
						unsettled[i] = false;
					else {
						if (!binrank(Arref, J, A[i], n))
							unsettled[i] = false;
					}
				}
			}
			vec_min_uns(epsi, excess, unsettled, s);
			if (disp && iter > 0)
				cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
			if (disp)
				cout << "Epsilon: " << epsi << endl;
			t_size_current = t_size;
			tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size, Atight, A, n);
			t_size_current = t_size - t_size_current;
			if (disp) {
				cout << "Number of tight coalitions: " << t_size_current << " (" << t_size << " out of " << s << " in total)" << endl;
				if (t_size_current < 100) {
					for (unsigned short int i = t_size-t_size_current; i < t_size; i++)
						cout << T_coord[i] + 1 << endl;
				}
			}
		}
		else {
			z = false;
			t = cpuTime() - t1;
			cout << "The solution is NOT the Nucleolus." << endl;
			cout << "Time needed: " << t << " seconds" << endl;
			cout << "Iterations needed: " << iter << endl;
			cout << "Subroutine solves needed: " << sr << endl;
			return;
		}
	}
	z = true;
	t = cpuTime() - t1;
	cout << "The solution IS the Nucleolus." << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Subroutine solves needed: " << sr << endl;
	return;
}

void SKA_mem(vector<double> &excess, double &prec, unsigned int &s, vector<bool> &unsettled, vector<bool> &a, unsigned short int &n, bool &disp, unsigned int &sr, double &t, bool &z, unsigned short int &iter, unsigned short int &t2_size, vector<unsigned int> &T2_coord, vector<vector<bool>> &Atight2, double &t1) {
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	vector<bool>J(n, true);
	double epsi = 0;
	if (disp) {
		cout << "On the imputation boundary: ";
		for (unsigned short int i = 0; i < t2_size; i++) {
			cout << log2(T2_coord[i] + 1) + 1;
			if (i < t2_size - 1)
				cout << ", ";
		}
		cout << endl;
	}
	if (disp)
		cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
		cout << "Epsilon: " << epsi << endl;
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	unsigned int t_size = 0;
	unsigned int t_size_current = 0;
	vector<vector<bool>> Atight(0, vector<bool>(n, false));
	Arref[0] = vector<double>(n, 1);
	J[0] = false;
	bool jmax = true;
	tight_coal_mem(T, excess, epsi, prec, s, T_coord, unsettled, t_size, Atight, a, n);
	if (disp) {
		cout << "Number of tight coalitions: " << t_size << " (out of " << s << ")" << endl;
		if (t_size < 100) {
			for (unsigned short int i = 0; i < t_size; i++)
				cout << T_coord[i] << endl;
		}
	}
	while (jmax) {
		vector<bool> U(t_size, false);
		balanced(t_size, t2_size, Atight, Atight2, sr, Arref, J, U, n, prec, disp);
		if (disp)
			cout << endl << "  --=  subroutine finished  ==--  " << endl << endl;
		jmax = false;
		vec_maxb(jmax, J);
		iter++;
		bool u = true;
		vec_minb(u, U);
		if (!jmax && u)
			break;
		if (u) {
			for (unsigned int i = 0; i < s; i++) {
				if (unsettled[i]) {
					if (T[i])
						unsettled[i] = false;
					else {
						de2bi(i, a, n);
						if (!binrank(Arref, J, a, n))
							unsettled[i] = false;
					}
				}
			}
			vec_min_uns(epsi, excess, unsettled, s);
			if (disp && iter > 0)
				cout << endl << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl << endl;
			if (disp)
				cout << "Epsilon: " << epsi << endl;
			t_size_current = t_size;
			tight_coal_mem(T, excess, epsi, prec, s, T_coord, unsettled, t_size, Atight, a, n);
			t_size_current = t_size - t_size_current;
			if (disp) {
				cout << "Number of tight coalitions: " << t_size_current << " (" << t_size << " out of " << s << " in total)" << endl;
				if (t_size_current < 100) {
					for (unsigned short int i = t_size - t_size_current; i < t_size; i++)
						cout << T_coord[i] + 1 << endl;
				}
			}
		}
		else {
			z = false;
			t = cpuTime() - t1;
			cout << "The solution is NOT the Nucleolus." << endl;
			cout << "Time needed: " << t << " seconds" << endl;
			cout << "Iterations needed: " << iter << endl;
			cout << "Subroutine solves needed: " << sr << endl;
			return;
		}
	}
	z = true;
	t = cpuTime() - t1;
	cout << "The solution IS the Nucleolus." << endl;
	cout << "Time needed: " << t << " seconds" << endl;
	cout << "Iterations needed: " << iter << endl;
	cout << "Subroutine solves needed: " << sr << endl;
	return;
}

void balanced(unsigned int &t_size, unsigned short int &t2_size, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &sr_count, vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &U, unsigned short int &n, double &prec, bool &disp) {
	unsigned int sumt = 0;
	IloEnv env;
	IloModel model(env);
	IloNumVarArray lambda(env, t_size + t2_size, 0, IloInfinity);
	IloExpr obj(env);
	for (unsigned short int i = 0; i < n; i++) {
		IloExpr p(env);
		for (unsigned int j = 0; j < t2_size; j++) {
			if (Atight2[j][i] == true)
				p += lambda[j];
		}
		for (unsigned int j = 0; j < t_size; j++) {
			if (Atight[j][i] == true)
				p += lambda[j + t2_size];
			if (i == 0)
				obj += lambda[j + t2_size];
		}
		IloConstraint r = (p == 1);
		model.add(r);
	}
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
	sr_count++;
	unsigned int s = 0;
	unsigned int i;
	while (feas) {
		if (sr.getObjValue() < prec) {
			if (disp)
				cout << "exiting because of zero objective function value" << endl;
			env.end();
			return;
		}
		if (disp)
			cout << "subroutine update starts" << endl;
		bal_upd(t_size, t2_size, sumt, lambda, U, Arref, J, Atight, i, obj, sr, prec, n, disp);
		sum_vecb(s, J);
		if (s == 0) {
			for (unsigned int j = 0; j < t_size; j++)
				U[j] = true;
			return;
		}
		else {
			for (unsigned int j = 0; j < t_size; j++) {
				if (!U[j] && !binrank(Arref, J, Atight[j], n)) {
					U[j] = true;
					sumt++;
					obj -= lambda[j + t2_size];
				}
			}
			if (sumt == t_size)
				return;
			else {
				model.remove(OBJ);
				OBJ = IloMaximize(env, obj);
				model.add(OBJ);
				cout << endl << "  --==  solving subroutine LP again  ==--  " << endl << endl;
				feas = sr.solve();
				if (disp)
					cout << "subroutine feasibility: " << feas << endl;
				sr_count++;
			}
		}
	}
	env.end();
	return;
}

void bal_upd(unsigned int &t_size, unsigned short int &t2_size, unsigned int &sumt, IloNumVarArray &lambda, vector<bool> &U, vector<vector<double>> &Arref, vector<bool> &J, vector<vector<bool>> &Atight, unsigned int &i, IloExpr &obj, IloCplex &sr, double &prec, unsigned short int &n, bool &disp) {
	vector<vector<bool>> At(0, vector<bool>(n, false));
	i = 0;
	while (i < t_size && sumt < t_size) {
		if (!U[i] && sr.getValue(lambda[i + t2_size]) > prec) {
			U[i] = true;
			obj -= lambda[i + t2_size];
			At.push_back(Atight[i]);
			sumt++;
			if (disp)
				cout << "lambda_" << bi2de(Atight[i], n) << " > 0" << endl;
		}
		i++;
	}
	rowechform_subr(Arref, J, At, n);
	if (disp)
		cout << "subroutine update finished" << endl;
}