#include "common.h"

void excess_init(vector<double> &exc, vector<bool> &unsettled, vector<vector<bool>> &A, vector<double> &x, vector<double> &v, unsigned int &s, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			exc[i] = -v[i];
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					exc[i] += x[j];
			}
		}
		else
			exc[i] = DBL_MAX;
	}
}

void excess_init_sg(vector<double> &exc, vector<bool> &unsettled, vector<vector<bool>> &A, vector<double> &x, vector<bool> &v, unsigned int &s, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (v[i])
				exc[i] = -1;
			for (unsigned short int j = 0; j < n; j++) {
				if (A[i][j])
					exc[i] += x[j];
			}
		}
		else
			exc[i] = DBL_MAX;
	}
}

void excess_init_mem(vector<double> &exc, vector<bool> &unsettled, vector<bool> &a, vector<double> &x, vector<double> &v, unsigned int &s, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			exc[i] = -v[i];
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					exc[i] += x[j];
			}
		}
		else
			exc[i] = DBL_MAX;
	}
}

void excess_init_sg_mem(vector<double> &exc, vector<bool> &unsettled, vector<bool> &a, vector<double> &x, vector<bool> &v, unsigned int &s, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (v[i])
				exc[i] = -1;
			de2bi(i, a, n);
			for (unsigned short int j = 0; j < n; j++) {
				if (a[j])
					exc[i] += x[j];
			}
		}
		else
			exc[i] = DBL_MAX;
	}
}

void tight_coal(vector<bool>&T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size, vector<vector<bool>> &Atight, vector<vector<bool>> &A, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (abs(excess[i] - epsi) < prec) {
				t_size++;
				T[i] = true;
				T_coord.push_back(i);
				Atight.push_back(A[i]);
			}
		}
	}
}

void tight_coal_mem(vector<bool>&T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size, vector<vector<bool>> &Atight, vector<bool> &a, unsigned short int &n) {
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i]) {
			if (abs(excess[i] - epsi) < prec) {
				t_size++;
				T[i] = true;
				T_coord.push_back(i);
				de2bi(i, a, n);
				Atight.push_back(a);
			}
		}
	}
}

void rowechform_subr2(vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n) {
	double prec = pow(10, -10);
	unsigned int m1 = 0;
	bool size = true;
	while (size) {
		if (nonz_vec(Arref[m1], prec))
			m1++;
		else
			size = false;
	}
	vector<vector<double>> rref(m1 + 1, vector<double>(n, 0));
	for (unsigned int i = 0; i < m1; i++)
		rref[i] = Arref[i];
	for (unsigned short int j = 0; j < n; j++) {
		if (b[j])
			rref[m1][j] = 1;
	}
	for (unsigned int i = 1; i < m1 + 1; i++) {
		if (rref[i][0] > prec || rref[i][0] < -prec) {
			if (rref[i][0] < 1 + prec && rref[i][0] > 1 - prec)
				vec_subtract(rref[i], rref[0], rref[i]);
			else {
				rowechform_piv2(rref, i, n);
			}
		}
	}//first column done
	unsigned int l = 1;
	unsigned short int j = 1;
	while (l < m1 + 1 && j < n) {
		rowechform_loop_sr2(rref, J, l, j, m1, prec, n);
	}
	if (m1 + 1 < n) {
		for (unsigned short int l = 1; l < m1 + 1; l++) {
			Arref[l] = rref[l];
		}
	}
	else {
		for (unsigned short int l = 1; l < n; l++) {
			Arref[l] = rref[l];
		}
	}
}

void rowechform_subr(vector<vector<double>>&Arref, vector<bool> &J, vector<vector<bool>> &B, unsigned short int &n) {
	double prec = pow(10, -10);
	unsigned int m2 = B.size();
	unsigned int m1 = 0;
	bool size = true;
	while (size) {
		if (nonz_vec(Arref[m1], prec))
			m1++;
		else
			size = false;
	}
	vector<vector<double>> rref(m1 + m2, vector<double>(n, 0));
	rref[0] = Arref[0];// first row done
	for (unsigned int i = 1; i < m1 + m2; i++) {
		for (unsigned short int j = 0; j < n; j++) {
			if (i < m1) {
				if (Arref[i][j] > prec || Arref[i][j] < -prec)
					rref[i][j] = Arref[i][j];
			}
			else {
				if (B[i - m1][j])
					rref[i][j] = 1;
			}
		}
	}
	for (unsigned int i = 1; i < m1 + m2; i++) {
		if (rref[i][0] > prec || rref[i][0] < -prec) {
			if (rref[i][0] < 1 + prec && rref[i][0] > 1 - prec)
				vec_subtract(rref[i], rref[0], rref[i]);
			else {
				rowechform_piv2(rref, i, n);
			}
		}
	}//first column done
	unsigned int l = 1;
	unsigned short int j = 1;
	while (l < m1 + m2 && j < n) {
		rowechform_loop_sr(rref, J, l, j, m1, m2, prec, n);
	}
	if (m1 + m2 < n) {
		for (unsigned short int l = 1; l < m1 + m2; l++) {
			Arref[l] = rref[l];
		}
	}
	else {
		for (unsigned short int l = 1; l < n; l++) {
			Arref[l] = rref[l];
		}
	}
}

void rowechform_loop_sr2(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned int &m1, double &prec, unsigned short int &n) {
	vector<int> nonz(0, 0);
	vector<int> ones(0, 0);
	for (unsigned int k = i; k < m1 + 1; k++) {
		if (rref[k][j] > prec || rref[k][j] < -prec) {
			if (rref[k][j] < 1 + prec && rref[k][j] > 1 - prec) {
				ones.push_back(k);
			}
			else {
				nonz.push_back(k);
			}
		}
	}
	if (ones.size() == 0 && nonz.size() == 0) {
		j++; // all zero column, we can move on
	}
	else {
		if (ones.size() == 0) { // there are no 1s, but have non-zeros
			if (nonz[0] != i) { // if the first non-zero is not in the i-th row => Arref[i][j]=0
				swap_ith_and_firstnz(rref, nonz, i);// swap i-th and first non-zero row
			}
			sc_vec_prod(rref[i], 1 / rref[i][j], rref[i]);
		}
		else { // there are 1s => if first 1 is in the 1-th row, there's nothing to do
			if (ones[0] != i) { // if it's not in the i-th row, we swap the i-th and the first 1
				swap_ith_and_firstone(rref, ones, nonz, i);
			}
		}
		// eliminate all the pos with i-th row, then i++, j++ and move on
		if (ones.size() > 0) {
			if (ones[0] == i) {
				for (unsigned int k = 1; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
			else {
				for (unsigned int k = 0; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
		}
		if (nonz.size() > 0) {
			if (nonz[0] == i) {
				for (unsigned int k = 1; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
			else {
				for (unsigned int k = 0; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
		}
		i++;
		J[j] = false;
		j++;
	}
}

void rowechform_loop_sr(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned int &m1, unsigned int &m2, double &prec, unsigned short int &n) {
	vector<int> nonz(0, 0);
	vector<int> ones(0, 0);
	for (unsigned int k = i; k < m1 + m2; k++) {
		if (rref[k][j] > prec || rref[k][j] < -prec) {
			if (rref[k][j] < 1 + prec && rref[k][j] > 1 - prec) {
				ones.push_back(k);
			}
			else {
				nonz.push_back(k);
			}
		}
	}
	if (ones.size() == 0 && nonz.size() == 0) {
		j++; // all zero column, we can move on
	}
	else {
		if (ones.size() == 0) { // there are no 1s, but have non-zeros
			if (nonz[0] != i) { // if the first non-zero is not in the i-th row => Arref[i][j]=0
				swap_ith_and_firstnz(rref, nonz, i);// swap i-th and first non-zero row
			}
			sc_vec_prod(rref[i], 1 / rref[i][j], rref[i]);
		}
		else { // there are 1s => if first 1 is in the 1-th row, there's nothing to do
			if (ones[0] != i) { // if it's not in the i-th row, we swap the i-th and the first 1
				swap_ith_and_firstone(rref, ones, nonz, i);
			}
		}
		// eliminate all the pos with i-th row, then i++, j++ and move on
		if (ones.size() > 0) {
			if (ones[0] == i) {
				for (unsigned int k = 1; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
			else {
				for (unsigned int k = 0; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
		}
		if (nonz.size() > 0) {
			if (nonz[0] == i) {
				for (unsigned int k = 1; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
			else {
				for (unsigned int k = 0; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
		}
		i++;
		J[j] = false;
		j++;
	}
}

void vec_maxb(bool &m, vector<bool> &U) {
	for (unsigned int i = 0; i < U.size(); i++) {
		if (U[i]) {
			m = true;
			return;
		}
	}
}

void vec_minb(bool &m, vector<bool> &U) {
	for (unsigned int i = 0; i < U.size(); i++) {
		if (!U[i]) {
			m = false;
			return;
		}
	}
}

double cpuTime() {
	return (double)clock() / CLOCKS_PER_SEC;
}

void rowechform(vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &B, unsigned short int &n, unsigned short int &rank) {
	double prec = pow(10, -10);
	vector<vector<double>> rref(rank + 1, vector<double>(n, 0));
	rref[0] = Arref[0];// first row done
	for (unsigned int i = 1; i < rank + 1; i++) {
		for (unsigned short int j = 0; j < n; j++) {
			if (i < rank) {
				if (Arref[i][j] > prec || Arref[i][j] < -prec)
					rref[i][j] = Arref[i][j];
			}
			else {
				if (B[j])
					rref[i][j] = 1;
			}
		}
	}
	for (unsigned int i = 1; i < rank + 1; i++) {
		if (rref[i][0] > prec || rref[i][0] < -prec) {
			if (rref[i][0] < 1 + prec && rref[i][0] > 1 - prec)
				vec_subtract(rref[i], rref[0], rref[i]);
			else {
				rowechform_piv2(rref, i, n);
			}
		}
	}//first column done
	unsigned int l = 1;
	unsigned short int j = 1;
	while (l < rank + 1 && j < n) {
		rowechform_loop(rref, J, l, j, rank, prec, n);
	}
	if (rank + 1 < n) {
		for (unsigned short int l = 1; l < rank + 1; l++) {
			Arref[l] = rref[l];
		}
	}
	else {
		for (unsigned short int l = 1; l < n; l++) {
			Arref[l] = rref[l];
		}
	}
}

void rowechform_loop(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned short int &rank, double &prec, unsigned short int &n) {
	vector<int> nonz(0, 0);
	vector<int> ones(0, 0);
	for (unsigned int k = i; k < rank + 1; k++) {
		if (rref[k][j] > prec || rref[k][j] < -prec) {
			if (rref[k][j] < 1 + prec && rref[k][j] > 1 - prec) {
				ones.push_back(k);
			}
			else {
				nonz.push_back(k);
			}
		}
	}
	if (ones.size() == 0 && nonz.size() == 0) {
		j++; // all zero column, we can move on
	}
	else {
		if (ones.size() == 0) { // there are no 1s, but have non-zeros
			if (nonz[0] != i) { // if the first non-zero is not in the i-th row => Arref[i][j]=0
				swap_ith_and_firstnz(rref, nonz, i);// swap i-th and first non-zero row
			}
			sc_vec_prod(rref[i], 1 / rref[i][j], rref[i]);
		}
		else { // there are 1s => if first 1 is in the 1-th row, there's nothing to do
			if (ones[0] != i) { // if it's not in the i-th row, we swap the i-th and the first 1
				swap_ith_and_firstone(rref, ones, nonz, i);
			}
		}
		// eliminate all the pos with i-th row, then i++, j++ and move on
		if (ones.size() > 0) {
			if (ones[0] == i) {
				for (unsigned int k = 1; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
			else {
				for (unsigned int k = 0; k < ones.size(); k++) {
					vec_subtract(rref[ones[k]], rref[ones[k]], rref[i]);
				}
			}
		}
		if (nonz.size() > 0) {
			if (nonz[0] == i) {
				for (unsigned int k = 1; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
			else {
				for (unsigned int k = 0; k < nonz.size(); k++) {
					rowechform_piv(rref, nonz, i, j, k, n);
				}
			}
		}
		i++;
		J[j] = false;
		j++;
	}
}

void rowechform_piv(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i, unsigned short int &j, unsigned int &k, unsigned short int &n) {
	vector<double> aux(n, 0);
	sc_vec_prod(aux, rref[nonz[k]][j], rref[i]);
	vec_subtract(rref[nonz[k]], rref[nonz[k]], aux);
}

void rowechform_piv2(vector<vector<double>> &rref, unsigned int &i, unsigned short int &n) {
	vector<double> aux(n, 0);
	sc_vec_prod(aux, rref[i][0], rref[0]);
	vec_subtract(rref[i], aux, rref[i]);
}

void swap_ith_and_firstnz(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i) {
	vector<double> aux = rref[nonz[0]]; // swap i-th and first non-zero row
	rref[nonz[0]] = rref[i];
	rref[i] = aux;
	nonz[0] = i;
}

void swap_ith_and_firstone(vector<vector<double>> &rref, vector<int> &ones, vector<int> &nonz, unsigned int &i) {
	vector<double> aux = rref[ones[0]];
	rref[ones[0]] = rref[i];
	rref[i] = aux;
	if (nonz.size() > 0) {
		if (nonz[0] == i) {
			nonz[0] = ones[0];
		}
	}
	ones[0] = i;
}

bool binrank(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n) {
	double prec = pow(10, -10);
	vector<double> B(n, 0);
	for (unsigned short int i = 0; i < n; i++) {
		if (b[i] == true)
			B[i] = 1;
	}
	unsigned int m = 0;
	bool size = true;
	while (size) {
		if (nonz_vec(Arref[m], prec))
			m++;
		else
			size = false;
	}
	if (m >= n)
		return false;
	else {
		vector<bool> pivot_col(n, false);
		for (unsigned short int i = 0; i < n; i++) {
			if (J[i] == false)
				pivot_col[i] = true;
		}
		unsigned short int j = 0;
		vector<bool> piv(n, false);
		vector<double> aux(n, 0);
		unsigned short int k = 0;
		unsigned short int I = 0;
		unsigned int s = 0;
		unsigned short int ind = 0;
		unsigned short int count = 0;
		while (j < n) {
			for (unsigned short i = 0; i < n; i++) {
				if (B[i] > prec || B[i] < -prec)
					piv[i] = true;
			}
			sum_vecb(s, piv);
			if (s == 0)
				return false;
			else {
				while (k == 0) {
					if (piv[I] == true)
						k = I + 1;
					I++;
				}
				k--;
				I = 0;
				if (J[k] == true)
					return true;
				else {
					while (count < k + 1) {
						if (pivot_col[count])
							ind++;
						count++;
					}
					ind--;
					count = 0;
					sc_vec_prod(aux, B[k] / Arref[ind][k], Arref[ind]);
					vec_subtract(B, B, aux);
					j++;
				}
			}
			for (unsigned short int l = 0; l < n; l++)
				piv[l] = false;
			k = 0;
			ind = 0;
		}
		return false;
	}
}

void A_mx(vector<vector<bool>>&A, unsigned short int &n, unsigned int &s) {
	// creates boolean matrix A containing all the possible n-length boolean vectors (except for full zeros)
	for (unsigned int k = 0; k != s + 1; k++) {
		unsigned int i = 2;
		for (unsigned short int c = 0; c < n - 2; c++)
			i += i;
		unsigned int j = k + 1;
		unsigned short int l = n - 1;
		while (j > 0) {
			if (j >= i) {
				A[k][l] = true;
				j -= i;
			}
			i /= 2;
			l--;
		}
	}
}

void de2bi(unsigned int &k, vector<bool>&a, unsigned short int &n) {
	vector<bool> zero(n, false);
	a = zero;
	unsigned int i = 2;
	for (unsigned short int c = 0; c < n - 2; c++)
		i += i;
	unsigned int j = k + 1;
	unsigned short int l = n - 1;
	while (j > 0) {
		if (j >= i) {
			a[l] = true;
			j -= i;
		}
		i /= 2;
		l--;
	}
}

void vec_min_uns(double&m, vector<double>&x, vector<bool> &unsettled, unsigned int &s) {
	m = DBL_MAX;
	for (unsigned int i = 0; i < s; i++) {
		if (unsettled[i] && x[i] < m)
			m = x[i];
	}
}

void sum_vecb(unsigned int&s, vector<bool>&x) {
	// sums up the values of boolean x
	s = 0;
	for (unsigned int i = 0; i < x.size(); i++)
		s += x[i];
}

unsigned int bi2de(vector<bool> &bin, unsigned short int &n) {
	unsigned int dec = 0;
	for (unsigned short int i = 0; i < n; i++) {
		if (bin[i])
			dec += pow(2, i);
	}
	return dec;
}

void vec_maxb(bool &m, vector<bool> &U, unsigned int &t_size) {
	for (unsigned int i = 0; i < t_size; i++) {
		if (U[i]) {
			m = true;
			return;
		}
	}
}

bool nonz_vec(vector<double>&x, double &prec) {
	for (unsigned int i = 0; i < x.size(); i++) {
		if (x[i] > prec || x[i] < -prec)
			return true;
	}
	return false;
}

void vec_subtract(vector<double>&z, vector<double>&x, vector<double>&y) {
	// subtracts vector (double) y from vector (double) x
	for (unsigned int i = 0; i != x.size(); i++)
		z[i] = x[i] - y[i];
}

void sc_vec_prod(vector<double>& y, double a, vector<double> &x) {
	for (unsigned int i = 0; i < x.size(); i++)
		y[i] = a*x[i];
}