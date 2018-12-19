#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <time.h>
#include <ilcplex/ilocplex.h>
#include <fstream>
#include "windows.h"
#include "psapi.h"

using namespace std;

void rowechform_loop_sr2(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned int &m1, double &prec, unsigned short int &n);
void rowechform_subr2(vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n);
void vec_minb(bool &m, vector<bool> &U);
void rowechform_loop_sr(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned int &m1, unsigned int &m2, double &prec, unsigned short int &n);
void rowechform_subr(vector<vector<double>>&Arref, vector<bool> &J, vector<vector<bool>> &B, unsigned short int &n);
void tight_coal(vector<bool>&T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size, vector<vector<bool>> &Atight, vector<vector<bool>> &A, unsigned short int &n);
void tight_coal_mem(vector<bool>&T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size, vector<vector<bool>> &Atight, vector<bool> &a, unsigned short int &n);
void swap_ith_and_firstone(vector<vector<double>> &rref, vector<int> &ones, vector<int> &nonz, unsigned int &i);
void swap_ith_and_firstnz(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i);
void rowechform_piv2(vector<vector<double>> &rref, unsigned int &i, unsigned short int &n);
void rowechform_piv(vector<vector<double>> &rref, vector<int> &nonz, unsigned int &i, unsigned short int &j, unsigned int &k, unsigned short int &n);
void rowechform_loop(vector<vector<double>> &rref, vector<bool> &J, unsigned int &i, unsigned short int &j, unsigned short int &rank, double &prec, unsigned short int &n);
void rowechform(vector<vector<double>>&Arref, vector<bool> &J, vector<bool> &B, unsigned short int &n, unsigned short int &rank);
bool binrank(vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &b, unsigned short int &n);
void A_mx(vector<vector<bool>>&A, unsigned short int &n, unsigned int &s);
void de2bi(unsigned int &k, vector<bool>&a, unsigned short int &n);
void excess_init(vector<double> &exc, vector<bool> &unsettled, vector<vector<bool>> &A, vector<double> &x, vector<double> &v, unsigned int &s, unsigned short int &n);
void excess_init_mem(vector<double> &exc, vector<bool> &unsettled, vector<bool> &a, vector<double> &x, vector<double> &v, unsigned int &s, unsigned short int &n);
void excess_init_sg(vector<double> &exc, vector<bool> &unsettled, vector<vector<bool>> &A, vector<double> &x, vector<bool> &v, unsigned int &s, unsigned short int &n);
void excess_init_sg_mem(vector<double> &exc, vector<bool> &unsettled, vector<bool> &a, vector<double> &x, vector<bool> &v, unsigned int &s, unsigned short int &n);
void vec_min_uns(double&m, vector<double>&x, vector<bool> &unsettled, unsigned int &s);
void sum_vecb(unsigned int&s, vector<bool>&x);
unsigned int bi2de(vector<bool> &bin, unsigned short int &n);
void vec_maxb(bool &m, vector<bool> &U);
bool nonz_vec(vector<double>&x, double &prec);
void sc_vec_prod(vector<double>& y, double a, vector<double> &x);
void vec_subtract(vector<double>&z, vector<double>&x, vector<double>&y);
double cpuTime();