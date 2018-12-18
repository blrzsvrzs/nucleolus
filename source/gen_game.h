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

double nchoosek(unsigned short int n, unsigned short int k);
void type7(vector<double> &v, unsigned short int &n);
void type6(vector<bool> &v, unsigned int &s, unsigned short int &n);
void wvg(unsigned short int &n, vector<double> &weights, double &q, vector<bool> &v, unsigned int &s);
void type5(vector<bool> &v, unsigned int &s, unsigned short int &n);
void type4(vector<double> &v, unsigned int &s, unsigned short int &n);
void type3(vector<bool> &v, unsigned int &s, unsigned short int &n);
void type2(vector<double> &v, unsigned int &s, unsigned short int &n);
void type1(vector<double> &v, unsigned int &s, unsigned short int &n);