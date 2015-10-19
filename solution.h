/*
 * solution.h
 *
 *  Created on: Oct 10, 2015
 *      Author: afarasat
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_
#include <math.h>

class solution {
public:
	solution();
	solution(int**,double &,int &,int &,int &);
	double ** getW();
	double ** getC();
	double objectiveFunction(const double &, const double &,const double &);
	double objectiveFunctionW(double &, double &);
	double objectiveFunctionC(double &);
	double updateWij(const int &,const int &,const double &,const double &);
	double updateCij(const int &,const int &, const double &,const double &);
	double * updateRowWij(int &,double &,double &);
	double * updateColumnCij(int &,double &,double &);

	virtual ~solution();
private:
	double ** _W;
	double ** _C;
	int _Q, _N, _K;
	int ** _observation;
};

#endif /* SOLUTION_H_ */
