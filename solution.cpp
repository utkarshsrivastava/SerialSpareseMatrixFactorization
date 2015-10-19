/*
 * solution.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: afarasat
 */

#include "solution.h"
# include <iostream>
#include <cstdlib>
#include <math.h>

using namespace std;

solution::solution() {


}
solution::solution(int ** observation, double & sparsityThreshold, int & Q, int  & N, int & K) {
	_observation = observation;
	cout<<"Initial W: " << '\n';
	_Q = Q; _N = N; _K = K;
	_W = new double *[_Q];
	_C = new double *[_K];
	for (int i = 0; i < _Q; i++) {
		_W[i] = new double [_K];
		for (int j = 0; j < _K; j++) {
			_W[i][j] = 0.0;
			double randt = ((double) rand() / (RAND_MAX));
			if (randt < sparsityThreshold) {
				_W[i][j] = 5 * ((double) rand() / (RAND_MAX));
			}
			cout << _W[i][j]<< " ";
		}
		cout << '\n';
	}
	cout<<"Initial C: "<<endl;
	for (int i = 0; i < _K; i++) {
		_C[i] = new double[_N];
		for (int j = 0; j < _N; j++) {
			_C[i][j] = 2 * ((double) rand() / (RAND_MAX)) - 1;
			cout << _C[i][j]<< " ";
		}
		cout << '\n';
	}

}
solution::~solution() {
	// TODO Auto-generated destructor stub
}
double ** solution::getW(){
	return _W;
}
double ** solution::getC(){
	return _C;
}
double solution::objectiveFunction (const double & lambda, const double & mu, const double & gamma){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumWl1(0.0);  double sumWl2(0.0); double sumC(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_observation[i][j] != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i][k]*_C[k][j];
				}
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),_observation[i][j])*pow(1-(1/(1+exp(-WiCj))),(1-_observation[i][j])));
				sumLog += _observation[i][j]* log(0.001+1/(1+exp(-WiCj))) + (1-_observation[i][j])*log(1.001-1/(1+exp(-WiCj)));
			}
		}
		double sumTemp(0.0);
		for (int k = 0; k < _K; k++){

			sumWl1 += (_W[i][k] < 0) ? -_W[i][k] : _W[i][k];
			sumTemp += _W[i][k]*_W[i][k];
		}
		sumWl2 += pow(sumTemp,0.5);
	}
	for (int i = 0; i < _K; i++){
		for (int j = 0; j < _N; j++){
			sumC += _C[i][j]*_C[i][j];
		}
	}
	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2 + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::objectiveFunctionW (double  & lambda, double & mu){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumWl1(0.0);  double sumWl2(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_observation[i][j] != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i][k]*_C[k][j];
				}
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),_observation[i][j])*pow(1-(1/(1+exp(-WiCj))),(1-_observation[i][j])));
				sumLog += _observation[i][j]* log(0.001+1/(1+exp(-WiCj))) + (1-_observation[i][j])*log(1.001-1/(1+exp(-WiCj)));
			}
		}
		double sumTemp(0.0);
		for (int k = 0; k < _K; k++){

			sumWl1 += (_W[i][k] < 0) ? -_W[i][k] : _W[i][k];
			sumTemp += _W[i][k]*_W[i][k];
		}
		sumWl2 += pow(sumTemp,0.5);
	}

	objectiveFunction = -sumLog + lambda*sumWl1 + mu*sumWl2;
	return objectiveFunction;
}
double solution::objectiveFunctionC (double & gamma){
	double objectiveFunction(0.0);
	double sumLog (0.0); double sumC(0.0);
	for (int i = 0; i < _Q; i++){
		for (int j = 0; j < _N; j++){
			if (_observation[i][j] != -1){
				double WiCj(0.0);
				for (int k=0; k < _K; k++){
					WiCj += _W[i][k]*_C[k][j];
				}
				sumLog+= (_observation[i][j]* log(0.001+1/(1+exp(-WiCj)))+(1-_observation[i][j])*log(1.001-(1/(1+exp(-WiCj)))));
				//sumLog+= -log(pow(1/(1+exp(-WiCj)),observation[i][j])   *    pow(1-(1/(1+exp(-WiCj))),(1-observation[i][j])));
			}
		}
	}
	for (int i = 0; i < _K; i++){
		for (int j = 0; j < _N; j++){
			sumC += _C[i][j]*_C[i][j];
		}
	}
	objectiveFunction = -sumLog + gamma*sqrt(sumC);
	return objectiveFunction;
}
double solution::updateWij(const int & index_i, const int  & index_k, const double & mu, const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_N];
	double sum (0.0);
	for (int j = 0; j < _N; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[index_i][k]*_C[k][j];
			}
		if (_observation[index_i][j] != -1){
			yipi[j] =_observation[index_i][j] - (1/(1+exp(-WikCkj)));
		}else{
			yipi[j] =(rand()%2) - (1/(1+exp(-WikCkj)));
		}

		sum += _C[index_k][j]*yipi[j];
	}
	//cout << "Sum->W: " << sum << endl;

	delF = -sum + mu*_W[index_i][index_k];
	//cout << "DeltaF_w: " << delF << endl;
	_W[index_i][index_k] =(_W[index_i][index_k]-stepSize*delF) < 0 ? 0 : (_W[index_i][index_k]-stepSize*delF);

	return _W[index_i][index_k];
}
double solution::updateCij(const int & index_k, const int & index_j, const double & gamma,const double & stepSize){
	double delF(0.0);
	double WikCkj (0.0);
	double yipi[_Q];
	double sum (0.0);
	for (int j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[j][k]*_C[k][index_j];
			}
		if (_observation[j][index_j] != -1){
			yipi[j] = _observation[j][index_j]-(1/(1+exp(-WikCkj)));
		}else{
			yipi[j] = (rand()%2)-(1/(1+exp(-WikCkj)));
		}
		sum += _W[j][index_k]*yipi[j];
	}
	//cout << "Sum->C: " << sum << endl;
	delF = -sum + gamma*_C[index_k][index_j];
	//cout << "DeltaF_w: " << delF << endl;

	_C[index_k][index_j] = (1/(1+gamma*stepSize))*(_C[index_k][index_j]-stepSize*delF);
	//_C[index_k][index_j] = _C[index_k][index_j]-stepSize*delF;

	return _C[index_k][index_j];
}
double * solution::updateRowWij(int & index_i, double & mu,double & stepSize){
	double delF [_K];
	double WikCkj (0.0);
	double yipi[_N];
	for (int j = 0; j < _N; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[index_i][k]*_C[k][j];
			}
		if (_observation[index_i][j] != -1){
			yipi[j] =_observation[index_i][j] - (1/(1+exp(-WikCkj)));
		}else{
			yipi[j] =(rand()%2) - (1/(1+exp(-WikCkj)));
		}
	}
	double sum [_K];
	for (int k = 0; k < _K; k++){
		for (int n = 0; n < _N; n++){
			sum[k] += _C[k][n]*yipi[n];
		}
		delF[k] = -sum[k] + mu*_W[index_i][k];
		_W[index_i][k] =(_W[index_i][k]-stepSize*delF[k]) < 0 ? 0 : (_W[index_i][k]-stepSize*delF[k]);
	}

	return _W[index_i];
}

double * solution::updateColumnCij(int & index_j, double & gamma,double & stepSize){
	double delF[_K];
	double WikCkj (0.0);
	double yipi[_Q];
	for (int j = 0; j < _Q; j++){
		WikCkj = 0.0;
		for (int k = 0; k < _K; k++){
				WikCkj += _W[j][k]*_C[k][index_j];
			}
		if (_observation[j][index_j] != -1){
			yipi[j] = _observation[j][index_j]-(1/(1+exp(-WikCkj)));
		}else{
			yipi[j] = (rand()%2)-(1/(1+exp(-WikCkj)));
		}
	}
	double sum [_K];
	for (int k = 0; k < _K; k++){
		for (int n = 0; n < _Q; n++){
			sum[k] += _W[n][k]*yipi[n];
		}
		delF[k] = -sum[k] + gamma*_C[k][index_j];
		_C[k][index_j] = (1/(1+gamma*stepSize))*(_C[k][index_j]-stepSize*delF[k]);
	}

	//_C[index_k][index_j] = _C[index_k][index_j]-stepSize*delF;

	return _C[index_j];
}
