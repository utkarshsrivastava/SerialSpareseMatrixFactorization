/*
 * mainDriver.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: afarasat
 */
#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <cstdlib>
#include "solution.h"
#include "dataPreparation.h"

using namespace std;




int main(int arg, char* args[]) {
	// *****************************
	// Parameters Setting
	double mu (10.0);
	double lambda (1.0);
	double gamma (10.0);
	int numOfIteration (100);
	double sparcityParameter = 0.8;
	double stepSize (10);
	//*****************************
	int ** data; // full data
	int ** observation; // observation = data-missing values
	int Q(5),N(4), K(2);
	double **W;
	double ** C;

	dataPreparation dataObj(Q,N,K);
	dataObj.initialization();
	data = dataObj.getMainData();
	observation = dataObj.getObservation();
	observation = data;
	cout<<"****************"<<" OBSERVATION "<<endl;
	for (int i=0; i < Q; i++){
		for (int j =0; j < N; j++){
			cout << observation[i][j] << " ";
		}
		cout << '\n';
	}
	solution solutionObj(observation, sparcityParameter, Q, N, K);
	W = solutionObj.getW();
	C = solutionObj.getC();
	cout << "Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;
	int currentIteration = 0;
	cout<<"****************"<<" Optimization "<<endl;

	while (currentIteration++ < numOfIteration){
		int iW =(rand() % (Q-1));
		int kW = (rand() % (K-1));
		W[iW][kW] = solutionObj.updateWij(iW,kW,mu,stepSize);
		//W[iW]  = solutionObj.updateRowWij(iW,mu,stepSize);
		//cout << "W[i][k]: " << W[iW][kW] << endl;
		int kC =(rand() % (K-1)); int jC = (rand() % (N-1));
		C[kC][jC] = solutionObj.updateCij(kC,jC,gamma,stepSize);
		//C[jC] = solutionObj.updateColumnCij(jC,gamma,stepSize);
		//cout << "C[i][j]: " << C[kC][jC] << endl;
		cout << "Objective function: " << solutionObj.objectiveFunction(lambda, mu,gamma) << endl;
		stepSize = stepSize/2;
		mu = mu/2;
		gamma = gamma/2;
	}
	for (int i=0; i < Q; i++){
			for (int j =0; j < K; j++){
				cout << W[i][j] << " ";
			}
			cout << '\n';
		}
return -1;
}



