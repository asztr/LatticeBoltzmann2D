#ifndef CELL_H
#define CELL_H
 
#include<iostream>
#include "Common.h"

#define Q_DIM 9

using namespace std;

class Cell {
	public:
		double f[Q_DIM];
		double rho;
		vec<double> u;
		bool isObstacle;

		Cell(void);
		double getFeq(int l);
		void updateRhoAndU(void);
		void collide(double omega);
		void collide(double omega, vec<double>& _u);
		void setU(vec<double>& _u, double _rho);
		void addU(vec<double>& _u);
		void multU(double factor);
};

#endif
