#include "Cell.h"

extern const vec<double> e[Q_DIM] = {{0, 0}, {1.,0}, {0,1.}, {-1.,0}, {0,-1.}, {1.,1.}, {-1.,1.}, {-1.,-1.}, {1.,-1.}};
extern const int opositeDirection[Q_DIM] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
extern const double w[Q_DIM] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};

Cell::Cell(void) {
	rho = 0.0;
	u = {0.0, 0.0};
	isObstacle = false;
	for (int l=0; l<Q_DIM; l++)
		f[l] = 0.0;
}

void Cell::updateRhoAndU(void) {
	rho = 0.0;
	u = {0.0, 0.0};
	for(int l=0; l<Q_DIM; l++) {
		rho = rho + f[l];
		u = u + e[l]*f[l];
	}
	if (rho > 0.0)
		u = u*(1.0/rho);
	else
		u = {0.0, 0.0};
}

double Cell::getFeq(int l) {
	double _eu = e[l]*u;
	return w[l]*rho*(1.0 - 1.5*(u*u) + 3.0*(_eu) + 4.5*(_eu*_eu));
}

//es medio raro pasar _rho por aca.
void Cell::setU(vec<double>& _u, double _rho) {
	rho = _rho;
	u = _u;
	for(int l=0; l<Q_DIM; l++) //guarda que saqué esto
		f[l] = getFeq(l);
}

//el uso de u es medio raro. se modifica para calcular getFeq.
//pero despues se recalcula en updateRhoAndU (medio al pedo...)
void Cell::addU(vec<double>& _u) {
	u = u + _u;
	for(int l=0; l<Q_DIM; l++) //guarda que saqué esto
		f[l] = getFeq(l);
}

void Cell::multU(double factor) {
	u = u*factor;
	for(int l=0; l<Q_DIM; l++) //guarda que saqué esto
		f[l] = getFeq(l);
}

void Cell::collide(double omega) {
	updateRhoAndU();
	for (int l = 0; l < Q_DIM; l++)
		f[l] = (1.0-omega)*f[l] + omega*getFeq(l);
}

void Cell::collide(double omega, vec<double>& _u) {
// 	updateRhoAndU();
	for (int l = 0; l < Q_DIM; l++) {
		double _eu = e[l]*_u;
		double _feq = w[l]*rho*(1.0 - 1.5*(u*u) + 3.0*(_eu) + 4.5*(_eu*_eu));
		f[l] = (1.0-omega)*f[l] + omega*_feq;
	}
}
