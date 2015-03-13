#include "Lattice.h"

extern const vec<double> e[Q_DIM];
extern const int opositeDirection[Q_DIM];
extern const double w[Q_DIM];
extern const Complex I(0.0, 1.0);

Lattice::Lattice() {
}

void Lattice::init(int nX, int nY, bool periodicX, bool periodicY, double _rho, double _omega, int _nParticles) {
	dim = {nX, nY};
	periodicBoundaries = {periodicX, periodicY};

	wFourierSpace.Allocate(nX, nY/2 + 1, sizeof(Complex));
	wCoordSpace.Allocate(nX, nY, sizeof(double));
	uCoordSpaceX.Allocate(nX, nY, sizeof(double));
	uCoordSpaceY.Allocate(nX, nY, sizeof(double));
	uFourierSpaceX.Allocate(nX, nY/2 + 1, sizeof(Complex));
	uFourierSpaceY.Allocate(nX, nY/2 + 1, sizeof(Complex));
	kx.Allocate(nX, sizeof(double));
	ky.Allocate(nY/2 + 1, sizeof(double));
	Ek.Allocate(nY/2 + 1, sizeof(double));

	//kx = 0,1,2,...,nx/2,-(nx/2-1),...,-2,-1
	//ky = 0,1,2,...,ny/2           
	for(int i=0; i < nX; i++)
		kx(i) = double(((i+nX/2) % nX) - nX/2);
	kx(nX/2) = double(nX/2);
	for(int j=0; j < nY/2+1; j++)
		ky(j) = double(j);

	grid1 = vector< vector<Cell> >(dim.x, vector<Cell>(dim.y));
	grid2 = vector< vector<Cell> >(dim.x, vector<Cell>(dim.y));
	grid = &grid1;
	gridTmp = &grid2;

	if (_nParticles > 0)
		particles = vector< vec<double> >(_nParticles);
	omega = _omega;
	nParticles = _nParticles;
	particlesAreShown = false;
	resetGrid(_rho);
}

void Lattice::resetGrid(double _rho) {
	vec<double> uNull = {0.0, 0.0};

	if (periodicBoundaries.x == false)
		for(int j=0; j<dim.y; j++) {
			(*grid)[0][j].isObstacle = true;
			(*grid)[dim.x-1][j].isObstacle = true;
		}
	if (periodicBoundaries.y == false)
		for(int i=0; i<dim.x; i++) {
			(*grid)[i][0].isObstacle = true;
			(*grid)[i][dim.y-1].isObstacle = true;
		}

	for (int i = 0; i < dim.x; i++)
		for (int j = 0; j < dim.y; j++)
			(*grid)[i][j].setU(uNull, _rho);

	if (nParticles > 0) {
		srand(time(NULL));
		for (int n = 0; n < nParticles; n++) {
			double _x = (int)rand()%(dim.x-1);
			double _y = (int)rand()%(dim.y-1);
			particles[n] = {_x, _y};
		}
	}
	*gridTmp = *grid;
}

double Lattice::divergence() {
	double _rmsdiv = 0.0;
	for (int i = 1; i < dim.x-1; i++)
		for (int j = 1; j < dim.y-1; j++) {
			double _dx = ((*grid)[i+1][j].u.x - (*grid)[i-1][j].u.x)/2.0;
			double _dy = ((*grid)[i][j+1].u.y - (*grid)[i][j-1].u.y)/2.0;
			double _div = _dx + _dy;
			_rmsdiv += abs(_div);
		}
	return _rmsdiv/((dim.x-2)*(dim.y-2));
}

// void Lattice::relaxRho() {
// 	array2< vec<double> > u0;
// 	//array2<double> rho0;
// 	u0.Allocate(dim.x, dim.y, sizeof(vec<double>));
// 	//rho0.Allocate(dim.x, dim.y, sizeof(double));
// 
// 	//inicializacion de rho y u
// 	for (int i = 0; i < dim.x; i++)
// 		for (int j = 0; j < dim.y; j++) {
// 			(*grid)[i][j].rho = 1.0;
// 			u0[i][j] = (*grid)[i][j].u;
// 		}
// 
// 	double rhoDiff = 1000;
// // 	while (rhoDiff > 1) {
// // 		for (int i = 0; i < dim.x; i++)
// // 			for (int j = 0; j < dim.y; j++)
// // 				if ((*grid)[i][j].isObstacle == false) {
// // 					(*grid)[i][j].collide(omega, u0[i][j]);
// // 					streamCell(i, j);
// // 					(*gridTmp)[i][j].updateRhoAndU();
// // 				}
// // 		swapGrids();
// // 		rhoDiff -= 1.0;
// // 	}
// 	while (rhoDiff > 5e-2) {
// 		//collision
// 		for (int i=0; i<dim.x; i++)
// 			for (int j=0; j<dim.y; j++)
// 				(*grid)[i][j].collide(1.0/(0.5 + 1.0/sqrt(6.0)), u0[i][j]);
// 
// 		//stream
// 		for (int i = 0; i < dim.x; i++)
// 			for (int j = 0; j < dim.y; j++)
// 				streamCell(i,j);
// 		swapGrids();
// 
// 		//calc total rho
// 		rhoDiff = 0.0;
// 		for (int i = 0; i < dim.x; i++)
// 			for (int j = 0; j < dim.y; j++) {
// 				double _rho = (*grid)[i][j].rho;
// 				(*grid)[i][j].updateRhoAndU();
// 				rhoDiff += abs(_rho - (*grid)[i][j].rho);
// 			}
// // 		cout << rhoDiff << endl;
// 	}
// 	for (int i = 0; i < dim.x; i++)
// 		for (int j = 0; j < dim.y; j++)
// 			(*grid)[i][j].setU((*grid)[i][j].u, (*grid)[i][j].rho);
// }

void Lattice::relaxRho() {
	array2< vec<double> > u0;
	array2<double> rho0;
	u0.Allocate(dim.x, dim.y, sizeof(vec<double>));
	rho0.Allocate(dim.x, dim.y, sizeof(double));

	//inicializacion de rho y u
	for (int i = 0; i < dim.x; i++)
		for (int j = 0; j < dim.y; j++) {
			rho0[i][j] = 0.0;
			u0[i][j] = (*grid)[i][j].u;
		}

	double rhoDiff = 10.0;
	while (rhoDiff > 1e-5) {
		//collision
		for (int i=0; i<dim.x; i++)
			for (int j=0; j<dim.y; j++)
				for (int l=0; l<Q_DIM; l++) {
					double _eu = e[l]*u0[i][j];
					//double _feql = w[l]*rho0[i][j]*(1.0 - 1.5*(u0[i][j]*u0[i][j]) + 3.0*(_eu) + 4.5*(_eu*_eu));
					double _feql = w[l]*(rho0[i][j] - 1.5*(u0[i][j]*u0[i][j]) + 3.0*(_eu) + 4.5*(_eu*_eu));
					double _fl = (*grid)[i][j].f[l];
					(*grid)[i][j].f[l] = _fl + omega*(_feql - _fl);//(1.0/0.908248290464)*(_feql - _fl);
				}

		//stream
		for (int i = 0; i < dim.x; i++)
			for (int j = 0; j < dim.y; j++)
				streamCell(i,j);
		swapGrids();

		//calc total rho
		rhoDiff = 0.0;
		for (int i = 0; i < dim.x; i++)
			for (int j = 0; j < dim.y; j++) {
				double _rho = 0.0;
				for (int l=0; l<Q_DIM; l++)
					_rho += (*grid)[i][j].f[l];
				rhoDiff += abs(_rho - rho0[i][j]);
				rho0[i][j] = _rho;
			}
		//cout << rhoDiff << endl;
	}
	//en teoria no necesito modificar los valores de f
	for (int i = 0; i < dim.x; i++)
		for (int j = 0; j < dim.y; j++) {
			(*grid)[i][j].rho = 1.0 + rho0[i][j];
			(*grid)[i][j].u = u0[i][j];
			for (int l=0; l<Q_DIM; l++)
				(*grid)[i][j].f[l] = (*grid)[i][j].getFeq(l);
		}
}

void Lattice::streamCell(int i, int j) {
	for (int l = 0; l < Q_DIM; l++) {
		int _opdir = opositeDirection[l];
		vec<int> fromCell = {i + (int)e[_opdir].x, j + (int)e[_opdir].y};

		if (periodicBoundaries.x == true) {
			if (fromCell.x < 0)
				fromCell.x = dim.x-1;
			else if (fromCell.x == dim.x)
				fromCell.x = 0;
		}
		if (periodicBoundaries.y == true) {
			if (fromCell.y < 0)
				fromCell.y = dim.y-1;
			else if (fromCell.y == dim.y)
				fromCell.y = 0;
		}

		if ((*grid)[fromCell.x][fromCell.y].isObstacle)
			(*gridTmp)[i][j].f[l] = (*grid)[i][j].f[_opdir];
		else
			(*gridTmp)[i][j].f[l] = (*grid)[fromCell.x][fromCell.y].f[l];
	}
}

void Lattice::advectParticles() {
// 	#pragma omp parallel for
	for (int n = 0; n < nParticles; n++) {
		int _x = (int)particles[n].x;
		int _y = (int)particles[n].y;

		if (_x <= 0 || _x >= dim.x-1 || _y <= 0 || _y >= dim.y-1) {
			particles[n].x = _x = 1 + (int)(rand() % dim.x-1);
			particles[n].y = _y = 1 + (int)(rand() % dim.y-1);
		}

		particles[n].x += 50*((*grid)[_x][_y].u.x);
		particles[n].y += 50*((*grid)[_x][_y].u.y);
	}
}

void Lattice::swapGrids() {
	//swap de los punteros a grids
	CellGrid* gridAuxPointer;
	gridAuxPointer = grid;
	grid = gridTmp;
	gridTmp = gridAuxPointer;
}

void Lattice::update() {
// 	#pragma omp parallel for
	for (int i = 0; i < dim.x; i++)
		for (int j = 0; j < dim.y; j++)
			if ((*grid)[i][j].isObstacle == false) {
				streamCell(i, j);
				(*gridTmp)[i][j].updateRhoAndU();
				(*gridTmp)[i][j].collide(omega);
			}

	swapGrids();

	if (particlesAreShown) advectParticles();
}

void Lattice::calcVorticity(Space space) {
	vec<double> _u;

	rcfft2d forwardX(dim.x, dim.y, uCoordSpaceX, uFourierSpaceX);
	rcfft2d forwardY(dim.x, dim.y, uCoordSpaceY, uFourierSpaceY);
	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y; j++) {
			_u = getU(i, j);
			uCoordSpaceX(i, j) = _u.x;
			uCoordSpaceY(i, j) = _u.y;
		}
	forwardX.fft(uCoordSpaceX, uFourierSpaceX);
	forwardY.fft(uCoordSpaceY, uFourierSpaceY);

	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y/2+1; j++) {
			uFourierSpaceX(i,j) = uFourierSpaceX(i,j)*(1./double(dim.x*dim.y));
			uFourierSpaceY(i,j) = uFourierSpaceY(i,j)*(1./double(dim.x*dim.y));
		}

	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y/2+1; j++)
			wFourierSpace(i,j) = I*(kx(i)*uFourierSpaceY(i,j) - ky(j)*uFourierSpaceX(i,j));

	if (space == REAL_SPACE) {
		crfft2d backwardW(dim.x, dim.y, wFourierSpace, wCoordSpace);
		backwardW.fftNormalized(wFourierSpace, wCoordSpace);
	}
}

void Lattice::writeVorticity(string fileName, Space space) {
	ofstream file(fileName);
	file.precision(15);

	calcVorticity(space);

	if (space == REAL_SPACE)
		for(int i=0; i < dim.x; i++) { //habria que ver si se puede mantener el orden del for de f90
			for(int j=0; j < dim.y; j++)
				file << real(wCoordSpace(i,j)) << " ";
			file << endl;
		}

	if (space == FOURIER_SPACE)
		for(int i=0; i < dim.x; i++) //habria que ver si se puede mantener el orden del for de f90
			for(int j=0; j < dim.y/2 + 1; j++)
				file << real(wFourierSpace(i,j)) << " " << imag(wFourierSpace(i,j)) << endl;
	file.close();
}

void Lattice::writeSpectrum(string fileName) {
	ofstream file(fileName);
	file.precision(15);

	calcVorticity(FOURIER_SPACE);

	for(int k=0; k < dim.y/2+1; k++)
		Ek(k) = 0.0;

	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y/2+1; j++)
			if ( ((i > 0) && (i < dim.x/2)) || (j > 0)) {
				double k2 = pow(kx(i),2.) + pow(ky(j),2.);
				double rk = sqrt(k2);
				if (rk < dim.y/2) {
					int k = (int)(sqrt(k2) + 0.5);
					Ek(k) = Ek(k) + norm(wFourierSpace(i,j))/k2;
				}
			}

	for(int k=0; k < dim.y/2+1; k++)
		file << k << " " << Ek(k) << endl;

	file.close();
}

void Lattice::introRandomFourierNoise(double cValue, int kmax) {
	vec<double> _u;
	vec<double> _k;
	vec<double> theta;

	crfft2d backwardX(dim.x, dim.y, uFourierSpaceX, uCoordSpaceX);
	crfft2d backwardY(dim.x, dim.y, uFourierSpaceY, uCoordSpaceY);

	srand(time(NULL));
	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y/2+1; j++)
			if ((i <= kmax) && (j <= kmax)) {
				_k.define((double)i/dim.x, (double)j/dim.y);
				double _uAbs = cValue/(1.0 + pow(_k.norm(),2.));
				theta.define(rand(), rand());
				theta = theta*(2.0*M_PI/RAND_MAX);
				uFourierSpaceX(i,j) = _uAbs*exp(I*theta.x);
				uFourierSpaceY(i,j) = _uAbs*exp(I*theta.y);
			}
			else {
				uFourierSpaceX(i,j) = Complex(0.0, 0.0);
				uFourierSpaceY(i,j) = Complex(0.0, 0.0);
			}

	backwardX.fft(uFourierSpaceX, uCoordSpaceX); //VER. CAPAZ QUE ES NORMALIZED
	backwardY.fft(uFourierSpaceY, uCoordSpaceY);

	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y; j++) {
			_u.define(uCoordSpaceX(i,j), uCoordSpaceY(i,j));
			addU(i, j, _u);
		}
}

void Lattice::normalizeLattice(double factor) {
	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y; j++) {
			(*grid)[i][j].multU(factor);
// 			vec<double> _u = (*grid)[i][j].u*factor;
// 			double _rho = (*grid)[i][j].rho;
// 			(*grid)[i][j].setU(_u, _rho);
		}
}

double Lattice::enstrophy() {
	calcVorticity(FOURIER_SPACE);

	double ret = 0.0;
	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y/2+1; j++)
			ret = ret + norm(wFourierSpace(i,j));
	return ret;
}

double Lattice::palinstrophy() {
	calcVorticity(FOURIER_SPACE);

	double ret = 0.0;
	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y/2+1; j++) {
			double k2 = pow(kx(i),2.) + pow(ky(j),2.);
			ret = ret + norm(wFourierSpace(i,j))*k2;
		}
	return ret;
}

double Lattice::energy() {
	double ret = 0.0;
	for(int i=0; i < dim.x; i++)
		for(int j=0; j < dim.y; j++)
			ret = ret + (*grid)[i][j].u.norm2();
	return ret/2.0;
}

double Lattice::energy(Space space) {
	//hay un factor 4 entre las energias en los espacios real y fourier
	if (space == FOURIER_SPACE) {
		calcVorticity(FOURIER_SPACE);
		double ret = 0.0;
		for(int i=0; i < dim.x; i++) {
			double _kx2 = kx(i)*kx(i);
			for(int j=0; j < dim.y/2+1; j++) {
				double _ky2 = ky(j)*ky(j);
				double _k2 = _kx2 + _ky2;
				if ( ((i > 0) && (i < dim.x/2)) || (j > 0))
					ret += norm(wFourierSpace(i,j))/_k2;
			}
		}
		return ret;
	}
	else return energy();
}

double Lattice::RMSvelocity() {
	return sqrt(2.*energy());
}

const vec<int>& Lattice::getDim() {
	return dim;
}

const vec<double>& Lattice::getU(int i, int j) {
	return (*grid)[i][j].u;
}

void Lattice::setObstacle(int i, int j) {
	(*grid)[i][j].isObstacle = true;
	(*gridTmp)[i][j].isObstacle = true;
}

void Lattice::putRectangleObstacle(int x, int y, int dx, int dy) {
	for(int i=x; i<x+dx; i++)
		for(int j=y; j<y+dy; j++) {
			(*grid)[i][j].isObstacle = true;
			(*gridTmp)[i][j].isObstacle = true;
		}
}

void Lattice::setU(int i, int j, vec<double> _u) {
	if ((i >= 0) && (i < dim.x) && (j >= 0) && (j < dim.y) && ((*grid)[i][j].isObstacle == false))
		(*grid)[i][j].setU(_u, (*grid)[i][j].rho);
}

void Lattice::setU(int i, int j, vec<double> _u, double _rho) {
	if ((i >= 0) && (i < dim.x) && (j >= 0) && (j < dim.y) && ((*grid)[i][j].isObstacle == false))
		(*grid)[i][j].setU(_u, _rho);
}
void Lattice::addU(int i, int j, vec<double> _u) {
	if ((i >= 0) && (i < dim.x) && (j >= 0) && (j < dim.y))
		(*grid)[i][j].addU(_u);
}

bool Lattice::isObstacle(int i, int j) {
	return (*grid)[i][j].isObstacle;
}

vec<double>& Lattice::getParticlePosition(int n) {
	return particles[n];
}

int Lattice::getNParticles() {
	return nParticles;
}

double Lattice::getUMax() {
	double uMax = 0;
	for (int i = 0; i < dim.x; i++)
		for (int j = 0; j < dim.y; j++)
			uMax = max(uMax, (*grid)[i][j].u.norm());
	return max(uMax, 0.0);
}

double Lattice::getPopulation() {
	double ret = 0.0;
	for (int i = 0; i < dim.x; i++)
		for (int j = 0; j < dim.y; j++)
			for (int l = 0; l < Q_DIM; l++)
				if ((*grid)[i][j].isObstacle == false) ret += (*grid)[i][j].f[l];
	return ret;
}

void Lattice::showParticles() {
	particlesAreShown = true;
}

void Lattice::hideParticles() {
	particlesAreShown = false;
}
