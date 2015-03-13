
#ifndef LATTICE_H
#define LATTICE_H
 
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>

//lbe2d
#include "Common.h"
#include "Cell.h"

//fftw++
#include "fftw++/Array.h"
#include "fftw++/fftw++.h"

using namespace std;
using namespace Array;
using namespace fftwpp;

typedef vector<Cell> cellVector;
typedef vector<cellVector> CellGrid;

enum Space {
	REAL_SPACE=1,
	FOURIER_SPACE
};

class Lattice {
	private:
		CellGrid grid1;
		CellGrid grid2;
		CellGrid* grid;
		CellGrid* gridTmp;

		array2<Complex> uFourierSpaceX;
		array2<Complex> uFourierSpaceY;
		array2<Complex> wFourierSpace;

		array2<double> uCoordSpaceX;
		array2<double> uCoordSpaceY;
		array2<double> wCoordSpace;
		array1<double> kx;
		array1<double> ky;
		array1<double> Ek;

		int nParticles;
		vector< vec<double> > particles;
		bool particlesAreShown;

		vec<int> dim;
		double omega;
		vec<bool> periodicBoundaries;

		void streamCell(int i, int j);
		void advectParticles();
		void swapGrids();

		//fourier
		void calcVorticity(Space space);

	public:
		Lattice();
		void init(int nX, int nY, bool periodicX, bool periodicY, double _rho, double _omega, int _nParticles);

		void update();
		void resetGrid(double _rho);
		void relaxRho();

		//seters
		void putRectangleObstacle(int x, int y, int dx, int dy);
		void setObstacle(int i, int j);
		void setU(int i, int j, vec<double> _u);
		void setU(int i, int j, vec<double> _u, double _rho);
		void addU(int i, int j, vec<double> _u);
		void showParticles();
		void hideParticles();

		//geters (hace falta que digan get?)
		const vec<int>& getDim();
		const vec<double>& getU(int i, int j);
		bool isObstacle(int i, int j);
		double getUMax();
		double getPopulation();
		double divergence();
		vec<double>& getParticlePosition(int n);
		int getNParticles();

		//fourier
		double enstrophy();
		double palinstrophy();
		double energy();
		double energy(Space space);
		double RMSvelocity();
		void normalizeLattice(double factor);
		void introRandomFourierNoise(double cValue, int kmax);
		void writeVorticity(string fileName, Space space);
		void writeSpectrum(string fileName);
};

#endif
