#ifndef INTERP_H
#define INTERP_H


//standard library includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <string>
#include <cctype>
#include <ctime>
#include <sys/stat.h>
#include "stringtools.h"
#include "pTable.h"
#include "icoord.h"
#include "eckart.h"

class Interp { 

public: 
	int init(string xyzfile); 
	void structure_read(string xyzfile); 
	void calc_interp();
	void com_rotate_move(int iR, int iP, int iN, double ff);
 	double* amasses; 
	double* amasses3; 
	string* anames;               //array of atomic symbols 
  int* anumbers;                //array of atomic indices 
  int natoms;
  double** coords;
  double** allcoords;
  ICoord* icoords;
  ICoord newic;
  ICoord intic;
  ICoord int2ic;
  double* coords0;
  int* coordn;                  //coordination number
  int nimptor;
  double farBond;
 	int** bonds;
  int nbonds;
  int** angles;
  int nangles;
  int** torsions;
  int ntor;
  int** imptor;

  int nfrags;
  int* frags;
	int max_bonds;
  int max_angles;
  int max_torsions;
  int max_imptor;
  int max_nonbond;
  int n_nonbond;
  int** nonbond;
  double* nonbondd;


	double* bondd;
  double* anglev;
  double* torv;
  double* torv0;
  double* torfix;
  double* imptorv;

	void print_xyz(); 
	void freemem(); 
	void alloc_mem();
};


#endif
