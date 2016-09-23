#ifndef CALLCHANGER_H
#define CALLCHANGER_H

#include "NewChanger.h"
#include <Eigen/Sparse>
#include "TorusSolver.h"
#include "GeneralTorus.h"
#include "shrinker.h"

#ifdef USE_CLUSTER
#include "matprod3.h"
#endif
#include "time.h"

void energy_variance();
double ph_overlap2(int Ne, int NPhi, string type, vector< vector<int> > cfl_ds, const NewChanger &control,  const Eigen::VectorXcd &vec0);
void CFL_berry();
void batch_overlap();
void orthogonality();
void PH_R2(Eigen::SparseMatrix<complex <double> > &sym);

#endif
