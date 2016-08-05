#ifndef NEW_CHANGER_H
#define NEW_CHANGER_H

#include <Eigen/Dense>
#include <Eigen/QR>
#include <vector>
#include <complex>
#include <algorithm>
#include <numeric>
#include "utils.h"

using namespace std;

extern"C"{
	void setup_z_function_table_();
	void set_l_(int *NPhi, complex<double> *l1, complex <double> *l2);
	void z_function_(double *x, double *y, complex<double> *l1, complex<double> *l2, int * rationalize, int *denom, complex<double> *z);
	complex<double> lattice_z_(int *NPhi, int *x, int *y, complex<double> *l1, complex<double> *l2, int * use_table);
}

class NewChanger{

public:
	NewChanger();
	void run();
	void test();
	void symmetry_checks();
private:
	void make_manybody_vector();
	void make_manybody_symmetry_x();
	void make_manybody_symmetry_y();
	void make_landau_symmetry_x();
	void make_landau_symmetry_y();
	void make_Amatrix();
	void make_landau_table();
	
	complex<double> get_wf(const vector< vector<int> > &zs);
	complex<double> landau_basis( int ix, int iy, int index);
	double det_helper(int z1, int z2, int d, double dbar);
	
	Eigen::VectorXcd manybody_vector;
	Eigen::MatrixXcd Amatrix, Tx_mb, Tx_lnd;
	Eigen::MatrixXcd Ty_mb, Ty_lnd;
	
	vector<double> lnd_zeros;
	vector< vector<double> > mb_zeros;
	vector< vector<int> >cfl_ds;
	vector<int> dsum;
	vector< vector< vector<int> > > mb_zs,lnd_zs;
	vector<unsigned int> lnd_states,mb_states;

	vector< vector< vector< complex<double> > > > lnd_table;
	vector< vector< complex<double> > > shifted_ztable;
	complex<double> modded_lattice_z(int x, int y);
	
	int NPhi,Ne,invNu;
	int zero,one;
	int n_mb, n_lnd,copies, ystart,ystep;
	complex<double> L1,L2;
	string type;
	int lnd_charge;
	int manybody_COM;
};
#endif

