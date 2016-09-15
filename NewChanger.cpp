#include "NewChanger.h"

NewChanger::NewChanger(){

}
NewChanger::NewChanger(int NPhi_t, int Ne_t, int manybody_COM_t, string type_t, vector< vector<int> >ds, double ddbarx_t, double ddbary_t, string zs_type_t):
	NPhi(NPhi_t),Ne(Ne_t),manybody_COM(manybody_COM_t),type(type_t){

	invNu=NPhi/Ne;
	double Lx=sqrt(2*M_PI*NPhi);
	double Ly=Lx;
	L1=complex<double>(Lx/sqrt(2.),0);
	L2=complex<double>(0,Ly/sqrt(2.));
	
	ddbarx=ddbarx_t;
	ddbary=ddbary_t;
	
	zero=0; //for calls to duncan's functions
	one=1; 
	
	//zeros for mb state
	mb_zeros=vector< vector<double> >(invNu,vector<double>(2));
	for(int i=0;i<invNu;i++){
		mb_zeros[i][0]=-0.5+(2*i+1)/(2.*invNu);	
//		mb_zeros[i][0]=(i+0.5)/(1.*invNu);	
		mb_zeros[i][1]=manybody_COM/(1.*invNu);
	}
	
	//CFL ds
	cfl_ds=ds;
	//vector< vector<int> > (Ne,vector<int>(2));
	dsum=vector<int>(2,0);
	for(int i=0;i<Ne;i++){
		dsum[0]+=cfl_ds[i][0]*invNu;
		dsum[1]+=cfl_ds[i][1]*invNu;
	}
//	cout<<"dsum: "<<dsum[0]<<" "<<dsum[1]<<endl;
//	cout<<endl;

	complex<double> temp;
	double x,y;
	shifted_ztable=vector< vector< complex<double> > > (NPhi,vector< complex<double> >(NPhi,0));

	for(int ix=0;ix<NPhi;ix++){
		x=(ix+dsum[0]/(1.*Ne)+ddbarx)/(1.*NPhi);
		for(int iy=0;iy<NPhi;iy++){
			y=(iy+dsum[1]/(1.*Ne)+ddbary)/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			shifted_ztable[ix][iy]=temp;
		}
	}

	
	//make the set of all possible positions, on the lnd side
	int xcharge;
	lnd_charge=Ne*manybody_COM;
	if(Ne%2==0) lnd_charge+=NPhi/2; //for some reason you dont get the charge sector you expect for even Ne
	if(type=="CFL") lnd_charge+=dsum[1]/invNu;
	lnd_charge=supermod(lnd_charge,NPhi);
//	cout<<"charge: "<<lnd_charge<<endl;
	for(int i=0;i<pow(2,NPhi);i++){
		if(count_bits(i)==Ne){
			xcharge=0;
			for(int x=0;x<NPhi;x++)
				if(i & 1<<x) xcharge+=x;
				
			if(xcharge%NPhi==lnd_charge) lnd_states.push_back(i);
//			lnd_states.push_back(i);
		}
	}
	n_lnd=lnd_states.size();
//	cout<<"lnd states"<<endl;
//	for(int i=0;i<lnd_states.size();i++) cout<<(bitset<NBITS>)lnd_states[i]<<endl;
	
	//zeros for landau basis
	lnd_zeros=vector<double>(NPhi);
	for(int i=0;i<NPhi;i++) lnd_zeros[i]=-0.5+(2*i+1)/(2.*NPhi);
	
	zs_type=zs_type_t;
	setup_mbl_zs();
	//make stuff
	set_l_(&NPhi, &L1, &L2);
	setup_z_function_table_();

}
void NewChanger::symmetry_checks(){
	make_manybody_symmetry_x();
	make_manybody_symmetry_y();
	make_landau_symmetry_x();
	make_landau_symmetry_y();
	cout<<n_mb<<" "<<n_lnd<<endl;

	complex<double> dfactor;
	if(type=="laughlin") dfactor=1.;
	else if(type=="CFL") dfactor=polar(1.,2*M_PI*dsum[1]/(1.*NPhi*invNu));
///***CHECKING X SYMMETRIES***////
	cout<<"*****CHECKING X SYMMETRIES****"<<endl;
	Eigen::VectorXcd tempX=Tx_mb*manybody_vector;
	for(int i=0;i<n_mb;i++){ 
		if(i%mb_states.size()==0) cout<<"----------y="<<i/mb_states.size()<<"---------------"<<endl;
		if(norm(manybody_vector(i)-tempX(i)*polar(1.,2*manybody_COM/(1.*invNu)*M_PI)*dfactor)>1e-10) {
			cout<<manybody_vector(i)<<" \t\t"<<tempX(i)*dfactor<<" "<<(bitset<NBITS>)mb_states[i%mb_states.size()]<<" ";
//		cout<<manybody_vector(i)<<" \t\t"<<tempX(i)<<" "<<(bitset<NBITS>)mb_states[i%mb_states.size()]<<" "<<manybody_vector(i)-tempX(i)<<endl;
			cout<<arg(manybody_vector(i)/tempX(i)/dfactor)/M_PI<<endl;
		}
	}
	cout<<endl; 

	cout<<Tx_mb.rows()<<" "<<Tx_mb.cols()<<" "<<Amatrix.rows()<<" "<<Amatrix.cols()<<" "<<Tx_lnd.rows()<<" "<<Tx_lnd.cols()<<endl;
	Eigen::MatrixXcd tempM=Tx_mb*Amatrix*(Tx_lnd.adjoint());
	for(int i=0;i<n_mb;i++){
		if(i%mb_states.size()==0) cout<<"----------y="<<i/mb_states.size()<<"---------------"<<endl;
		for(int j=0;j<n_lnd;j++){
			if(norm(Amatrix(i,j)-tempM(i,j))>1e-10){
				cout<<Amatrix(i,j)<<" \t\t"<<tempM(i,j)<<" \t\t";
				cout<<(bitset<NBITS>)mb_states[i%mb_states.size()]<<" "<<(bitset<NBITS>)lnd_states[j]<<" ";
				cout<<arg(Amatrix(i,j)/tempM(i,j))/M_PI<<endl;
			}
		}
	}

///***CHECKING Y SYMMETRIES***///
	cout<<"*****CHECKING Y SYMMETRIES****"<<endl;

	if(type=="laughlin") dfactor=1.;
	else if(type=="CFL") dfactor=polar(1.,-2*M_PI*dsum[0]/(1.*Ne*invNu));

	Eigen::VectorXcd tempY=Ty_mb*manybody_vector;
	for(int i=0;i<n_mb;i++){ 
		if(i%mb_states.size()==0) cout<<"----------y="<<i/mb_states.size()<<"---------------"<<endl;
		if(norm(manybody_vector(i)-tempY(i)*dfactor)>1e-10) {
			cout<<manybody_vector(i)<<" \t\t"<<tempY(i)*dfactor<<" "<<(bitset<NBITS>)mb_states[i%mb_states.size()]<<" ";
//		cout<<manybody_vector(i)<<" \t\t"<<tempX(i)<<" "<<(bitset<NBITS>)mb_states[i%mb_states.size()]<<" "<<manybody_vector(i)-tempX(i)<<endl;
			cout<<arg(manybody_vector(i)/tempY(i)/dfactor)/M_PI<<endl;
		}
	}
	cout<<endl; 
	
	tempM=Ty_mb*Amatrix*(Ty_lnd.adjoint());
	for(int i=0;i<n_mb;i++){
		if(i%mb_states.size()==0) cout<<"----------y="<<i/mb_states.size()<<"---------------"<<endl;
		for(int j=0;j<n_lnd;j++){
			if(norm(Amatrix(i,j)-tempM(i,j))>1e-10){
				cout<<Amatrix(i,j)<<" \t\t"<<tempM(i,j)<<" \t\t";
				cout<<(bitset<NBITS>)mb_states[i%mb_states.size()]<<" "<<(bitset<NBITS>)lnd_states[j]<<" ";
				cout<<arg(Amatrix(i,j)/tempM(i,j))/M_PI<<endl;
			}
		}
	}
}
Eigen::VectorXcd NewChanger::run(bool print, bool compute_A){
	make_manybody_vector();
	if(compute_A){
		make_landau_table();
		make_Amatrix();
	}
//	cout<<Amatrix<<endl;
//	for(int i=0;i<n_mb;i++){
//		for(int j=0; j<n_lnd; j++) cout<<abs(Amatrix(i,j))/M_PI<<" ";
//		cout<<endl;
//	}
//	cout<<endl;
//	for(int i=0;i<n_mb;i++){
//		for(int j=0; j<n_lnd; j++) cout<<arg(Amatrix(i,j))/M_PI<<" ";
//		cout<<endl;
//	}
//*****SOLUTION***///
	double outnorm;
	Eigen::VectorXcd out;
//	Eigen::VectorXcd tempB=Eigen::VectorXcd(mb_states.size());;
//	Eigen::MatrixXcd tempA;
//	for(int i=0;i<copies;i++){
//		tempA=Eigen::MatrixXcd(mb_states.size(),n_lnd);
//		for(int x=0;x<mb_states.size();x++){
//			tempB(x)=manybody_vector(x+i*mb_states.size());
//			for(int y=0;y<n_lnd;y++){
//				tempA(x,y)=Amatrix(x+i*mb_states.size(),y);
//			}
//		}
////		tempA=Amatrix.middleRows(mb_states.size()*i,mb_states.size()*(i+1));
////		tempB=manybody_vector.segment(mb_states.size()*i,mb_states.size()*(i+1));
//		cout<<tempA.rows()<<" "<<tempA.cols()<<" "<<tempB.size()<<endl;
////		out=Amatrix.middleRows(mb_states.size()*i,mb_states.size()*(i+1)).colPivHouseholderQr().solve(manybody_vector.segment(mb_states.size()*i,mb_states.size()*(i-1)));	
//		out=tempA.colPivHouseholderQr().solve(tempB);	

//		outnorm=out.norm();
//		Eigen::VectorXcd outSym=Ty_lnd*out;
//		cout<<outnorm<<endl;
//		if(!manybody_vector.isApprox(Amatrix*out,1e-9)) cout<<"bad solution!"<<endl;
//		for(int i=0;i<n_lnd;i++)
//			if(norm(out(i))>1e-16) cout<<out(i)/outnorm<<" "<<(bitset<NBITS>)lnd_states[i]<<" "<<outSym(i)/outnorm<<endl;
//	}

//	make_landau_symmetry_y();
//	make_manybody_symmetry_y();
//	
//	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es_lnd(Ty_lnd), es_mb(Ty_mb);

//	Eigen::MatrixXcd sym_A=es_mb.eigenvectors()*Amatrix*es_lnd.eigenvectors().adjoint();
//	cout<<"block diagonal?"<<endl;
//	cout<<sym_A<<endl;
//	cout<<endl;

	out=Amatrix.colPivHouseholderQr().solve(manybody_vector);	
	if(!manybody_vector.isApprox(Amatrix*out,1e-4)) cout<<"*******************bad solution!****************"<<endl;
	if(zs_type=="conserve_y" or zs_type=="conserve_y_zs" or zs_type=="conserve_y_parallel") out=shrinkMatrix.adjoint()*out;
	outnorm=out.norm();
	if(print){
		cout<<"final version"<<endl;
		cout<<outnorm<<endl;
		for(int i=0;i<out.size();i++)
			if(norm(out(i))>-1) cout<<abs(out(i))/outnorm<<" "<<arg(out(i))/M_PI<<" "<<(bitset<NBITS>)lnd_states[i]<<endl;
	}
	return out/outnorm;
	
//	Eigen::FullPivHouseholderQR<Eigen::MatrixXcd> qr(Amatrix);

}

//computes elements of the many body wavefunction for a set of positions zs, which are on an integer grid
complex<double> NewChanger::get_wf(const vector< vector<int> > &zs){
	complex<double> out=1.;
	double x,y;
	int ix, iy;
	complex<double> temp,temp2;
    
	//vandermonde piece
	int vandermonde_exponent=invNu;
	if(type=="CFL") vandermonde_exponent-=2;
	complex<double> z;
	for( int i=0;i<Ne;i++){
		for( int j=i+1;j<Ne;j++){
			x=(zs[i][0]-zs[j][0])/(1.*NPhi);
			y=(zs[i][1]-zs[j][1])/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			//temp=jies_weierstrass(x,y);
			out*=pow(temp,vandermonde_exponent);
		}
	}
	
	//com p complex<double> 
	int COM[2]={0,0};
	for( int i=0;i<Ne;i++){
		COM[0]+=zs[i][0];
		COM[1]+=zs[i][1];
	}

	for( int i=0;i<invNu;i++){
		if(type=="CFL"){
			x=(COM[0]-dsum[0]/invNu)/(1.*NPhi)-mb_zeros[i][0];
			y=(COM[1]-dsum[1]/invNu)/(1.*NPhi)-mb_zeros[i][1];
		}
		else{
			x=COM[0]/(1.*NPhi)-mb_zeros[i][0];
			y=COM[1]/(1.*NPhi)-mb_zeros[i][1];
		
		}
		z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
		out*=temp;
	}
	double mb_zeros_sumx=0;
	if(type=="laughlin"){
	}
//	out*=polar(1.,-M_PI*COM[0]*mb_zeros[0][1]/(1.*Ne));
	for(int i=0;i<invNu;i++) mb_zeros_sumx+=mb_zeros[i][0];
//	out*=polar(1.,-M_PI*COM[1]*mb_zeros_sumx/(1.*NPhi));
	
	//determinant p complex<double> 
	Eigen::FullPivLU<Eigen::MatrixXcd> detSolver;
	if(type=="CFL"){
		complex<double> product;
		Eigen::MatrixXcd M(Ne,Ne);
		for(int i=0;i<Ne;i++){
			for(int j=0;j<Ne;j++){
				product=1;
				for(int k=0;k<Ne;k++){
					if(k==i) continue;
					ix=zs[i][0]-zs[k][0]-invNu*cfl_ds[j][0];
					iy=zs[i][1]-zs[k][1]-invNu*cfl_ds[j][1];
                    x=(ix+dsum[0]/(1.*Ne))/(1.*NPhi);
                    y=(iy+dsum[1]/(1.*Ne))/(1.*NPhi);
					//z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
					temp=modded_lattice_z(ix,iy);
					product*=temp;
				}
				//this p complex<double>  is only valid on a square torus!
				M(i,j)=product*polar(1., M_PI*(zs[i][1]*cfl_ds[j][0] - zs[i][0]*cfl_ds[j][1])/(1.*NPhi) );
			}
		}
//		cout<<M<<endl;
		detSolver.compute(M);
		temp=detSolver.determinant(); 
//		cout<<temp<<endl;
	  	out*=temp;
	}
        
	return out;
} 
void NewChanger::reset_ds(vector< vector<int> > ds, double ddbarx_t, double ddbary_t){
	ddbarx=ddbarx_t;
	ddbary=ddbary_t;

	//CFL ds
	cfl_ds=ds;
	int old_dsum=dsum[1];
	//vector< vector<int> > (Ne,vector<int>(2));
	dsum=vector<int>(2,0);
	for(int i=0;i<Ne;i++){
		dsum[0]+=cfl_ds[i][0]*invNu;
		dsum[1]+=cfl_ds[i][1]*invNu;
	}
	if(old_dsum%NPhi!=dsum[1]%NPhi){
		cout<<"you reset the ds into a different charge sector!"<<endl;
		exit(0);
	}
//	cout<<"dsum: "<<dsum[0]<<" "<<dsum[1]<<endl;
//	cout<<endl;

	complex<double> temp;
	double x,y;
	shifted_ztable=vector< vector< complex<double> > > (NPhi,vector< complex<double> >(NPhi,0));

	for(int ix=0;ix<NPhi;ix++){
		x=(ix+dsum[0]/(1.*Ne)+ddbarx)/(1.*NPhi);
		for(int iy=0;iy<NPhi;iy++){
			y=(iy+dsum[1]/(1.*Ne)+ddbary)/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			shifted_ztable[ix][iy]=temp;
		}
	}


}
inline double NewChanger::det_helper(int z1, int z2, int d, double dbarp){ return z1-z2-d*invNu+dbarp;}

complex<double> NewChanger::modded_lattice_z(int x, int y){// why need this function? because want to use z_function table, which is given in first BZ.
	int modx=supermod(x,NPhi);
	int mody=supermod(y,NPhi);
//	complex<double> out=lattice_z_(&NPhi,&modx,&mody,&L1,&L2,&one);
	complex<double> out=shifted_ztable[modx][mody];
	int j=(modx-x)/NPhi, k=(mody-y)/NPhi;
//	out*=omega[supermod(-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k,2*NPhi)];
	out*=polar( 1., (-(mody+dsum[1]/(1.*Ne)+ddbary)*j+(modx+dsum[0]/(1.*Ne)+ddbarx)*k)*M_PI/(1.*NPhi));
//	out*=polar(1.,-M_PI/(1.*NPhi)*(mody*j-modx*k));
//	cout<<polar(1.,-M_PI/(1.*NPhi)*(y*j-x*k))<<endl;
	if(j%2 || k%2) return -out;
	else return out;
}

void NewChanger::setup_mbl_zs(){

	bool symmetry_contract=true;
	bool found,found2;
	int seed=0;
	MTRand ran;
	ran.seed(seed);

	if(zs_type=="lines"){	
		//make the set of all possible positions, for mb side
		int tempbit;
		vector<unsigned int>::iterator it;
		for(int i=0;i<pow(2,NPhi);i++){
			if(count_bits(i)==Ne){
				tempbit=i;
				found=false;
				for(int x=0;x<NPhi;x++){
					tempbit=cycle_bits(tempbit,NPhi);
					it=find(mb_states.begin(),mb_states.end(),tempbit);
					if(it != mb_states.end()){
						found=true;
						break;
					} 
				}
				if(!found || !symmetry_contract) mb_states.push_back(i);
			} 
		}
	//	for(int i=0;i<mb_states.size();i++) cout<<(bitset<6>)mb_states[i]<<endl;
		ystart=0; ystep=1;
		copies=NPhi; //amount of redundancy 
	//	for(int i=0;i<mb_states.size();i++) cout<<(bitset<NBITS>)mb_states[i]<<endl;
		n_mb=mb_states.size()*copies;

		vector<int> bits;
		mb_zs=vector< vector< vector<int> > >(n_mb, vector< vector<int> > (Ne, vector<int>(2)));
		for(int i=0;i<(signed)mb_states.size();i++){
			bits=bitset_to_pos(mb_states[i],NPhi);	
		
			for(int y=0;y<copies;y++){
				for(int x=0;x<Ne;x++){
					mb_zs[i+y*mb_states.size()][x][0]=bits[x];
					mb_zs[i+y*mb_states.size()][x][1]=ystart+y*ystep;
				}
			}
		}
	}else if(zs_type=="random" or zs_type=="parallel"){
		n_mb=lnd_states.size();
		mb_zs=vector< vector< vector<int> > >(n_mb, vector< vector<int> > (Ne, vector<int>(2)));
	
		vector<int> tempr(2);
		found=true;
		for(int i=0;i<n_mb;i++){
			while(found){ 
				found=false;
				found2=true;
				//attempt to generate new set of zs
				for(int x=0;x<Ne;x++){
					//loop to avoid duplicates
					while(found2){
						found2=false;
						mb_zs[i][x][0]=ran.randInt(NPhi-1);
						mb_zs[i][x][1]=ran.randInt(NPhi-1);
						for(int y=0;y<x;y++){
							if(mb_zs[i][x]==mb_zs[i][y]){
								found2=true;
								break;
							}
						}
					}
					found2=true;
				}
				//check to make sure we haven't made the same set of zs twice
				for(int j=0;j<i;j++){
					if(mb_zs[i]==mb_zs[j]){
						found=true;
						break;
					}
				}
			}
			found=true;
		}
					
	}else if(zs_type=="conserve_y" or zs_type=="conserve_y_zs" or zs_type=="conserve_y_parallel"){
		//for this type, we randomly generate only n_lnd/Ne positions. Then we shift the y components of these by 1 Ne times
		//the resulting thing *should* have a symmetry upon shifting all positions in the y direction
		
		//the landau basis have some states which transform into themselves under T_y, so here we make slightly more states than the landau basis
		makeShrinker(supermod(dsum[0]/invNu+Ne,NPhi));

		int extra_sites=4;
		int n_sweeps=n_lnd/Ne;
		if (n_sweeps*Ne < n_lnd) n_sweeps++;		
		n_mb=n_sweeps*Ne*extra_sites;
		
		mb_zs=vector< vector< vector<int> > >(n_mb, vector< vector<int> > (Ne, vector<int>(2)));
		
		vector<int> tempr(2);
		for(int i=0;i<n_mb;i+=Ne){
			while(found){ 
				found=false;
				found2=true;
				//attempt to generate new set of zs
				for(int x=0;x<Ne;x++){
					//loop to avoid duplicates
					while(found2){
						found2=false;
						mb_zs[i][x][0]=ran.randInt(NPhi-1);
						mb_zs[i][x][1]=ran.randInt(NPhi-1);
						for(int y=0;y<x;y++){
							if(mb_zs[i][x]==mb_zs[i][y]){
								//cout<<"match found"<<endl;
								found2=true;
								break;
							}
						}
					}
					found2=true;
				}
				//check to make sure we haven't made the same set of zs twice
				for(int j=0;j<i;j++){
					if(mb_zs[i]==mb_zs[j]){
						found=true;
						break;
					}
				}

				//we've found a good one, now make Ne copies of it
				for(int j=0;j<Ne;j++){
					for(int x=0;x<Ne;x++){
						mb_zs[i+j][x][1]=supermod(mb_zs[i][x][1]+j*(NPhi/Ne),NPhi);
						mb_zs[i+j][x][0]=mb_zs[i][x][0];
					}
					//check to make sure we didn't accidentally make smth like 0101 
					for(int k=i;k<j;k++){
						if(mb_zs[i+j]==mb_zs[k]){
							found=true;
							break;
						}
					}
				}
			}
			found=true;
			
		}
	}else{
		cout<<"invalid mbl zs type"<<endl;
		exit(0);
	}
//	cout<<"n_mb, n_lnd: "<<n_mb<<" "<<n_lnd<<endl;
//	for(int i=0;i<n_mb;i++){
//		for(int x=0;x<Ne;x++) cout<<"("<<mb_zs[i][x][0]<<","<<mb_zs[i][x][1]<<") ";
//		cout<<endl;
//	}	
}
complex<double> NewChanger::landau_basis(int ix, int iy, int index){
	complex<double> out=1., temp;
	double x,y,yzero=index/(1.*NPhi),xzerosum=0;
	int twoNPhi=2*NPhi;
	for(int i=0;i<NPhi;i++){
		x=ix/(1.*NPhi)-lnd_zeros[i];
		y=iy/(1.*NPhi)-yzero;
		z_function_(&x,&y,&L1,&L2,&one,&twoNPhi,&temp);
		out*=temp;
		xzerosum+=lnd_zeros[i];
	}
//	if(index%2) out*=-1.;
	out*=polar(1.,-M_PI*(yzero)*ix);
	out*=polar(1.,M_PI*(xzerosum/(1.*NPhi))*iy);
	
	return out;
}

void NewChanger::make_manybody_vector(){
	manybody_vector=Eigen::VectorXcd::Zero(n_mb);
	for(int i=0;i<n_mb;i++){
		manybody_vector(i)=get_wf(mb_zs[i]);
//		cout<<setprecision(10)<<mb_zs[i][0][0]<<" "<<mb_zs[i][0][1]<<" "<<mb_zs[i][1][0]<<" "<<mb_zs[i][1][1]<<" "<<manybody_vector(i)<<endl;
	}
	
	//normalize
	double norm=manybody_vector.norm();
	manybody_vector=manybody_vector/norm;
}

void NewChanger::make_Amatrix(){
	Amatrix=Eigen::MatrixXcd::Zero(n_mb,n_lnd);
	Eigen::MatrixXcd detMatrix(Ne,Ne);
	Eigen::FullPivLU<Eigen::MatrixXcd> LUsolver;
	vector<int> lnd_zs;
	double temp;

	if(zs_type=="conserve_y"){
		int xsum,count,sign,newstate,oldstate,newindex,jsize;
		vector<unsigned int>::iterator it;
		vector< vector<int> > j_copies(n_lnd), j_signs(n_lnd);
		vector<int> sites;

		for(int j=0;j<n_lnd;j++){
			//first, make a vector of all the states which are just this state, shifted by invNu
			//also keep track of any signs coming from permutations
			
			//cout<<(bitset<NBITS>)lnd_states[j]<<endl;
			newstate=lnd_states[j];
			j_copies[j].push_back(j);
			j_signs[j].push_back(1);
			for(int i=0;i<Ne;i++){
				temp=1;
				for(int k=0;k<invNu;k++){
					oldstate=newstate;
					newstate=0;
					sites=bitset_to_pos(oldstate,NPhi);
					for(int j=0;j<(signed)sites.size(); j++) newstate=newstate | 1<<((sites[j]+1)%NPhi);
					if (newstate<oldstate && Ne%2==0) temp*=-1;
				}
//				if(newstate==(signed)lnd_states[j]) break;
				it=find(lnd_states.begin(),lnd_states.end(),newstate);
				if(it!=lnd_states.end()){
					newindex=it-lnd_states.begin();
					if(invNu%2) temp*=-1.;
					j_copies[j].push_back(newindex);
					//cout<<newindex<<" "<<(bitset<NBITS>)newstate<<endl;
					j_signs[j].push_back(temp);
				}else{
					cout<<"state: "<<(bitset<NBITS>)newstate<<" not found"<<endl;
					exit(0);
				}
			}
			jsize=j_copies[j].size();
			//cout<<jsize<<endl<<endl;
		
			//now fill in the A matrix			
			lnd_zs=bitset_to_pos(lnd_states[j],NPhi);
			for(int i=0;i<n_mb;i+=Ne){
				xsum=0;
				//for the first row of each set of Ne rows, have to calculate the determinants explicitly
				for(int x=0;x<Ne;x++){
					xsum+=mb_zs[i][x][0];
					for(int y=0;y<Ne;y++){
						detMatrix(x,y)=lnd_table[lnd_zs[y]] [mb_zs[i][x][0]] [mb_zs[i][x][1]];
					}
				}										
				LUsolver.compute(detMatrix);
				Amatrix(i,j)=LUsolver.determinant()*pow(0.2,Ne*Ne);

				//for the other rows, just multiply the first row by a phase!
				for(int k=1;k<Ne;k++){
					count=0;
					for(int x=0;x<Ne;x++){
						if(mb_zs[i+k-1][x][1]>=NPhi-invNu) count+=mb_zs[i+k-1][x][0];
						//if(xsum%2==0 && mb_zs[i+k-1][x][1]>=NPhi-2*invNu) count++;
					}
					sign=1;
					if(count%2) sign=-1;
					//cout<<sign<<" "<<xsum*2/(1.*NPhi)<<" "<<j_signs[j][k%jsize]<<endl;
					Amatrix(i+k,j_copies[j][k%jsize])=(double)(j_signs[j][k%jsize]*sign)*polar(1.,-xsum*2.*M_PI/(1.*NPhi))*Amatrix(i+k-1,j_copies[j][(k-1)%jsize]);
				}
			}
		}
	}else if(zs_type=="conserve_y_parallel"){
		int newstate,oldstate,newindex,nthreads;
		vector<unsigned int>::iterator it;
		vector< vector<int> > j_copies(n_lnd), j_signs(n_lnd);
		vector<int> sites;

		//how many cores can we use
		nthreads=omp_get_max_threads();
		if(n_mb/Ne<nthreads) nthreads=n_mb/Ne;
		cout<<"number of threads used: "<<nthreads<<endl;
		omp_set_num_threads(nthreads);
	    vector<Eigen::PartialPivLU<Eigen::MatrixXcd> > LUsolverPar(nthreads);
	    vector<Eigen::MatrixXcd> detMatrixPar(nthreads,Eigen::MatrixXcd(Ne,Ne));
			
		for(int j=0;j<n_lnd;j++){
			//first, make a vector of all the states which are just this state, shifted by invNu
			//also keep track of any signs coming from permutations
			
			//cout<<(bitset<NBITS>)lnd_states[j]<<endl;
			newstate=lnd_states[j];
			j_copies[j].push_back(j);
			j_signs[j].push_back(1);
			for(int i=0;i<Ne;i++){
				temp=1;
				for(int k=0;k<invNu;k++){
					oldstate=newstate;
					newstate=0;
					sites=bitset_to_pos(oldstate,NPhi);
					for(int j=0;j<(signed)sites.size(); j++) newstate=newstate | 1<<((sites[j]+1)%NPhi);
					if (newstate<oldstate && Ne%2==0) temp*=-1;
				}
//				if(newstate==(signed)lnd_states[j]) break;
				it=find(lnd_states.begin(),lnd_states.end(),newstate);
				if(it!=lnd_states.end()){
					newindex=it-lnd_states.begin();
					if(invNu%2) temp*=-1.;
					j_copies[j].push_back(newindex);
					//cout<<newindex<<" "<<(bitset<NBITS>)newstate<<endl;
					j_signs[j].push_back(temp);
				}else{
					cout<<"state: "<<(bitset<NBITS>)newstate<<" not found"<<endl;
					exit(0);
				}
			}
		}

#pragma omp parallel for
		for(int j=0;j<n_lnd;j++){
	        int coren = omp_get_thread_num();
			//now fill in the A matrix			
			vector<int> lnd_zs=bitset_to_pos(lnd_states[j],NPhi);
			for(int i=0;i<n_mb;i+=Ne){
				//for the first row of each set of Ne rows, have to calculate the determinants explicitly
				for(int x=0;x<Ne;x++){
					for(int y=0;y<Ne;y++){
						detMatrixPar[coren](x,y)=lnd_table[lnd_zs[y]] [mb_zs[i][x][0]] [mb_zs[i][x][1]];
					}
				}										
				LUsolverPar[coren].compute(detMatrixPar[coren]);
				Amatrix(i,j)=LUsolverPar[coren].determinant()*pow(0.2,Ne*Ne);
			}
		}
		int count,xsum,jsize;
		for(int j=0;j<n_lnd;j++){
			jsize=j_copies[j].size();
			for(int i=0;i<n_mb;i+=Ne){
				//for the other rows, just multiply the first row by a phase!

				xsum=0;
				for(int x=0;x<Ne;x++) xsum+=mb_zs[i][x][0];
														
				for(int k=1;k<Ne;k++){
					count=0;
					for(int x=0;x<Ne;x++){
						if(mb_zs[i+k-1][x][1]>=NPhi-invNu) count+=mb_zs[i+k-1][x][0];
						//if(xsum%2==0 && mb_zs[i+k-1][x][1]>=NPhi-2*invNu) count++;
					}
					int sign=1;
					if(count%2) sign=-1;
					//cout<<sign<<" "<<xsum*2/(1.*NPhi)<<" "<<j_signs[j][k%jsize]<<endl;
					Amatrix(i+k,j_copies[j][k%jsize])=(double)(j_signs[j][k%jsize]*sign)*polar(1.,-xsum*2.*M_PI/(1.*NPhi))*Amatrix(i+k-1,j_copies[j][(k-1)%jsize]);
				}
			}
		}
				
	}else if (zs_type=="parallel"){
		int nthreads=omp_get_max_threads();
	    omp_set_num_threads(nthreads);
	    vector<Eigen::PartialPivLU<Eigen::MatrixXcd> > LUsolverPar(nthreads);
	    vector<Eigen::MatrixXcd> detMatrixPar(nthreads,Eigen::MatrixXcd(Ne,Ne));
#pragma omp parallel for			
		for(int j=0;j<n_lnd;j++){
	        int coren = omp_get_thread_num();
			vector<int> lnd_zsPar=bitset_to_pos(lnd_states[j],NPhi);

			for(int i=0;i<n_mb;i++){
				for(int x=0;x<Ne;x++){
					for(int y=0;y<Ne;y++){
	//					detMatrix(x,y)=landau_basis(mb_zs[i][x],lnd_zs[y]);
						detMatrixPar[coren](x,y)=lnd_table[lnd_zsPar[y]] [mb_zs[i][x][0]] [mb_zs[i][x][1]];
					}
				}
				LUsolverPar[coren].compute(detMatrixPar[coren]);
				Amatrix(i,j)=LUsolverPar[coren].determinant()*pow(0.2,Ne*Ne);
			}
		}
	
	}else{
	
		for(int j=0;j<n_lnd;j++){
			lnd_zs=bitset_to_pos(lnd_states[j],NPhi);
			for(int i=0;i<n_mb;i++){
				for(int x=0;x<Ne;x++){
					for(int y=0;y<Ne;y++){
	//					detMatrix(x,y)=landau_basis(mb_zs[i][x],lnd_zs[y]);
						detMatrix(x,y)=lnd_table[lnd_zs[y]] [mb_zs[i][x][0]] [mb_zs[i][x][1]];
					}
				}
				LUsolver.compute(detMatrix);
				Amatrix(i,j)=LUsolver.determinant()*pow(0.2,Ne*Ne);
			}
		}
	}
	if(zs_type=="conserve_y" or zs_type=="conserve_y_zs" or zs_type=="conserve_y_parallel") Amatrix=Amatrix*shrinkMatrix.adjoint();		
}

void NewChanger::make_landau_table(){
	lnd_table=vector< vector< vector< complex<double> > > > (NPhi,vector< vector< complex<double> > >(NPhi, vector< complex<double> > (NPhi,0.) ) );
	for(int k=0;k<NPhi;k++){
		for(int x=0;x<NPhi;x++){
			for(int y=0;y<NPhi;y++){
				lnd_table[k][x][y]=landau_basis(x,y,k);
			}
		}
	}
}
				
void NewChanger::make_landau_symmetry_x(){
	Tx_lnd=Eigen::MatrixXcd::Zero(n_lnd,n_lnd);
	int xcount;
	double sign=1;
	if(Ne%2==0) sign=-1;
	for(int i=0;i<n_lnd;i++){
		xcount=0;
		for(int x=0;x<NPhi;x++){
			if(lnd_states[i] & 1<<x) xcount+=x;
		}
		Tx_lnd(i,i)=sign*polar(1.,-2*M_PI*xcount/(1.*NPhi));
	}
}
void NewChanger::make_landau_symmetry_y(){
	Eigen::MatrixXcd tempOut;
	Ty_lnd=Eigen::MatrixXcd::Zero(n_lnd,n_lnd);
	vector<int> sites;
	int newstate=0,newindex,oldstate;
	vector<unsigned int>::const_iterator it;
	double temp;

	for(int i=0;i<n_lnd;i++){

		temp=1;
		//translation matrix
		newstate=lnd_states[i];
		for(int k=0;k<invNu;k++){
			oldstate=newstate;
			newstate=0;
			sites=bitset_to_pos(oldstate,NPhi);
			for(int j=0;j<(signed)sites.size(); j++) newstate=newstate | 1<<((sites[j]+1)%NPhi);
			if (newstate<oldstate && Ne%2==0) temp*=-1;
		}
		it=find(lnd_states.begin(),lnd_states.end(),newstate);
		if(it!=lnd_states.end()){
			newindex=it-lnd_states.begin();
			if(invNu%2) temp*=-1.;
			Ty_lnd(i,newindex)=temp;//*polar(1.,-2*M_PI*(ystart+y)/(1.*NPhi));
		}else{
			cout<<"state: "<<(bitset<NBITS>)newstate<<" not found"<<endl;
			exit(0);
		}
		
		
	}
}

void NewChanger::make_manybody_symmetry_x(){
	Tx_mb=Eigen::MatrixXcd::Zero(n_mb,n_mb);

	vector<int> sites;
	int newstate=0,newindex;
	vector<unsigned int>::const_iterator it;
	complex<double> temp;
	double sign;

	for(int i=0;i<(signed)mb_states.size();i++){

		//translation matrix
		newstate=0;
		sites=bitset_to_pos(mb_states[i],NPhi);
		for(int j=0;j<(signed)sites.size(); j++) newstate=newstate | 1<<((sites[j]+1)%NPhi);
		it=find(mb_states.begin(),mb_states.end(),newstate);
		if(it!=mb_states.end()){
			newindex=it-mb_states.begin();
			temp=-1.;
			for(int y=0;y<copies;y++){
				sign=1;
				if(mb_states[newindex] <= mb_states[i]){
					if (type=="laughlin" ) sign*=-1.;
					if (type== "CFL" && Ne%2==0) sign*=-1;
					if (y*ystep % 2) sign*=-1; //this one is from the gauge term below
					//if (type=="laughlin" && y%2) sign*=-1;
				}
				Tx_mb(i+y*mb_states.size(),newindex+y*mb_states.size())=temp*sign*polar(1.,-M_PI*y*ystep/(1.*invNu));
			}
		}else{
			cout<<"state: "<<(bitset<NBITS>)newstate<<" not found in many body symmetry x"<<endl;
			exit(0);
		}
		
		
	}
}
void NewChanger::make_manybody_symmetry_y(){
	Ty_mb=Eigen::MatrixXcd::Zero(n_mb,n_mb);
	int sum_x;
	vector<int> xs;
	double sign=-1.;
	for(int i=0;i<(signed)mb_states.size();i++){
		xs=bitset_to_pos(mb_states[i],NPhi);
		sum_x=accumulate(xs.begin(),xs.end(),0);
		for(int j=0;j<copies;j++){
			if(invNu%2==0) sign=1;
			else sign=-1;
			if( j*ystep+invNu>=NPhi){
				if(sum_x%2 ) sign*=-1;
				//if(Ne%2) sign*=-1;
			}
			Ty_mb(i+j*mb_states.size(),i+( (j+invNu/ystep)%copies )*mb_states.size() )=sign*polar(1.,M_PI*sum_x/(1.*Ne));
		}
	}
}
int NewChanger::get_dsum(int dir){ return dsum[dir]; }
void NewChanger::test(){
	double x=0.5;
	double y=0.5;
	complex<double> temp;
	z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
	cout<<L1<<" "<<L2<<endl;
	cout<<setprecision(10)<<"testing duncan's z function"<<endl;
	cout<<temp<<endl;
}	

//construct denstiy operator rho(kx,ky)
Eigen::SparseMatrix< complex<double> > NewChanger::density_operator(int my, int mx){
	//this operator doesn't conserve charge, and so it maps between different bases
	//this p complex<double>  constructs the basis to map to, it works just like make_states
	//the new basis has its charge INCREASED by kx
	int verbose=0,xcharge;
	vector<unsigned int> new_states;
	if(mx==0) new_states=lnd_states;
	else{
		for(int i=0;i<intpow(2,NPhi);i++){
			xcharge=0;
			for(int x=0;x<NPhi;x++)
				if(i & 1<<x) xcharge+=x;
			if(count_bits(i)==Ne && xcharge%NPhi==supermod(lnd_charge+mx,NPhi) )
				new_states.push_back(i);
		}
	}
	int newNStates=new_states.size();

	if(verbose>0){	
		for(int i=0;i<(signed)lnd_states.size();i++) cout<<(bitset<16>)lnd_states[i]<<endl;
		cout<<endl;
		for(int i=0;i<(signed)new_states.size();i++) cout<<(bitset<16>)new_states[i]<<endl;
	}
	Eigen::SparseMatrix< complex<double> > rho(newNStates,n_lnd);
	vector<Eigen::Triplet< complex<double> > > triplets;

	double Lx=sqrt(2*M_PI*NPhi), Ly=Lx;
	//for every state, find all the states which you get from moving one electron kx, and add those matrix elements to rho
	 complex<double>  prefactor;
	double sign,kx,ky;
	if(mx>NPhi/2) kx=2*M_PI/Lx*(mx-NPhi);
	else kx=2*M_PI/Lx*mx;
	if(my>NPhi/2) ky=2*M_PI/Ly*(my-NPhi);
	else ky=2*M_PI/Ly*my;
	vector<unsigned int>::iterator it;
	unsigned int newstate;
	prefactor=polar(exp(-0.25*(pow(kx,2)+pow(ky,2))), 0.5*kx*ky);
	for(int i1=0; i1<n_lnd; i1++){
		if(verbose>1) cout<<(bitset<6>)lnd_states[i1]<<endl;
		for(int x=0; x<NPhi; x++){
			try{
				newstate=move_bit(lnd_states[i1],NPhi,x,mx);
			}catch(int e){
				continue;
			}
			if(verbose>1) cout<<x<<" "<<(bitset<6>)newstate<<endl;
			it=find(new_states.begin(),new_states.end(),newstate);
			if(it==new_states.end()) continue;
			//a minus sign from normal ordering if there's an even number of electrons
			if( newstate<lnd_states[i1] && Ne%2==0 ) sign=-1.;
			else sign=1.;
			//another minus sign for every electron this electron hops over
			if (x+mx>NPhi){
				for(int y=0;y<(x+mx)%NPhi;y++)
					if(lnd_states[i1] & 1<<y) sign*=-1;
				for(int y=x+1;y<NPhi;y++)
					if(lnd_states[i1] & 1<<y) sign*=-1;
			}else{
				for(int y=x+1;y<x+mx;y++) 
					if(lnd_states[i1] & 1<<y) sign*=-1;
			}
			triplets.push_back(Eigen::Triplet< complex<double> >( it-new_states.begin(),i1,sign*prefactor*polar(1.,ky*2*M_PI/Lx*x) ) );
		}
	}
	rho.setFromTriplets(triplets.begin(),triplets.end());
	if(verbose>1) cout<<rho<<endl;
	return rho;				
}

//make shrinking matrix
void NewChanger::makeShrinker(int nx){
	vector<Eigen::Triplet< complex<double> > > triplets;
	vector<Eigen::Triplet< complex<double> > > temptrips;
	vector<int> found_states;

	int temp,index;
	vector<unsigned int>::iterator it;
	int col=0, phase,sign;

//	for(int i=0;i<nStates;i++) cout<<(bitset<12>)states[i]<<endl;
//	cout<<endl;
	for(int i=0;i<n_lnd;i++){
		if(find(found_states.begin(),found_states.end(),i)!=found_states.end()) continue;
		
		it=lnd_states.begin()+i;
		temp=lnd_states[i];
		phase=0;
		temptrips.clear();
		sign=1;
		while(true){
			index=it-lnd_states.begin();
			temptrips.push_back( Eigen::Triplet< complex<double> >( col,index,polar(1.*sign,phase*2*M_PI/(1.*Ne)*nx ) ) );
			found_states.push_back(index);
//			cout<<(bitset<12>)states[index]<<" "<<phase*nx%Ne<<" "<<sign<<endl;
			temp=cycle_M(temp,NPhi,invNu,sign);
			if(temp==(signed)lnd_states[i]) break;
			it=lower_bound(lnd_states.begin(),lnd_states.end(),temp);
			phase++;
			
		}
		//some cases only exist in certain momentum sectors (eg 0101) 
		if( (Ne%2!=0 && nx%(Ne/temptrips.size()) ) || (Ne%2==0 && (nx-Ne/2)%(Ne/temptrips.size()) ) ) {
//			cout<<"cancelled"<<endl;
			 continue;
		}

		for(unsigned int j=0;j<temptrips.size();j++) 
			temptrips[j]=Eigen::Triplet< complex<double> >(temptrips[j].col(), temptrips[j].row(), temptrips[j].value()/sqrt(temptrips.size()));
		triplets.insert(triplets.end(),temptrips.begin(),temptrips.end());
//		cout<<endl;
		col++;
	}
	shrinkMatrix=Eigen::SparseMatrix< complex<double> >(n_lnd,col);
//	cout<<"triplet size "<<triplets.size()<<endl;
//	for(unsigned int j=0;j<triplets.size();j++)
//		cout<<triplets[j].row()<<" "<<triplets[j].col()<<" "<<triplets[j].value()<<" "<<(bitset<10>)states[triplets[j].row()]<<endl;
	shrinkMatrix.setFromTriplets(triplets.begin(),triplets.end());
	shrinkMatrix=shrinkMatrix.adjoint();
}

