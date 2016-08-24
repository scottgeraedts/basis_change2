#include "NewChanger.h"
#include <Eigen/Sparse>
#define NBITS 10

NewChanger::NewChanger(){

}
NewChanger::NewChanger(int NPhi_t, int Ne_t, int manybody_COM_t, string type_t, vector< vector<int> >ds, bool contract_t, double ddbarx_t=0., double ddbary_t=0.):
	NPhi(NPhi_t),Ne(Ne_t),manybody_COM(manybody_COM_t),type(type_t),symmetry_contract(contract_t){

	invNu=NPhi/Ne;
	double Lx=sqrt(2*M_PI*NPhi);
	double Ly=Lx;
	L1=complex<double>(Lx/sqrt(2.),0);
	L2=complex<double>(0,Ly/sqrt(2.));
	
	ddbarx=ddbarx_t;
	ddbary=ddbary_t;
	
	zero=0; //for calls to duncan's functions
	one=1; 
	
	//make the set of all possible positions, for mb side
	int tempbit;
	bool found;
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
	copies=NPhi/ystep/Ne; //amount of redundancy 
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
Eigen::VectorXcd NewChanger::run(bool print, bool compute_A=true){
	make_manybody_vector();
	if(compute_A){
		make_landau_table();
		make_Amatrix();
	}
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
	outnorm=out.norm();
	if(print){
		cout<<"final version"<<endl;
		cout<<outnorm<<endl;
		if(!manybody_vector.isApprox(Amatrix*out,1e-9)) cout<<"bad solution!"<<endl;
		for(int i=0;i<n_lnd;i++)
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
	
	//com part
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
	out*=polar(1.,-M_PI*COM[0]*mb_zeros[0][1]/(1.*Ne));
	for(int i=0;i<invNu;i++) mb_zeros_sumx+=mb_zeros[i][0];
	out*=polar(1.,-M_PI*COM[1]*mb_zeros_sumx/(1.*NPhi));
	
	//determinant part
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
				//this part is only valid on a square torus!
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
void NewChanger::reset_ds(vector< vector<int> > ds, double ddbarx_t=0., double ddbary_t=0.){
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
	if(old_dsum!=dsum[1]){
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
			cout<<"state: "<<(bitset<NBITS>)newstate<<" not found"<<endl;
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
int main(){

//	int NPhi,Ne,manybody_COM;
//	string type;
//	ifstream params("params");
//	params>>NPhi;
//	params>>Ne;
//	params>>manybody_COM;
//	params>>type;
//	params.close();

//	vector<vector<int> > cfl_ds=vector<vector<int> >(Ne,vector<int>(2,0));
//	ifstream kfile("kfile");
//	cout<<"ks"<<endl;
//	for(int i=0;i<Ne;i++){
//		kfile>>cfl_ds[i][0]>>cfl_ds[i][1];		
//		cout<<cfl_ds[i][0]<<" "<<cfl_ds[i][1]<<endl;
//	}
//	kfile.close();
//	NewChanger control(NPhi,Ne,manybody_COM,type,cfl_ds,true);
//	Eigen::VectorXcd vec0=control.run(false);

//	void ph_overlap2(int Ne, int NPhi, string type, vector< vector<int> > cfl_ds, const NewChanger &control,  const Eigen::VectorXcd &vec0);
//	ph_overlap2(Ne, NPhi, type, cfl_ds, control, vec0);

//	void batch_overlap();
//	batch_overlap();

	void orthogonality();
	orthogonality();
}
void ph_overlap(int Ne, int NPhi, string type, vector< vector<int> > cfl_ds, const NewChanger &control, const Eigen::VectorXcd &vec0){
	vector<vector<int> > new_cfl_ds=vector<vector<int> >(Ne,vector<int>(2,0));
	for(int i=0;i<Ne;i++){
		new_cfl_ds[i][0]=-cfl_ds[i][0];
		new_cfl_ds[i][1]=-cfl_ds[i][1]+1;
	}
//	new_cfl_ds=cfl_ds;
//	new_cfl_ds[Ne-1][0]=1-new_cfl_ds[Ne-1][0];
//	new_cfl_ds[Ne-1][1]=1-new_cfl_ds[Ne-1][1];
	NewChanger test1(NPhi,Ne,0,type,new_cfl_ds,true);	
	Eigen::VectorXcd vec1=test1.run(false);

	///***MAKE PH SYMMETRY TRANSLATION MATRIX***///
	vector<Eigen::Triplet<complex<double> > > ph_triplets;
	vector<unsigned int>::iterator it;
	int partner,sign,xcharge,j;
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if(! (control.lnd_states[i] & 1<<x)) partner=partner | 1<<x;

		it=find(test1.lnd_states.begin(),test1.lnd_states.end(),partner);
		if(it != test1.lnd_states.end()){
			j=it-test1.lnd_states.begin();
			xcharge=0;
			for(int x=0;x<NPhi;x++) 
				if(control.lnd_states[i] & 1<<x) xcharge+=x;
			
			if(xcharge%2) sign=-1;
			else sign=1;
			ph_triplets.push_back(Eigen::Triplet<complex<double> >(i,j,sign) ); 		
		}else{
			cout<<(bitset<NBITS>)partner<<" not found!"<<endl;
		}
	}
	Eigen::SparseMatrix<complex<double> > ph_sym(control.lnd_states.size(),test1.lnd_states.size());
	ph_sym.setFromTriplets(ph_triplets.begin(),ph_triplets.end());
	
	cout<<"actual dot products: "<<endl;
//	cout<<vec0<<endl<<endl;
//	cout<<ph_sym*vec1<<endl<<endl;
	complex<double> overlap=vec0.dot( (ph_sym*vec1).conjugate());
	cout<<overlap<<endl;
	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;
	
	cout<<"absolute value dot products: "<<endl;
	int n_lnd=control.lnd_states.size();
	Eigen::VectorXcd abs0(n_lnd),abs1(n_lnd);
	for(int i=0;i<(signed)n_lnd;i++){
		abs0(i)=abs(vec0(i));
		abs1(i)=abs(vec1(i));
	}
//	cout<<abs0<<endl<<endl;
//	cout<<ph_sym*abs1<<endl<<endl;
	cout<<abs0.dot(ph_sym*abs1)<<endl;
	
		
}
double ph_overlap2(int Ne, int NPhi, string type, vector< vector<int> > cfl_ds, const NewChanger &control, const Eigen::VectorXcd &vec0){
	vector<vector<int> > new_cfl_ds=vector<vector<int> >(Ne,vector<int>(2,0));
	for(int i=0;i<Ne;i++){
		new_cfl_ds[i][0]=-cfl_ds[i][0];
		new_cfl_ds[i][1]=-cfl_ds[i][1]+1;
	}
	NewChanger test1(NPhi,Ne,0,type,new_cfl_ds,true);	
	Eigen::VectorXcd vec1=test1.run(false);

	///***MAKE PH SYMMETRY TRANSLATION MATRIX***///
	vector<Eigen::Triplet<complex<double> > > ph_triplets;
	vector<unsigned int>::const_iterator it;
	int partner,sign,xcharge,j;
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if(! (control.lnd_states[i] & 1<<x)) partner=partner | 1<<x;

		it=find(test1.lnd_states.begin(),test1.lnd_states.end(),partner);
		if(it != test1.lnd_states.end()){
			j=it-test1.lnd_states.begin();
			xcharge=0;
			for(int x=0;x<NPhi;x++) 
				if(control.lnd_states[i] & 1<<x) xcharge+=x;
			
			if(xcharge%2) sign=-1;
			else sign=1;
			ph_triplets.push_back(Eigen::Triplet<complex<double> >(i,j,sign) ); 		
		}else{
			cout<<(bitset<NBITS>)partner<<" not found in ph!"<<endl;
		}
	}
	Eigen::SparseMatrix<complex<double> > ph_sym(control.lnd_states.size(),control.lnd_states.size());
	ph_sym.setFromTriplets(ph_triplets.begin(),ph_triplets.end());

	///***MAKE INVERSION MATRIX***///
	ph_triplets.clear();
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if( control.lnd_states[i] & 1<<x) partner=partner | 1<<(NPhi-1-x);

		it=find(test1.lnd_states.begin(),test1.lnd_states.end(),partner);
		if(it != test1.lnd_states.end()){
			j=it-test1.lnd_states.begin();
			ph_triplets.push_back(Eigen::Triplet<complex<double> >(i,j,1) ); 		
		}else{
			cout<<(bitset<NBITS>)partner<<" not found in inv!"<<endl;
		}
	}
	Eigen::SparseMatrix<complex<double> > inv_sym(control.lnd_states.size(),control.lnd_states.size());
	inv_sym.setFromTriplets(ph_triplets.begin(),ph_triplets.end());	

	for(int i=0;i<Ne;i++){
		new_cfl_ds[i][0]=cfl_ds[i][0];
		new_cfl_ds[i][1]=cfl_ds[i][1]+1;
	}
	NewChanger test2(NPhi,Ne,0,type,new_cfl_ds,true);	
	Eigen::VectorXcd vec2=test2.run(false);

	///***MAKE COM INVERSION MATRIX***///
	ph_triplets.clear();
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if( control.lnd_states[i] & 1<<x) partner=partner | 1<<supermod(x+1,NPhi);

		it=find(test2.lnd_states.begin(),test2.lnd_states.end(),partner);
		if(it != test2.lnd_states.end()){
			j=it-test2.lnd_states.begin();
			ph_triplets.push_back(Eigen::Triplet<complex<double> >(i,j,1) ); 		
		}else{
			cout<<(bitset<NBITS>)partner<<" not found in com!"<<endl;
		}
	}
	Eigen::SparseMatrix<complex<double> > com_sym(control.lnd_states.size(),control.lnd_states.size());
	com_sym.setFromTriplets(ph_triplets.begin(),ph_triplets.end());	

	for(int i=0;i<Ne;i++){
		new_cfl_ds[i][0]=-cfl_ds[i][0];
		new_cfl_ds[i][1]=cfl_ds[i][1];
	}
	NewChanger test3(NPhi,Ne,0,type,new_cfl_ds,true);	
	Eigen::VectorXcd vec3=test3.run(false);


	complex<double> overlap;
	overlap=(vec0).dot( inv_sym*vec1);
	cout<<"R2:"<<endl;
	cout<<overlap<<endl;
	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;
	
	overlap=(vec0).dot( ph_sym*vec1.conjugate());
	cout<<"PH:"<<endl;
	cout<<overlap<<endl;
	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

	overlap=(vec0).dot( inv_sym*ph_sym.transpose()*vec0.conjugate());
	cout<<"PH R2:"<<endl;
	cout<<overlap<<endl;
	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

	overlap=(vec0).dot( com_sym*vec2);
	cout<<"COM:"<<endl;
	cout<<overlap<<endl;
	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

	overlap=(vec0).dot( vec3.conjugate());
	cout<<"Rx:"<<endl;
	cout<<overlap<<endl;
	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

	return abs(overlap);	
}

void batch_overlap(){

	int Ne=8;
	int NPhi=16;
	vector< vector<int> > cfl_ds( Ne, vector<int>(2));
	ifstream kfile("batch_ds");
	vector<int>  Dbar(2);
	double Dvar,ddbarx,ddbary;
	Eigen::VectorXcd vec0;

	ifstream control_ks("kfile");
	for(int x=0;x<Ne;x++){
		control_ks>>cfl_ds[x][0]>>cfl_ds[x][1];
//			cout<<cfl_ds[x][0]<<" "<<cfl_ds[x][1]<<" ";
	}
		
	while(true){
//		for(int x=0;x<Ne;x++){
//			kfile>>cfl_ds[x][0]>>cfl_ds[x][1];
////			cout<<cfl_ds[x][0]<<" "<<cfl_ds[x][1]<<" ";
//		}
		kfile>>ddbarx>>ddbary;
//		cout<<endl;
		if (kfile.eof()){
//			cout<<"eof reached"<<endl;
			break;
		}
		
		Dbar=vector<int>(2,0);
		for(int x=0; x<Ne; x++){
			Dbar[0]+=cfl_ds[x][0];
			Dbar[1]+=cfl_ds[x][1];
		}
		Dvar=0.;
//		for(int x=0; x<Ne; x++){
//			Dvar+=sqrt( pow(cfl_ds[x][0]-(1.*Dbar[0])/(1.*Ne),2)+pow(cfl_ds[x][1]-(1.*Dbar[1])/(1.*Ne),2));
//		}
//		Dvar/=sqrt(Ne);
		Dvar+=sqrt( pow(ddbarx/(2.*Ne),2)+pow(ddbary/(2.*Ne),2));
		
		NewChanger control(NPhi,Ne,0,"CFL",cfl_ds,true,ddbarx/(1.*Ne),ddbary/(1.*Ne));
		vec0=control.run(false);
		cout<<Dvar<<" "<<ph_overlap2(Ne,NPhi,"CFL",cfl_ds,control,vec0)<<endl;
	}
	kfile.close();
}

void orthogonality(){
	int Ne=8;
	int NPhi=16;
	bool first=true;
	vector< vector<int> > cfl_ds( Ne, vector<int>(2));
	ifstream kfile("batch_ks");
	vector<Eigen::VectorXcd> vecs;
	NewChanger control;
	cout<<"(dx,dy) values "<<endl;
	while(true){
		for(int x=0;x<Ne;x++){
			kfile>>cfl_ds[x][0]>>cfl_ds[x][1];
		}
		if (kfile.eof()){
			break;
		}
		for(int x=0;x<Ne;x++){
			cout<<"("<<cfl_ds[x][0]<<" "<<cfl_ds[x][1]<<"), ";
		}
		cout<<endl;
		
		
		if(first) control=NewChanger(NPhi,Ne,0,"CFL",cfl_ds,true);
		else control.reset_ds(cfl_ds);
		vecs.push_back(control.run(true,first));
		first=false;
	}
	cout<<endl;
	kfile.close();
	
	int nvecs=vecs.size();
	Eigen::MatrixXcd overlaps(nvecs,nvecs);
	for(int i=0; i<nvecs; i++){
		for(int j=0; j<nvecs; j++){
			overlaps(i,j)=vecs[i].dot(vecs[j]);
		}
	} 
	cout<<"overlaps between the different states considered"<<endl;
	cout<<overlaps<<endl<<endl;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(overlaps);
	cout<<"eigenvalues of the overlap matrix  (to reveal rank)"<<endl;
	cout<<es.eigenvalues()<<endl<<endl;
	
	///***MAKE PH SYMMETRY TRANSLATION MATRIX***///
	vector<Eigen::Triplet<complex<double> > > ph_triplets;
	vector<unsigned int>::const_iterator it;
	int partner,sign,xcharge,j;
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if(! (control.lnd_states[i] & 1<<x)) partner=partner | 1<<supermod(NPhi-1-x,NPhi);

		it=find(control.lnd_states.begin(),control.lnd_states.end(),partner);
		if(it != control.lnd_states.end()){
			j=it-control.lnd_states.begin();
			xcharge=0;
			for(int x=0;x<NPhi;x++) 
				if(control.lnd_states[i] & 1<<x) xcharge+=x;
			
			if(xcharge%2) sign=-1;
			else sign=1;
			ph_triplets.push_back(Eigen::Triplet<complex<double> >(i,j,sign) ); 		
		}else{
			cout<<(bitset<NBITS>)partner<<" not found in ph!"<<endl;
		}
	}
	Eigen::SparseMatrix<complex<double> > ph_sym(control.lnd_states.size(),control.lnd_states.size());
	ph_sym.setFromTriplets(ph_triplets.begin(),ph_triplets.end());
	complex<double> overlap;
	cout<<"particle-hole (plus R2) symmetry of the various states"<<endl;
	for (int i=0; i<nvecs; i++){
		overlap=vecs[i].dot(ph_sym*vecs[i].conjugate());
		cout<<abs(overlap)<<" "<<endl;
	}
	cout<<endl;
	
	//***COMPARE TO ED STATE***///
	cout<<"overlaps with ED ground states"<<endl;
	double r,phi;
	int conf;
	stringstream filename;
	filename<<"eigen"<<supermod(control.get_dsum(0)/2,Ne)<<"_"<<supermod(control.get_dsum(1)/2,Ne)+Ne;
	ifstream eigen(filename.str().c_str());
	cout<<filename.str()<<endl;
	Eigen::VectorXcd EDstate(control.lnd_states.size());
	for(int i=0; i<(signed)control.lnd_states.size(); i++){
		eigen>>r>>phi>>conf;
		EDstate(i)=polar(r,phi*M_PI);
	}
	for(int i=0;i<nvecs;i++){
		overlap=vecs[i].dot(EDstate);
		cout<<abs(overlap)<<endl;
	}
}

