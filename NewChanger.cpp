#include "NewChanger.h"
#define NBITS 10

NewChanger::NewChanger(){

	ifstream params("params");
	params>>NPhi;
	params>>Ne;
	params>>manybody_COM;
	params>>type;
	invNu=NPhi/Ne;
	double Lx=sqrt(2*M_PI*NPhi);
	double Ly=Lx;
	L1=complex<double>(Lx/sqrt(2.),0);
	L2=complex<double>(0,Ly/sqrt(2.));
	
	zero=0; //for calls to duncan's functions
	one=1; 
	
	//make the set of all possible positions, for mb side
	for(int i=0;i<pow(2,NPhi);i++){
		if(count_bits(i)==Ne) mb_states.push_back(i);
	}
	ystart=0; ystep=1;
	copies=NPhi/ystep; //amount of redundancy 
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
	ifstream kfile("kfile");
	cfl_ds=vector< vector<int> > (Ne,vector<int>(2));
	cout<<"ks"<<endl;
	for(int i=0;i<Ne;i++){
		kfile>>cfl_ds[i][0]>>cfl_ds[i][1];		
		cout<<cfl_ds[i][0]<<" "<<cfl_ds[i][1]<<endl;
		cfl_ds[i][0]=supermod(cfl_ds[i][0],NPhi);
		cfl_ds[i][1]=supermod(cfl_ds[i][1],NPhi);
	}
	dsum=vector<int>(2,0);
	for(int i=0;i<Ne;i++){
		dsum[0]+=cfl_ds[i][0]*invNu;
		dsum[1]+=cfl_ds[i][1]*invNu;
	}
	dsum[0]=supermod(dsum[0],NPhi);
	dsum[1]=supermod(dsum[1],NPhi);
	cout<<"dsum: "<<dsum[0]<<" "<<dsum[1]<<endl;
	cout<<endl;

	complex<double> temp;
	double x,y;
	shifted_ztable=vector< vector< complex<double> > > (NPhi,vector< complex<double> >(NPhi,0));

	for(int ix=0;ix<NPhi;ix++){
		x=(ix+dsum[0]/(1.*Ne))/(1.*NPhi);
		for(int iy=0;iy<NPhi;iy++){
			y=(iy+dsum[1]/(1.*Ne))/(1.*NPhi);
			z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
			shifted_ztable[ix][iy]=temp;
		}
	}

	
	//make the set of all possible positions, on the lnd side
	int xcharge;
	lnd_charge=Ne*manybody_COM;
	if(Ne%2==0) lnd_charge+=NPhi/2; //for some reason you dont get the charge sector you expect for even Ne
	if(type=="CFL") lnd_charge+=dsum[1]/invNu;
	lnd_charge=lnd_charge%NPhi;
	cout<<"charge: "<<lnd_charge<<endl;
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
	
	make_manybody_vector();
	make_landau_table();
	make_Amatrix();
	make_manybody_symmetry_x();
	make_manybody_symmetry_y();
	make_landau_symmetry_x();
	make_landau_symmetry_y();
}
void NewChanger::symmetry_checks(){
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
void NewChanger::run(){
//*****SOLUTION***///
	Eigen::VectorXcd out,tempB=Eigen::VectorXcd(mb_states.size());;
	Eigen::MatrixXcd tempA;
	double outnorm;
	for(int i=0;i<copies;i++){
		tempA=Eigen::MatrixXcd(mb_states.size(),n_lnd);
		for(int x=0;x<mb_states.size();x++){
			tempB(x)=manybody_vector(x+i*mb_states.size());
			for(int y=0;y<n_lnd;y++){
				tempA(x,y)=Amatrix(x+i*mb_states.size(),y);
			}
		}
//		tempA=Amatrix.middleRows(mb_states.size()*i,mb_states.size()*(i+1));
//		tempB=manybody_vector.segment(mb_states.size()*i,mb_states.size()*(i+1));
		cout<<tempA.rows()<<" "<<tempA.cols()<<" "<<tempB.size()<<endl;
//		out=Amatrix.middleRows(mb_states.size()*i,mb_states.size()*(i+1)).colPivHouseholderQr().solve(manybody_vector.segment(mb_states.size()*i,mb_states.size()*(i-1)));	
		out=tempA.colPivHouseholderQr().solve(tempB);	

		outnorm=out.norm();
		Eigen::VectorXcd outSym=Ty_lnd*out;
		cout<<outnorm<<endl;
		if(!manybody_vector.isApprox(Amatrix*out,1e-9)) cout<<"bad solution!"<<endl;
		for(int i=0;i<n_lnd;i++)
			if(norm(out(i))>1e-16) cout<<out(i)/outnorm<<" "<<(bitset<NBITS>)lnd_states[i]<<" "<<outSym(i)/outnorm<<endl;
	}
	out=Amatrix.colPivHouseholderQr().solve(manybody_vector);	
	outnorm=out.norm();
	Eigen::VectorXcd outSym=Ty_lnd*out;
	cout<<"final version"<<endl;
	cout<<outnorm<<endl;
	if(!manybody_vector.isApprox(Amatrix*out,1e-9)) cout<<"bad solution!"<<endl;
	for(int i=0;i<n_lnd;i++)
		if(norm(out(i))>1e-16) cout<<abs(out(i)/outnorm)<<" "<<arg(out(i))/M_PI<<" "<<(bitset<NBITS>)lnd_states[i]<<endl;

	
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

//	if(type=="CFL"){
//		COM[0]-=dsum[0]/invNu;//this is correct since dsum lives on the same lattice as the zs
//		COM[1]-=dsum[1]/invNu;
//	}
	//cout<<COM[0]<<" "<<COM[1]<<endl;
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
	out*=polar(1.,-M_PI*COM[0]*mb_zeros[0][1]/(1.*Ne));
	double mb_zeros_sumx=0;
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
                    ix=zs[i][0]-zs[k][0]-cfl_ds[j][0];
                    iy=zs[i][1]-zs[k][1]-cfl_ds[j][1];
                    x=(ix+dsum[0]/(1.*Ne))/(1.*NPhi); y=(iy+dsum[1]/(1.*Ne))/(1.*NPhi);
//					temp=modded_lattice_z(ix,iy);
					z_function_(&x,&y,&L1,&L2,&zero,&NPhi,&temp);
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

inline double NewChanger::det_helper(int z1, int z2, int d, double dbarp){ return z1-z2-d*invNu+dbarp;}

//call's duncan's lattice_z function, if the arguments x or y are outside of the range (0,NPhi) it shifts them into that range
//and multiplies by the appropriate phase
//only works for a square torus
//on a square torus, the phase is always +- 1? 
complex<double> NewChanger::modded_lattice_z(int x, int y){// why need this function? because want to use z_function table, which is given in first BZ.
	int modx=supermod(x,NPhi);
	int mody=supermod(y,NPhi);
//	complex<double> out=lattice_z_(&NPhi,&modx,&mody,&L1,&L2,&one);
	complex<double> out=shifted_ztable[modx][mody];
	int j=(modx-x)/NPhi, k=(mody-y)/NPhi;
//	out*=omega[supermod(-(mody+dbar_parameter[1])*j+(modx+dbar_parameter[0])*k,2*NPhi)];
	out*=polar( 1., (-(mody+dsum[1]/(1.*Ne))*j+(modx+dsum[0]/(1.*Ne))*k)*M_PI/(1.*NPhi));
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
			Amatrix(i,j)=LUsolver.determinant()*pow(0.2,NPhi);
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
void NewChanger::test(){

}	
int main(){
	NewChanger object;
	object.symmetry_checks();
	object.run();	
}
