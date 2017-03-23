#include "callchanger.h"
#ifdef USE_CLUSTER
ArpackError::ErrorCode ArpackError::code = NO_ERRORS;
#endif


int main(){

//	int NPhi, Ne;
//	vector< vector<int> > cfl_ds;
//	ifstream kfile("batch_ds");

//	vector<NewChanger> wfs;	
//	string zs_type;

//	kfile>>NPhi>>Ne;
//	kfile>>zs_type;	
//	cfl_ds=vector< vector<int> >( Ne, vector<int>(2));
//	time_t t;
//	map<string,double> params;
//	while(true){
//		t=time(NULL);
//		for(int x=0;x<Ne;x++){
//			kfile>>cfl_ds[x][0]>>cfl_ds[x][1];
//		}
//		if (kfile.eof()){
////			cout<<"eof reached"<<endl;
//			break;
//		}
//		for(int x=0;x<Ne;x++){
//			cout<<"("<<cfl_ds[x][0]<<","<<cfl_ds[x][1]<<") ";
//		}
//		cout<<endl;
//				
//		NewChanger control(NPhi,Ne,0,"CFL",cfl_ds,params,"random");
//		wfs.push_back(control);
//	}
//	
//	int nzs=10;
//	Eigen::MatrixXcd z_coeffs(wfs.size(),nzs);
//	for(int nvec=0;nvec<(signed)wfs.size();nvec++){
//		for(int i=0;i<nzs;i++){
//			z_coeffs(nvec,i)=wfs[nvec].get_wf(wfs[0].mb_zs[i])/1e10;
//		}
//	}
//	cout<<z_coeffs<<endl;
//	Eigen::JacobiSVD<Eigen::MatrixXcd> es(z_coeffs);
//	cout<<es.singularValues()<<endl;
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
//	NewChanger control(NPhi,Ne,manybody_COM,type,cfl_ds);
//	Eigen::VectorXcd vec0=control.run(true);

//	ph_overlap2(Ne, NPhi, type, cfl_ds, control, vec0);

//	CFL_berry();
//	batch_overlap();

//	orthogonality();

	energy_variance();
}

//given a model wf, computes another model wf at an appropriate momentum so that the action of PH can be tested
//currently deprecated
/*
void ph_overlap(int Ne, int NPhi, string type, vector< vector<int> > cfl_ds, const NewChanger &control, const Eigen::VectorXcd &vec0){
	vector<vector<int> > new_cfl_ds=vector<vector<int> >(Ne,vector<int>(2,0));
	for(int i=0;i<Ne;i++){
		new_cfl_ds[i][0]=-cfl_ds[i][0];
		new_cfl_ds[i][1]=-cfl_ds[i][1]+1;
	}
//	new_cfl_ds=cfl_ds;
//	new_cfl_ds[Ne-1][0]=1-new_cfl_ds[Ne-1][0];
//	new_cfl_ds[Ne-1][1]=1-new_cfl_ds[Ne-1][1];
	NewChanger test1(NPhi,Ne,0,type,new_cfl_ds);	
	Eigen::VectorXcd vec1=test1.run(false);

	///MAKE PH SYMMETRY TRANSLATION MATRIX///
	vector<Eigen::Triplet<complex<double> > > ph_triplets;
	vector<unsigned int>::iterator it;
	int partner,sign,xcharge,j;
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if(! (control.lnd_states[i] & one<<x)) partner=partner | one<<x;

		it=find(test1.lnd_states.begin(),test1.lnd_states.end(),partner);
		if(it != test1.lnd_states.end()){
			j=it-test1.lnd_states.begin();
			xcharge=0;
			for(int x=0;x<NPhi;x++) 
				if(control.lnd_states[i] & one<<x) xcharge+=x;
			
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
//	cout<<ph_sym*vecone<<endl<<endl;
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
//	cout<<ph_sym*absone<<endl<<endl;
	cout<<abs0.dot(ph_sym*abs1)<<endl;
	
		
}
*/

//given a model wavefunction, finds the overlap with the action of PH + inversion
double ph_overlap2(int Ne, int NPhi, string type, vector< vector<int> > cfl_ds, const NewChanger &control, const Eigen::VectorXcd &vec0){
	state_int one=1;
	vector<vector<int> > new_cfl_ds=vector<vector<int> >(Ne,vector<int>(2,0));
	for(int i=0;i<Ne;i++){
		new_cfl_ds[i][0]=-cfl_ds[i][0];
		new_cfl_ds[i][1]=-cfl_ds[i][1]+1;
	}
	map<string,double> params;
	params["theta"]=60;
	NewChanger test1(NPhi,Ne,0,type,new_cfl_ds,params);	
//	Eigen::VectorXcd vec1=test1.run(false);

	///***MAKE PH SYMMETRY TRANSLATION MATRIX***///
	vector<Eigen::Triplet<complex<double> > > ph_triplets;
	vector<state_int>::iterator it;
	int partner,sign,xcharge,j;
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if(! (control.lnd_states[i] & one<<x)) partner=partner | one<<x;

		it=find(test1.lnd_states.begin(),test1.lnd_states.end(),partner);
		if(it != test1.lnd_states.end()){
			j=it-test1.lnd_states.begin();
			xcharge=0;
			for(int x=0;x<NPhi;x++) 
				if(control.lnd_states[i] & one<<x) xcharge+=x;
			
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
			if( control.lnd_states[i] & one<<x) partner=partner | one<<(NPhi-1-x);

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

//	for(int i=0;i<Ne;i++){
//		new_cfl_ds[i][0]=cfl_ds[i][0];
//		new_cfl_ds[i][1]=cfl_ds[i][1]+1;
//	}
//	NewChanger test2(NPhi,Ne,0,type,new_cfl_ds);	
//	Eigen::VectorXcd vec2=test2.run(false);

//	///***MAKE COM INVERSION MATRIX***///
//	ph_triplets.clear();
//	for(int i=0;i<(signed)control.lnd_states.size();i++){
//		partner=0;
//		for(int x=0;x<NPhi;x++)
//			if( control.lnd_states[i] & one<<x) partner=partner | one<<supermod(x+1,NPhi);

//		it=find(test2.lnd_states.begin(),test2.lnd_states.end(),partner);
//		if(it != test2.lnd_states.end()){
//			j=it-test2.lnd_states.begin();
//			ph_triplets.push_back(Eigen::Triplet<complex<double> >(i,j,1) ); 		
//		}else{
//			cout<<(bitset<NBITS>)partner<<" not found in com!"<<endl;
//		}
//	}
//	Eigen::SparseMatrix<complex<double> > com_sym(control.lnd_states.size(),control.lnd_states.size());
//	com_sym.setFromTriplets(ph_triplets.begin(),ph_triplets.end());	

//	for(int i=0;i<Ne;i++){
//		new_cfl_ds[i][0]=-cfl_ds[i][0];
//		new_cfl_ds[i][1]=cfl_ds[i][1];
//	}
//	NewChanger test3(NPhi,Ne,0,type,new_cfl_ds,true);	
//	Eigen::VectorXcd vec3=test3.run(false);


	complex<double> overlap;
//	overlap=(vec0).dot( inv_sym*vec1);
//	cout<<"R2:"<<endl;
//	cout<<overlap<<endl;
//	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;
//	
//	overlap=(vec0).dot( ph_sym*vec1.conjugate());
//	cout<<"PH:"<<endl;
//	cout<<overlap<<endl;
//	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

	overlap=(vec0).dot( inv_sym*ph_sym.transpose()*vec0.conjugate());
//	cout<<"PH R2:"<<endl;
//	cout<<overlap<<endl;
//	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

//	overlap=(vec0).dot( com_sym*vec2);
//	cout<<"COM:"<<endl;
//	cout<<overlap<<endl;
//	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

//	overlap=(vec0).dot( vec3.conjugate());
//	cout<<"Rx:"<<endl;
//	cout<<overlap<<endl;
//	cout<<abs(overlap)<<" "<<arg(overlap)<<endl;

	return abs(overlap);	
}

//calls ph overlap many times
void batch_overlap(){

	int Ne=4;
	int NPhi=16;
	vector< vector<int> > cfl_ds( Ne, vector<int>(2));
	ifstream kfile("batch_ds");
	vector<int>  Dbar(2);
	double Dvar,ddbarx,ddbary;
	Eigen::VectorXcd vec0;
	map<string,double> params;

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
		
		
		params["ddbarx"]=ddbarx/(1.*Ne);
		params["ddbary"]=ddbary/(1.*Ne);
		NewChanger control(NPhi,Ne,0,"CFL",cfl_ds,params);
		vec0=control.run(false);
		cout<<Dvar<<" "<<ph_overlap2(Ne,NPhi,"CFL",cfl_ds,control,vec0)<<endl;
	}
	kfile.close();
}

//finds the overlaps between different model wf, good for checking rank
void orthogonality(){
	state_int one=1;

	int Ne;
	int NPhi;
	string crap;
	ifstream params("params");
	params>>Ne>>crap;
	params>>NPhi>>crap;
	bool first=true;
	vector< vector<int> > cfl_ds( Ne, vector<int>(2));
	ifstream kfile("batch_ks");
	vector<Eigen::VectorXcd> vecs, rvecs;
	NewChanger control;
	cout<<"(dx,dy) values "<<endl;
	map<string,double> params_map;
	ifstream genfile("genparams");
	genfile>>params_map["alpha"];
	genfile>>params_map["theta"];
	genfile.close();
	
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
		 	
		
		if(first) control=NewChanger(NPhi,Ne,0,"CFL",cfl_ds,params_map,"conserve_y");
		else control.reset_ds(cfl_ds);
		vecs.push_back(control.run(false,first));
		rvecs.push_back(control.manybody_vector);
		first=false;
	}
	cout<<endl;
	kfile.close();
	
	int nvecs=vecs.size();
	//take symmetric and antisymmetric versions
//	vector<Eigen::VectorXcd> tempvecs(2);
//	for(int i=1;i<nvecs;i+=2){
//		tempvecs[0]=vecs[i];
//		tempvecs[1]=vecs[i+1];
//		vecs[i]=(tempvecs[0]+tempvecs[1]);
//		vecs[i]=vecs[i]/vecs[i].norm();
//		vecs[i+1]=(tempvecs[0]-tempvecs[1]);
//		vecs[i+1]=vecs[i+1]/vecs[i+1].norm();
//	}
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
	vector<state_int>::const_iterator it;
	int partner,sign,xcharge,j;
	for(int i=0;i<(signed)control.lnd_states.size();i++){
		partner=0;
		for(int x=0;x<NPhi;x++)
			if(! (control.lnd_states[i] & one<<x)) partner=partner | one<<supermod(NPhi-1-x,NPhi);

		it=find(control.lnd_states.begin(),control.lnd_states.end(),partner);
		if(it != control.lnd_states.end()){
			j=it-control.lnd_states.begin();
			xcharge=0;
			for(int x=0;x<NPhi;x++) 
				if(control.lnd_states[i] & one<<x) xcharge+=x;
			
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
	cout<<"particle-hole (plus R2) symmetry of the model states"<<endl;
	for (int i=0; i<nvecs; i++){
		overlap=vecs[i].dot(ph_sym*vecs[i].conjugate());
		cout<<abs(overlap)<<" "<<endl;
	}
	cout<<endl;
	
	//***COMPARE TO ED STATE***///
	//relation between these and NewChanger parameters dx,dy: 
	//kx = dy+8
	//ky = dx+8
	int kx,ky;
	if(Ne%2) kx=supermod(control.get_dsum(1)/2,NPhi);
	else kx=supermod(control.get_dsum(1)/2+Ne,NPhi);
	ky=supermod(control.get_dsum(0)/2+Ne,NPhi);
	
	TorusSolver< complex<double> > T(kx);
	T.store_sparse=true;
	T.make_Hnn();

	T.makeShrinker(ky);
	//cout<<H1.shrinkMatrix.rows()<<" "<<H1.shrinkMatrix.cols()<<" "<<H1.EigenDense.rows()<<endl;
//		Hnn=Eigen::MatrixXcd(T.shrinkMatrix * T.EigenSparse * T.shrinkMatrix.adjoint());
	Eigen::MatrixXcd Hnn;
	Eigen::SparseMatrix<complex <double> > tempMat=(T.shrinkMatrix * T.EigenSparse * T.shrinkMatrix.adjoint());

	vector<Eigen::VectorXcd> smallvecs(nvecs);
	for(int i=0; i<nvecs; i++) smallvecs[i]=T.shrinkMatrix*vecs[i];
	
	ofstream vout("vecs");
	for(int i=0;i<(signed)smallvecs[0].size();i++){
//	for(int i=0;i<(signed)control.lnd_states.size();i++){
		for(int j=0;j<nvecs;j++) vout<<abs(smallvecs[j](i))<<" "<<arg(smallvecs[j](i))<<" ";
		vout<<(bitset<NBITS>)control.lnd_states[i];
		vout<<endl;
	}
	vout.close();

	int n_exacts=tempMat.rows();
	if(Ne<=5) n_exacts=5;
	vector< Eigen::VectorXcd> EDout(n_exacts);
	if(Ne>10){
#ifdef USE_CLUSTER
		MatrixWithProduct3 mat2(tempMat.rows());
		mat2.set(tempMat);
		mat2.eigenvalues(n_exacts);
		for(int i=0;i<n_exacts;i++) EDout[i]=T.shrinkMatrix.adjoint()*Std_To_Eigen(mat2.eigvecs[i]);
		//ED_E=mat2.eigvals[0];
#else
		cout<<"you need to use the cluster to use arpack"<<endl;
#endif
	}else{
		Hnn=Eigen::MatrixXcd(tempMat);
		es.compute(Hnn);
		for(int i=0;i<n_exacts;i++) EDout[i]=T.shrinkMatrix.adjoint()*es.eigenvectors().col(i);
		//ED_E=es.eigenvalues()(0);
	}
//	for(int i=0;i<n_exacts;i++) EDout[i]=T.shrinkMatrix.adjoint()*EDout[i];

	Eigen::MatrixXcd coeffs(n_exacts,nvecs);
	cout<<"model states in terms of exact states"<<endl;
	vector<double> norms(nvecs,0);
	for(int i=0;i<n_exacts;i++){
		for(int j=0;j<nvecs;j++){
			coeffs(i,j)=EDout[i].dot(vecs[j]);
			cout<<coeffs(i,j)<<" ";
			norms[j]+=abs(coeffs(i,j));
		}
		cout<<endl;
	}
	for(int j=0;j<nvecs;j++){
		if(1-norms[j]>0.02) cout<<"vector "<<j<<" not well described by the computed exact states"<<endl;
	}
}


//computes the berry phase going on a path few many different wf
void CFL_berry(){
	//make a TorusSolver object that has x momentum 1
	int NPhi,nks,Ne;

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<complex<double>,-1,-1> > es;
	Eigen::MatrixXcd Hnn;

	ifstream kfile("batch_ds");
	kfile>>NPhi>>Ne;
	vector<vector<int> > temp(Ne,vector<int>(2));
	vector< vector< vector<int> > > ds,ds2;
	while(true){
		for(int i=0;i<Ne;i++) kfile>>temp[i][0]>>temp[i][1];
		if(kfile.eof()) break;
		ds.push_back(temp);
		for(int i=0;i<Ne;i++) temp[i][1]++;
		ds2.push_back(temp);
	}
	vector< vector<int> >charges(ds.size(),vector<int>(2,0));
	for(int i=0;i<(signed)ds.size();i++){
		for(int j=0;j<Ne;j++){
			for(int l=0;l<2;l++) charges[i][l]+=ds[i][j][l];
		}
		charges[i][0]=supermod(charges[i][0],NPhi);
		charges[i][1]=supermod(charges[i][1],NPhi);
	}
	kfile.close();		
	nks=ds.size();
	vector<Eigen::VectorXcd> vecs1(nks), vecs2(nks);
	complex<double> tempz,overlap;
	Eigen::MatrixXcd tempmat=Eigen::MatrixXcd::Zero(2,2),product;
	product=Eigen::MatrixXcd::Identity(2,2);
	vector<Eigen::SparseMatrix< complex<double> > > rhos11(nks), rhos12(nks), rhos21(nks), rhos22(nks);
	map<string,double> params;

	int dkx,dky; 
	for(unsigned int k=0;k<(unsigned)nks;k++){
		NewChanger H1(NPhi, Ne, 0, "CFL", ds[k], params), H2(NPhi, Ne, 0, "CFL", ds2[k], params);

		dkx=charges[(k+1)%ds.size()][0]-charges[k][0];
		dky=charges[(k+1)%ds.size()][1]-charges[k][1];
		if(dky<0) dky+=NPhi;
		if(dkx<0) dkx+=NPhi;
		for(int i=0;i<Ne;i++) cout<<"("<<ds[k][i][0]<<","<<ds[k][i][1]<<")" ;
		cout<<endl;
		for(int i=0;i<Ne;i++) cout<<"("<<ds2[k][i][0]<<","<<ds2[k][i][1]<<")" ;
		cout<<endl;
		cout<<charges[k][0]<<" "<<charges[k][1]<<" "<<dkx<<" "<<dky<<endl;

		rhos11[k]=H1.density_operator(dkx,dky);
//		cout<<rhos11[k]<<endl;
		rhos12[k]=H1.density_operator(dkx,(dky+NPhi/2)%NPhi);
		rhos22[k]=H2.density_operator(dkx,dky);
		rhos21[k]=H2.density_operator(dkx,(dky+NPhi/2)%NPhi);

		vecs1[k]=H1.run(false);
		vecs2[k]=H2.run(false);

		//***COMPARE TO ED STATE***///
		cout<<"overlaps with ED ground states"<<endl;
		double r,phi;
		int conf;
		stringstream filename;
		filename<<"eigen"<<k;
		ifstream eigen(filename.str().c_str());
		cout<<filename.str()<<endl;
		Eigen::VectorXcd EDstate(H1.lnd_states.size());
		for(int i=0; i<(signed)H1.lnd_states.size(); i++){
			eigen>>r>>phi>>conf;
			EDstate(i)=polar(r,phi*M_PI);
		}
		overlap=vecs1[k].dot(EDstate);
		cout<<abs(overlap)<<endl;

	}
	for(unsigned int k=0;k<ds.size();k++){
//		tempz=vecs[k].dot(vecs[(k+1)%ds.size()].conjugate());
		cout<<rhos11[k].rows()<<" "<<rhos11[k].cols()<<" "<<vecs1[k].size()<<" "<<vecs1[(k+1)%ds.size()].size()<<endl;
		tempz=vecs1[(k+1)%ds.size()].adjoint()*rhos11[k]*vecs1[k];
		cout<<abs(tempz)<<" "<<arg(tempz)<<endl;
		tempmat(0,0)=tempz;
		
		tempz=vecs2[(k+1)%ds.size()].adjoint()*rhos22[k]*vecs2[k];
		cout<<abs(tempz)<<" "<<arg(tempz)<<endl;
		tempmat(1,1)=tempz;
		
		tempz=vecs2[(k+1)%ds.size()].adjoint()*rhos12[k]*vecs1[k];
		cout<<abs(tempz)<<" "<<arg(tempz)<<endl;
		tempmat(0,1)=tempz;
		
		tempz=vecs1[(k+1)%ds.size()].adjoint()*rhos21[k]*vecs2[k];
		cout<<abs(tempz)<<" "<<arg(tempz)<<endl;
		tempmat(1,0)=tempz;
		cout<<endl;	
		product*=tempmat;
	}
	cout<<product<<endl;
}
//does a bunch of tests to asses how 'good' a model wf is
//these include: overlap with ED state, <H^2>-<H>^2, PH symmetry
void energy_variance(){

	int Ne,NPhi,kx,ky;
	vector< vector<int> > cfl_ds;
	ifstream kfile("batch_ds");
	vector<int>  Dbar(2), states;
	double Dvar;
//	Eigen::VectorXcd vec0;
//	TorusSolver<complex <double> > T(0);
	Eigen::MatrixXcd Hnn;
	double E2,E,ED_E;
	complex<double> temp_mult;
	Eigen::VectorXcd ev1,EDout;
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<complex<double>,-1,-1> > es;
	
	double self_energy;
	bool have_self_energy=false;	
	string zs_type;
	map<string,double> params;
	params["alpha"]=1;
	params["theta"]=90;
		
	kfile>>NPhi>>Ne;
	kfile>>zs_type;	
	cfl_ds=vector< vector<int> >( Ne, vector<int>(2));
	vector<Eigen::VectorXcd> vec0;	
	vector<NewChanger> controls;

	int nvec=0;
	time_t t;
	while(true){
		t=time(NULL);
		for(int x=0;x<Ne;x++){
			kfile>>cfl_ds[x][0]>>cfl_ds[x][1];
		}
		if (kfile.eof()){
//			cout<<"eof reached"<<endl;
			break;
		}
		for(int x=0;x<Ne;x++){
			cout<<"("<<cfl_ds[x][0]<<","<<cfl_ds[x][1]<<") ";
		}
		cout<<endl;
		
		Dbar=vector<int>(2,0);
		for(int x=0; x<Ne; x++){
			Dbar[0]+=cfl_ds[x][0];
			Dbar[1]+=cfl_ds[x][1];
		}
		Dvar=0.;
		for(int x=0; x<Ne; x++){
			Dvar+=sqrt( pow(cfl_ds[x][0]-(1.*Dbar[0])/(1.*Ne),2)+pow(cfl_ds[x][1]-(1.*Dbar[1])/(1.*Ne),2));
		}
		Dvar/=sqrt(Ne);
		cout<<"Dvar: "<<Dvar<<endl;
		
		NewChanger control(NPhi,Ne,0,"HLR",cfl_ds,params,zs_type);
		cout<<"constructed"<<endl;
		vec0.push_back(control.run(false));
		cout<<"ran"<<endl;
		if(zs_type=="lines") control.symmetry_checks();

		controls.push_back(control);
//	}
//	Eigen::VectorXcd tempvec0=vec0[0]+vec0[1];
//	Eigen::VectorXcd tempvec1=vec0[0]-vec0[1];
//	vec0[0]=tempvec0/tempvec0.norm();
//	vec0[1]=tempvec1/tempvec1.norm();
//	for(int nvec;nvec<(signed)vec0.size();nvec++){
		//relation between these and NewChanger parameters dx,dy: 
		//kx = dy+8
		//ky = dx+8
		if(Ne%2) kx=supermod(Dbar[1],NPhi);
		else if(NPhi/Ne==2) kx=supermod(Dbar[1]+Ne,NPhi);
		else kx=supermod(Dbar[1]+2*Ne,NPhi);
		ky=supermod(Dbar[0]+Ne,NPhi);

		TorusSolver<complex<double> > T(kx);
		T.store_sparse=true;
		T.make_Hnn();

		T.makeShrinker(ky);
		cout<<T.shrinkMatrix.rows()<<" "<<T.shrinkMatrix.cols()<<" "<<T.EigenSparse.rows()<<endl;
//		Hnn=Eigen::MatrixXcd(T.shrinkMatrix * T.EigenSparse * T.shrinkMatrix.adjoint());
		Eigen::SparseMatrix<complex <double> > tempMat=(T.shrinkMatrix * T.EigenSparse * T.shrinkMatrix.adjoint());

		if(Ne>10){
#ifdef USE_CLUSTER
			MatrixWithProduct3 mat2(tempMat.rows());
			mat2.set(tempMat);
			mat2.eigenvalues(5);
			EDout=Std_To_Eigen(mat2.eigvecs[0]);
			ED_E=mat2.eigvals[0];
#else
			cout<<"you need to use the cluster to use arpack"<<endl;
#endif
		}else{
//			Hnn=Eigen::MatrixXcd(T.EigenDense);
			Hnn=Eigen::MatrixXcd(tempMat);
			es.compute(Hnn);
			EDout=es.eigenvectors().col(0);
			ED_E=es.eigenvalues()(0);
		}	
		ev1=T.shrinkMatrix.adjoint()*EDout;
//		ev1=EDout;
		states=T.get_states();

		cout<<"Ed state"<<endl;
		Eigen::VectorXd absED(states.size()), absWF(states.size());
//		for(int i=0;i<(signed)states.size();i++){
//			cout<<abs(ev1(i))<<" "<<arg(ev1(i))/M_PI<<" "<<(bitset<NBITS>)states[i]<<" ";
//			cout<<abs(vec0[nvec](i))<<" "<<arg(vec0[nvec](i))<<endl;
//			absWF(i)=abs(vec0[nvec](i));
//			absED(i)=abs(ev1(i));
//			
//		}
//		cout<<endl;
//		
		cout<<"absolute overlaps"<<absED.dot(absWF)<<endl;
		
		if(!have_self_energy) self_energy=T.self_energy();
		cout<<"ED energy: "<<ED_E/(1.*Ne)+self_energy<<endl;
		//cout<<EDout.size()<<" "<<ev1.size()<<" "<<vec0.size()<<" "<<control.lnd_states.size()<<endl;
//		if(NPhi/Ne==2) cout<<"ED PH symmetry: "<<ph_overlap2(Ne,NPhi,"CFL",cfl_ds,controls[nvec],ev1)<<endl;

		temp_mult=vec0[nvec].adjoint()*T.EigenSparse*T.EigenSparse*vec0[nvec];
		E2=real(temp_mult);
		temp_mult=vec0[nvec].adjoint()*T.EigenSparse*vec0[nvec];
		E=real(temp_mult);
		temp_mult=vec0[nvec].dot(ev1);
		cout<<"overlap with ED state: "<<abs(temp_mult)<<endl;
		cout<<"energy, variance: "<<E/(1.*Ne)+self_energy<<" "<<(E2-E*E)/(1.*Ne*Ne)<<endl;
		if(NPhi/Ne==2) cout<<"PH symmetry: "<<ph_overlap2(Ne,NPhi,"CFL",cfl_ds,controls[nvec],vec0[nvec])<<endl;
		cout<<endl;	
		cout<<"time elapsed: "<<difftime(time(NULL),t)<<endl;
		cout<<endl;
		nvec++;
	}
	kfile.close();
}

