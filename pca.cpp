//Designed by Junbo ZHAO
//2013.10.29

#include <iostream>
#include <fstream>
#include "matrix.h"
#include "eigen.h"

using namespace std;

bool fexist(const string filename){
	ifstream fid(filename.c_str());
	return fid;
}

bool PCAcalc(const string datafile,const string resfile,int pcanum){
	if(fexist(resfile))
		return true;
	if(!fexist(datafile))
		return false;
	Matrix Fea,Fea_reduce;
	Fea.readfile(datafile.c_str());
	vector<vector<double>> mdata(1,Fea.mean(1));
	Matrix mu(1,Fea.Getcols());
	mu.setdata(mdata);
	if(Fea.Getrows()>Fea.Getcols()){
		double size = Fea.Getrows();
		Matrix St = Fea.trans()*Fea/size-mu.trans()*mu;
		Matrix V;
		vector<double> ei;
		if(!St.eig(ei,&V)){
			exit(-1);
		}
		V = V.trunc(pcanum,2);		
		Fea_reduce = Fea * V;
		Fea_reduce.writefile(resfile.c_str());
	}
	else{								//For High-Dimensional conditions
		Matrix *t = new Matrix;
		(*t) = ones(Fea.Getrows(),1)*mu;
		Matrix St = (Fea-(*t))*(Fea.trans()-(*t).trans());	//U*U' -> U'*U
		delete t;
		Matrix Q;
		vector<double> ei;
		if(!St.eig(ei,&Q)){				
			exit(-1);
		}
		ei.pop_back();
		Matrix I = diag(ei);
		I = I.m_sqrt();
		I = I.invert();
		Q = Q.trunc(Q.Getcols()-1,2);
		Matrix V = Fea.trans()*Q*I;
		V = V.trunc(pcanum,2);
		Fea_reduce = Fea * V;
		Fea_reduce.writefile(resfile.c_str());
	}
	return true;

}

void main(){
	if(PCAcalc("pcadata.txt","pca_reduce.txt",40))
		cout<<"PCA Finished"<<endl;
	else 
		cout<<"PCA Failed"<<endl;
}
