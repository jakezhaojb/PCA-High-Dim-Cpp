#include "matrix.h"

void Matrix::readfile(std::string fname){
	ifstream fin;
	fin.open(fname.c_str());
	if(!fin){
		cerr<<"No Input File!";
		exit(-1); 
	}
	fin>>rows;
	fin>>cols;
	vector<double> t;
	double ele;
	for(vector<vector<double>>::size_type i=0;i!=rows;i++){
		for(vector<double>::size_type j=0;j!=cols;j++){
			fin>>ele;
			t.push_back(ele);
		}
		data.push_back(t);
		t.clear();
	}
	fin.close();
}

void Matrix::writefile(std::string fname){
	ifstream fexist;
	ofstream fout;
	fexist.open(fname.c_str());
	if(fexist){
		cout<<"Are you sure to OVERWRITE file ["<<fname<<"] of a new matrix?[y/N]"<<endl;
		string exist;
		cin>>exist;
		while(1){
			if(exist=="N"){
				cout<<"Refuse to overwrite"<<endl;
				exit(0);
			}
			else if(exist=="y")
				break;
			else if(exist!="y"&&exist!="N")
				cout<<"WRONG answer! Please enter again"<<endl;
			else;
		}
	}
	fout.open(fname.c_str(),ios::out);
	fout.setf(ios::left);
	fout<<rows<<"    "<<cols<<endl;
	for(vector<vector<double>>::size_type i=0;i!=rows;i++){
		vector<double> rowvec = data[i];
		for(vector<double>::size_type j=0;j!=cols;j++){
			fout.width(14);	
			fout<<rowvec[j];
		}
		fout<<endl;
	}
	fout.close();
}

Matrix Matrix::operator+(Matrix &ad) const{
	Matrix res(rows,cols);
	if(this->rows!=ad.rows||this->cols!=ad.cols){
		cerr<<"Matrix dimensions must agree."<<endl;
		exit(-1);
	}
	else{
		for(int i=0;i<rows;i++)
			for(int j=0;j<cols;j++)
				res.data[i][j] = this->data[i][j]+ad.data[i][j];
	}
	return res;
}

Matrix Matrix::operator-(Matrix &su) const{
	Matrix res(rows,cols);
	if(this->rows!=su.rows||this->cols!=su.cols){
		cerr<<"Matrix dimensions must agree."<<endl;
		exit(-1);
	}
	else{
		for(int i=0;i<rows;i++)
			for(int j=0;j<cols;j++)
				res.data[i][j] = this->data[i][j]-su.data[i][j];
	}
	return res;
}

Matrix Matrix::operator*(Matrix &mul) const{
	if(mul.rows!=cols){
		cerr<<"Matrix dimensions must agree."<<endl;
		exit(-1);
	}
	Matrix res(rows,mul.Getcols());
	for(int i=0;i!=rows;i++){
		for(int j=0;j!=mul.Getcols();j++){
			for(int k=0;k!=cols;k++)
				res.data[i][j] = res.data[i][j]+data[i][k]*mul.data[k][j];
		}
	}
	return res;
}

Matrix Matrix::operator/(const double &div)const{
	Matrix res(rows,cols);
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
			res.data[i][j] = this->data[i][j]/div;
	return res;
}

Matrix Matrix::m_sqrt(){
	Matrix res(rows,cols);
	vector<vector<double>> data_t(rows,vector<double>(cols,0));
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
			data_t[i][j] = sqrt(data[i][j]);
	res.setdata(data_t);
	return res;
}

Matrix Matrix::trans(){
	Matrix res(cols,rows);
	vector<double> t;
	for(int i=0;i!=cols;i++)
		for(int j=0;j!=rows;j++)		
			res.data[i][j] = this->data[j][i];
	return res;
}

vector<double> Matrix::mean(int option){ 
	vector<double> res;
	switch(option){
		case 1:{
			//Mean on rows
			double s=0;
			for(int j=0;j<cols;j++){
				for(int i=0;i<rows;i++)
					s+=this->data[i][j];
				res.push_back(s/rows);
				s=0;
			}
			break;
			   }
		case 2:{
			//Mean on rows
			double s=0;
			for(int j=0;j<rows;j++){
				for(int i=0;i<cols;i++)
					s+=this->data[j][i];
				res.push_back(s/cols);
				s=0;
			}
			break;
			   }
		default:{
			cerr<<"WRONG Input arguments"<<endl;
			exit(-1);
				}	
	}
	return res;
}

Matrix Matrix::sort(vector<int>& id,int option){
	switch(option){
		case 1:{
			//on rows
			Matrix res(rows,cols);
			if(id.size()!=rows){
				cerr<<"Matrix dimensions must agree."<<endl;
				exit(-1);
			}
			for(vector<double>::size_type i=0;i!=id.size();i++)
				res.data[i] = data[id[i]];
			return res;
			   }
		case 2:{
			//on columns
			Matrix res = trans();
			res = res.sort(id,1);
			res = res.trans();
			return res;
			   }
		default:{
			cerr<<"WRONG Input arguments"<<endl;
			exit(-1);
				}
	}

}

Matrix Matrix::trunc(int num,int option){
	switch(option){
		case 1:{
			//on rows
			Matrix res(num,cols);
			vector<vector<double>> d;
			for(int i=0;i<num;i++)
				d.push_back(data[i]);
			res.setdata(d);
			return res; 
			   }
		case 2:{
			//on column
			Matrix res = trans();
			res = res.trunc(num,1);
			res = res.trans();
			return res;
			   }
		default:{
			cerr<<"WRONG Input arguments"<<endl;
			exit(-1);
				}
	}
	
}

Matrix Matrix::invert(){
	if(rows!=cols){
		cerr<<"Matrix dimensions must agree."<<endl;
		exit(-1);
	}
	Matrix res(rows,rows);
	int actualsize = rows;
	int maxsize = rows;
	vector<vector<double>> data_t = data;
	for (int i=1; i<actualsize;i++) 
		data_t[0][i]/= data_t[0][0];
	for (int i=1; i < actualsize;i++){ 
		for (int j=i; j < actualsize;j++){
			double sum = 0.0;
			for (int k = 0;k < i;k++)  
				sum += data_t[j][k] * data_t[k][i];
			data_t[j][i] -= sum;
		}
		if (i == actualsize-1) 
			continue;
		for (int j=i+1;j < actualsize;j++){  
			double sum = 0.0;
			for (int k = 0;k < i;k++)
				sum += data_t[i][k]*data_t[k][j];
			data_t[i][j] = (data_t[i][j]-sum) / data_t[i][i];
		}
	}
	for (int i = 0;i < actualsize;i++)
		for (int j = i;j < actualsize;j++){
			double x = 1.0;
			if ( i!= j ) {
				x = 0.0;
				for (int k = i;k < j;k++) 
					x-= data_t[j][k]*data_t[k][i];
			}
			data_t[j][i] = x/data_t[j][j];
		}
		for (int i = 0;i < actualsize;i++)
			for (int j = i;j < actualsize;j++)  {
				if (i == j) continue;
				double sum = 0.0;
				for (int k = i;k < j; k++)
					sum += data_t[k][j]*( (i==k) ? 1.0 : data_t[i][k] );
				data_t[i][j] = -sum;
			}
		for (int i = 0;i < actualsize;i++) 
			for (int j = 0;j < actualsize;j++)  {
				double sum = 0.0;
				for (int k = ((i>j)?i:j);k < actualsize;k++)  
					sum += ((j==k)?1.0:data_t[j][k])*data_t[k][i];
				data_t[j][i] = sum;
			}
		res.setdata(data_t);
		return res;
}
