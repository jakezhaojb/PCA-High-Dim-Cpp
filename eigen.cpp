#include "eigen.h"

inline void sortArray(valarray<double>& x, vector<int> &id)
{   
	int i,j,k,t1;
	int n=x.size();
	double t2;
	for(i=0;i<n;i++)
		id.push_back(i);
	for(i=0;i<n-1;i++)
	{
		k=i;
		for(j=i+1;j<n;j++) 
			if(x[j]>x[k]) 
				k=j;
		if(k!=i)
		{
			t2=x[i]; 
			x[i]=x[k]; 
			x[k]=t2;
			t1=id[i];
			id[i]=id[k];
			id[k]=t1;
		}
	}
}

inline bool Symmetry(const Matrix &m)
{
	bool judge = true;
	int row = m.Getrows();
	if(m.Getcols() == row)
	{
		for(int i = 1; i < row; i++)			
			for(int j = 0; j < i; j++)
				if(abs(m(i,j) - m(j,i)) >= 1.0e-30)
				{
					judge = false;
					return judge;
				}
	}
	else
		judge = false;
	return judge;	
}

int Householder(Matrix m, Matrix &q, valarray<double>& main_diag, valarray<double>& minor_diag)
{
	int i,j, k, Rank;
	double h, f, g, h2;
	if(Symmetry(m)!=true)
		return 0;
	Rank = m.Getcols();
	for(i=0; i<Rank; i++)
		for(j=0; j<Rank; j++)
			q(i,j)=m(i,j);
	for(i=Rank-1; i>=1; i--)
	{ 
		h=0.0;
		if(i>1)
			for(k=0; k<i; k++)
				h=h+q(i,k)*q(i,k);
		if((abs(h - 0) < 1.0e-15))
		{
			minor_diag[i]=0.0;
			if(i==1) minor_diag[i]=q(i,i-1);
			main_diag[i]=0.0;
		}
		else
		{ 
			minor_diag[i]=sqrt(h);
			if(q(i,i-1)>0.0) minor_diag[i]=-minor_diag[i];
			h=h-q(i,i-1)*minor_diag[i];
			q(i,i-1)=q(i,i-1)-minor_diag[i];
			f=0.0;
			for(j=0; j<i; j++)
			{ 
				q(j,i) = q(i,j) / h;
				g = 0.0;
				for(k=0; k<=j; k++)
					g = g + q(j,k) * q(i,k);
				if(j+1<i)
					for(k=j+1; k<i; k++)
						g = g + q(k,j) * q(i,k);
				minor_diag[j] = g / h;
				f = f + g * q(j,i);
			}
			h2 = f / (h+h);
			for(j=0; j<i; j++)
			{ 
				f = q(i,j);
				g = minor_diag[j] -h2 * f;
				minor_diag[j] = g;
				for(k=0; k<=j; k++)
					q(j,k)=q(j,k)-f*minor_diag[k]-g*q(i,k);
			}
			main_diag[i]=h;
		}
	}
	for(i=0; i<Rank-1; i++) minor_diag[i]=minor_diag[i+1];
	minor_diag[Rank-1]=0.0;
	main_diag[0]=0.0;
	for(i=0; i<Rank; i++)
	{ if((main_diag[i]!=0.0)&&(i-1>=0))
	for(j=0; j<i; j++)
	{ 
		g=0.0;
		for(k=0; k<i; k++)
			g=g+q(i,k)*q(k,j);
		for(k=0; k<i; k++)
			q(k,j)=q(k,j)-g*q(k,i);
	}
	main_diag[i]=q(i,i);
	q(i,i)=1.0;
	if(i-1>=0)
		for(j=0; j<i; j++)
		{ 
			q(i,j)=0.0;
			q(j,i)=0.0;
		}
	}
	return 1;
}


int QR(valarray<double>& main_diag, valarray<double>& minor_diag, Matrix &q, double accuracy, int l)
{
	int i, j, k, m, it, Rank;
	double h, g, p, r, e, s, d(0), f(0);
	Rank = q.Getcols();
	minor_diag[Rank-1]=0.0;
	for(int j=0; j<Rank; j++)
	{ 
		it=0;
		h=accuracy*(abs(main_diag[j])+abs(minor_diag[j]));
		if(h>d) d=h;
		m=j;
		while((m<Rank)&&(abs(minor_diag[m])>d)) m++;
		if(m!=j)
		{
			do
			{
				if(it==l)
				{ 
					return 0;
				}
				it++;
				g=main_diag[j];
				p=(main_diag[j+1]-g)/(2.0*minor_diag[j]);
				r=sqrt(p*p+1.0);
				if(p>0.0 || (abs(p - 0) < 1.0e-15)) main_diag[j]=minor_diag[j]/(p+r);
				else main_diag[j]=minor_diag[j]/(p-r);
				h=g-main_diag[j];
				for(i=j+1; i<Rank; i++)	main_diag[i]=main_diag[i]-h;
				f=f+h; 
				p=main_diag[m]; 
				e=1.0; 
				s=0.0;
				for(i=m-1; i>=j; i--)
				{
					g=e*minor_diag[i];
					h=e*p;
					if(abs(p)>=abs(minor_diag[i]))
					{
						e=minor_diag[i]/p;
						r=sqrt(e*e+1.0);
						minor_diag[i+1]=s*p*r;
						s=e/r; 
						e=1.0/r;
					}
					else
					{
						e=p/minor_diag[i]; 
						r=sqrt(e*e+1.0);
						minor_diag[i+1]=s*minor_diag[i]*r;
						s=1.0/r; 
						e=e/r;
					}
					p=e*main_diag[i]-s*g;
					main_diag[i+1]=h+s*(e*g+s*main_diag[i]);
					for(k=0; k<Rank; k++)
					{ 
						h=q(k,i+1);
						q(k,i+1)=s*q(k,i)+e*h;
						q(k,i)=e*q(k,i)-s*h;
					}
				}
				minor_diag[j]=s*p;
				main_diag[j]=e*p;
			}while(abs(minor_diag[j])>d);
		}
		main_diag[j]=main_diag[j]+f;
	}
	for(i=0; i<Rank; i++)
	{ 
		k=i;
		p=main_diag[i];
		if(i+1<Rank)
		{
			j=i+1;
			while((j<Rank)&&(main_diag[j]<=p))
			{
				k=j;
				p=main_diag[j];
				j=j+1;
			}
		}
		if(k!=i)
		{ 
			main_diag[k]=main_diag[i]; 
			main_diag[i]=p;
			for(j=0; j<Rank; j++)
			{
				p=q(j,i);
				q(j,i)=q(j,k);
				q(j,k)=p;
			}
		}
	}
	return 1;
}

bool Matrix::eig(vector<double> &eigval,Matrix *eigvec){
	bool flag = false;
	if(this->rows!=this->cols){
		cout<<"Inappropriate Input Matrix! Only SQUARE Matrix Permitted In EIGEN Calculating"<<endl;
		return false;
	}
	int size = rows;
	valarray<double> main_diag(size), minor_diag(size);
	Matrix q(size,size);
	Matrix digenvalue(size,size);
	digenvalue.setdata(data);
	double accuracy = 1.0e-8;	
	int k=Householder(digenvalue,q,main_diag,minor_diag);
	if(k==0){
		cout<<"The Matrix is not SYMMETRICAL"<<endl; 
		return false;
	}
	k = QR(main_diag,minor_diag,q,accuracy,160);
	if(k<0){
		cout<<"Calculating failed in Required Iteration Times"<<endl; 
		return false;
	}
	vector<int> id;
	sortArray(main_diag,id);			
	for(int j=0;j<size;j++)
		eigval.push_back(main_diag[j]);
	q = q.sort(id,2);
	*eigvec = q;
	flag = true;
	return flag;
}
