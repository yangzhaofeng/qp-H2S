#include <iostream>
#include <cstdlib>
#include <armadillo>
#include "qpcall.h"
//#include "transfer.h"

using namespace std;
using namespace arma;

/*
const double m_H = 1.0078250319; //amu, Pure and Applied Chemistry. 75 (6): 683–800
const double m_S = 31.97207073; //amu, Pure and Applied Chemistry. 75 (6): 683–800
const double De = 38667.2857; //cm^-1
const double alpha = 1.6627E10; //m^-1
const double frr = -9.8705E22; //cm^-1 m^-2
const double h = 6.626069934E-34; //J s
const double h_bar = 1.054571800E-34; //J s
*/


int main(){
	cout <<omega<<endl<<omegax<<endl<<grr<<endl<<frr<<endl;
	mat E0(8,8,fill::zeros);
	for(int _m=0;_m<=7;_m++){
		E0(_m,_m)=E(_m);
	}
	//cout<<E0<<endl;

	mat p1(8,8,fill::zeros);
	mat p2(8,8,fill::zeros);
	mat x1(8,8,fill::zeros);
	mat x2(8,8,fill::zeros);
	for (int _m=0;_m<=7;_m++){
		for (int _k=0;_k<=7;_k++){
			p1(_m,_k) = sqrt(_k+1) * delta(_m,_k+1) - sqrt(_k) * delta(_m,_k-1);
			p2(_m,_k) = sqrt(_k+1) * delta(_m,_k+1) - sqrt(_k) * delta(_m,_k-1);
			x1(_m,_k) = sqrt(_k+1) * delta(_m,_k+1) + sqrt(_k) * delta(_m,_k-1);
			x2(_m,_k) = sqrt(_k+1) * delta(_m,_k+1) + sqrt(_k) * delta(_m,_k-1);
		}
	}
	//cout <<p1<<"\n"<<p2<<"\n"<<x1<<"\n"<<x2<<endl;
	mat H1(8,8,fill::zeros);
	H1 = grr * -0.5 * amutokg(m) * omega * h_bar * p1 * p2 + cmtoj(frr) * 0.5 * h_bar / (amutokg(m) * omega) * x1 * x2;
	//cout<<H1<<endl;

	//mat H = E0 + H1;
	//cout<<diagmat(eig_sym(H))<<endl;
	mat etemp = diagmat(eig_sym(H1));
	mat E1(8,8,fill::zeros);
	
	for(int _k=0;_k<=7;_k++){ //reverse the diag
		E1(_k,_k)=etemp(7-_k,7-_k);
	}
	//cout<<E1<<endl;

	mat Ex(8,8);
	Ex = E0 + E1;
	cout<<Ex / 1.6021766208E-19 * 8065.5409<<endl;
	mat Ey(8,8,fill::zeros);
	for(int _m=0;_m<=7;_m++){
		Ey(_m,_m)=Ex(0,0);
	}
	cout<<(Ex - Ey) / 1.6021766208E-19 * 8065.5409<<endl;
	return 0;
}
