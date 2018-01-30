#include<cmath>
#include"transfer.h"

const double _pi = 3.1415926535897932384626;

const double m_H = 1.0078250319; //amu, Pure and Applied Chemistry. 75 (6): 683–800
const double m_S = 31.97207073; //amu, Pure and Applied Chemistry. 75 (6): 683–800
const double De = 38667.2857; //cm^-1
const double alpha = 1.6627E10; //m^-1
const double frr = -9.8705E22; //cm^-1 m^-2
const double h = 6.626069934E-34; //J s
const double h_bar = 1.054571800E-34; //J s
const double beta = 92.11; //deg
const double m = m_S * m_H / (m_S + m_H); //amu
const double omega = alpha * sqrt(2.0 * cmtoj(De) / amutokg(m));
const double omegax = alpha * alpha * h_bar / (2 * amutokg(m));
const double k = omega / omegax;
const double grr = cos(beta * _pi /180) / amutokg(m_S);

/*class _laguerre
{
private:
	int __n; double __alpha, __x, _result;
public:
	void init(){
		__n=0; __alpha=0; __x=0; _result= -1;
	}
	bool same(int _n,double _alpha,double _x){
		return (_n==__n && _alpha==__alpha && _x==__x && result != -1);
	}
	void flush(){
		result = -1;
	}
	void result_write(int _n; double _alpha,double _x,double result){
		__n=_n; __alpha=_alpha; __x=_x; _result=result;
	}
	double result_get(){
		return _result;
	}
}lagu[50];*/

//double omega(double,double,double);
//double omega();
//double omegax(double,double,double);
//double omegax();
//double k();
//double grr();
//double grr(double,double);
double Nn();
double Nn(double,double,int);
int factorial(int);
double laguerre(int,double,double);
void laguerre_init();
int delta(int,int);
double E(int);

/*
#if (defined alpha) && (defined de) && (defined m)
inline double omega(){
	return omega(alpha,De,m);
}
#endif

inline double omega(double _alpha, double _De, double _m){
	return (_alpha * sqrt(2.0 * _De / _m));
}

#if (defined alpha) && (defined hbar) && (defined m)
inline double omegax(){
	return omegax(alpha, hbar, _m);
}
#endif

inline double omegax(double _alpha, double _hbar, double _m){
	return (_alpha * _alpha * _hbar / ( 2 * _m ));
}

#if (defined beta) && (defined m_S)
inline double grr(){
	return grr(beta,m_S);
}
#endif


inline double grr(double _beta,double _m_S){
	return cos(_beta * _pi /180) / _m_S;
}
*/

inline double Nn(double _alpha, double _k, int _n){
	return sqrt(_alpha * factorial(_n) * (_k - 2 * _n - 1) / tgamma(_k - _n));
}

inline int factorial(int _n){
	return tgamma(_n + 1);
}

inline double E(int n){
	return h_bar * omega * (n + 0.5) - h_bar * omegax * (n + 0.5) * (n + 0.5);
}

/*double laguerre(int _n,double _alpha,double _x){
	if(_n==0){
		lagu[_n].result_write(_n,_alpha,_x,1);
	}
	if(_n==1){
		lagu[_n].result_write(_n,_alpha,_x,1 + _alpha - _x);
	}
	if(lagu[_n].same(_n,_alpha,_x)){
		return lagu[_n].result_get();
	}else{
		return ((2*(_n - 1) + 1 + _alpha - _x) * laguerre(_n - 1,_alpha,_x) - (_n - 1 + _alpha) * lagurre(_n - 2,_alpha,_x)) / _n;
	}
}*/

/*void laguerre_init(){
	for(int __i=0;__i<=49;__i++){
		lagu[i].init();
	}
}*/

inline int delta(int _n,int _k){
	return (_n==_k);
}
