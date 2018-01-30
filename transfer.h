#include <cstdlib>

double cmtoev(double);
double amutokg(double);
double cmtoj(double);
double jtocm(double);

inline double cmtoev(double _e0){
	return _e0 / 8065.5409;
}

inline double amutokg(double _m0){
	return 1.660539040E-27 * _m0;
}

inline double cmtoj(double _e0){
	return cmtoev(_e0) * 1.6021766208E-19;
}

inline double jtocm(double _e){
	return _e / 1.6021766208E-19 * 8065.5409;
}
