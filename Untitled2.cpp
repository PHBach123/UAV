#include <stdio.h>
#include <math.h>

double const H = 10;
double const Vmax = 20;            
double const Noise = 0.000000001;
double const S = 50;   
double const Ps = 16;              
double const Pwpt = 1000;            
double const w0 = 0.001;      
double const alpha = 2.3;         
double const Pb = 0.000001;
double const Pu = 8;     
double const DelT = 0.5;
double const muy = 0.84;
double const Nmax = 0.5;
double const E = 0.5772156649;
double const B = 1000000; //mhz

double const ws[2] = {5,0};
double const wd[2] = {15,0};

int const N = 10;

double q[N*N];
double Tn = 0.5;
double f22a = 0, f22b = 0, f22d1 = 0, f22d2 = 0;
double F22a = 0, F22b = 0, F22d1 = 0, F22d2 = 0;
double dsu, ddu;

void Calc(double q1, double q2, double d) {
	dsu = sqrt( H*H + (q1 - ws[0]) * (q1 - ws[0]) + (q2 - ws[1]) * (q2 - ws[1]));
	ddu = sqrt( H*H + (q1 - wd[0]) * (q1 - wd[0]) + (q2 - wd[1]) * (q2 - wd[1]));
	f22a = B * Tn / log(2) * log(1 + (exp(-E) * w0 / Noise) * (Nmax * w0 * Ps + Pu * ceil(sqrt(Noise)) * pow(dsu,alpha)) / pow(dsu,alpha) / pow(ddu,alpha) );
	f22b = B * Tn / log(2) * log(1 + (exp(-E) * w0 * Ps) / pow(dsu,alpha) / Noise);
	double Efly = 0;
	f22d1 = 0;
	f22d2 = muy * (1 - Tn) * DelT * w0 * Pwpt / pow(dsu, alpha);
}

int BScheck() {
	if (F22b + f22b + sqrt(Noise) * S < F22a + f22a)
	   return 0;
	   else return 1;
}

int Scheck() {
	if (F22a + f22a < S)
	   return 0;
	   else return 1;
}

int Echeck() {
	if (F22d1 + f22d1 > F22d2 + f22d2)
	   return 0;
	   else return 1;
}

int Lcheck(int d) {
	if (d > (Vmax * DelT))
	   return 0;
 	   else return 1;
}

int main () {
	
	// init
	
	q[0] = 0; 
	q[N] = 20;
	q[N-1] = 10; 
	q[(N-1)*(N-1)] = 10;
	Calc(q[0], q[N], 0);
	F22a = f22a;
	F22b = f22b;
	F22d1 = f22d1;
	F22d2 = f22d2;
	
	printf("%lf %lf %lf %lf",f22a, f22b, dsu, ddu);
	
	//
	
}
