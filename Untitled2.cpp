#include <stdio.h>
#include <math.h>

// Table 1, paper 1
double const H = 10;
double const Vmax = 20;            
double const Noise = pow(10,-6);
double const S = 70;   
double const Ps = 0.039810717055;              
double const Pwpt = pow(10,7);            
double const w0 = pow(10,-3);      
double const alpha = 2.3;         
double const Pb = pow(10,-3);
double const Pu = 5;     
double const DelT = 0.5;
double const muy = 0.84;
double const Nmax = 0.5;
double const E = 0.5772156649;
double const B = 20; //mhz
double const I = 0.1;
double const sigma = 0.5;

// Table 1, paper 2
double const W = 0.5;
double const p = 1.225;
double const R = 0.08;
double const A = 0.79;
double const Ohm = 100;
double const s = 0.05;
double const d0 = 0.0151;
double const k = 0.1;
double const del = 0.012;
double v0 = sqrt(W / (2*p*A));

struct Loca {
	double x;
	double y;
};
Loca ws, wd;
int const N = 10;  //s? khe th?i gian
double Tn = 0.5;   // 
Loca q[N];

double f22a = 0, f22b = 0, f22d1 = 0, f22d2 = 0;
double F22a = 0, F22b = 0, F22d1 = 0, F22d2 = 0;
double dsu, ddu;

void Calc(double q1, double q2, double d) {
	dsu = sqrt( H*H + (q1 - ws.x) * (q1 - ws.x) + (q2 - ws.y) * (q2 - ws.y));
	ddu = sqrt( H*H + (q1 - wd.x) * (q1 - wd.x) + (q2 - wd.y) * (q2 - wd.y));

	f22a = B * Tn * DelT / log(2) * log(1 + (exp(-E) * w0 / Noise) * (Nmax * w0 * Ps + Pu * ceil(sigma) * pow(dsu,alpha)) / pow(dsu,alpha) / pow(ddu,alpha) );
	f22b = B * Tn * DelT / log(2) * log(1 + (exp(-E) * w0 * Ps) / pow(dsu,alpha) / Noise);
	
	// Ham nang luong
	double P0 = del / 8 * p * s * A * pow(Ohm,3) * pow(R,3);
	double P1 = (1 + I) * pow(W,3/2) / sqrt(2 * p * A);
	double K1 = 3 / pow(Ohm,2) / pow(R,2);
	double K2 = 1 / 2 / pow(v0,2);
	double K3 = 0.5 * d0 * p * s * A;
	double Efly = P0 * (DelT + K1 * pow(d,2)) + P1 * sqrt(sqrt(pow(DelT,4) + pow(K2,2) * pow(d,4)) - K2 * pow(d,2) ) + K3 * pow(del,3) / pow(DelT,2);
	f22d1 = Efly + Tn * DelT * (Pb + Pu);
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
	q[0].x = 0; q[0].y = 10;
	q[N-1].x = 20; q[N-1].y = 10;
	ws.x = 5; ws.y = 0;
	wd.x = 15; wd.y = 0;
	Calc(q[0].x, q[0].y, 0);
	F22a = f22a;
	F22b = f22b;
	F22d1 = f22d1;
	F22d2 = f22d2;
	
	printf("%lf %lf %lf %lf\n",f22a, f22b, f22d1, f22d2);
	printf("%d %d %d %d",BScheck(), Scheck(), Echeck(), Lcheck(0));
	
	//
	
	
	
	
}
