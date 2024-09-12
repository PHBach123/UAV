#include <stdio.h>
#include <math.h>

// Table 1, paper 1
double const H = 10;
double const Vmax = 10;            
double const Noise = 0.000000001;
double const S = 50000000;   
double const Ps = 0.04;              
double const Pwpt = 1000;            
double const w0 = 0.001;      
double const alpha = 2.4;         
double const Pb = 0.000001;
double const Pu = 0.005;     
double const DelT = 0.5;
double const muy = 0.84;
double const Nmax = 0.5;
double const E = 0.5772156649;
double const B = 1500000; //mhz
double const I = 0.1;

// Table 1, paper 2
double const W = 100;
double const p = 1.225;
double const R = 0.08;
double const A = 0.79;
double const Ohm = 0.5;
double const s = 0.05;
double const d0 = 0.0151;
double const k = 0.1;
double const del = 0.012;
double const v0 = 7.2;

double const ws[2] = {5,0};   //v? tr� ngu?n
double const wd[2] = {15,0};  //v? tr� d�ch

int const N = 20;  //s? khe th?i gian
double Tn = 0.5;   // 
double x[100],xopt[100];
double y[100],yopt[100];
int jmin[100],jmax[100];

double f22a = 0, f22b = 0, f22d1 = 0, f22d2 = 0,Efly;
double F22a = 0, F22b = 0, F22d1 = 0, F22d2 = 0;
double dsu, ddu,f22amax;

void Calc(double q1, double q2, double d) {
	dsu = sqrt( H*H + (q1 - ws[0]) * (q1 - ws[0]) + (q2 - ws[1]) * (q2 - ws[1]));
	ddu = sqrt( H*H + (q1 - wd[0]) * (q1 - wd[0]) + (q2 - wd[1]) * (q2 - wd[1]));
	f22a = B *DelT*Tn / log(2) * log(1 + (exp(-E) * w0 / Noise) * (Nmax * w0 * Ps + Pu * ceil(sqrt(Noise)) * pow(dsu,alpha)) / pow(dsu,alpha) / pow(ddu,alpha) );
	f22b = B *DelT* Tn / log(2) * log(1 + (exp(-E) * w0 * Ps) / pow(dsu,alpha) / Noise);
	
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
void max()
{ if (F22a>f22amax) {
	int t;
	for(t=0;t<N;t++)
	{xopt[t]=x[t]; yopt[t]=y[t];}
	f22amax=F22a;
}
}
void find(int i,int j){
//double jmin,jmax;
int k;
	i=i+1;x[i]=i; //printf("%d",i);
	if (j-9<0) jmin[i]=0 ; else jmin[i]=j-9;
	if (j+9>20) jmax[i]=20;  else jmax[i]=j+9;
	for(k=jmin[i];k<jmax[i]+1;k++) {
	y[i]=k;
	Calc(x[i], y[i], 10);
	 //printf("\n%lf",f22a);
	//printf("\n%lf %d",sqrt((x[i]-x[N])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N])),i);
	if (i<N) {
	//	printf("\n%lf %d",sqrt((x[i]-x[N])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N])),10*(N-i));
	if	(sqrt((x[i]-x[N-1])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N]))>10*(N-i)) continue;
	if	(sqrt((x[i]-x[N-1])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N]))<10*(N-i)) {
	printf("\n%lf %d",sqrt((x[i]-x[N])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N])),10*(N-i));
	F22a=F22a+f22a; 
	find(x[i],y[i]);
} continue;
	if	(sqrt((x[i]-x[N-1])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N]))>10*(N-i)) continue;
}
	if (i=N) max();
	}
	}

int main () {
	int i,j;
	// init
	// v? tr� d?u
		x[0] = 0;  y[0] = 20;
	// v? tr� cu?i
		x[N] = 20; y[N] = 20;
	find(0,20);
	printf("\n%lf", f22amax);
	for(i=0;i<20;i++) { printf("\n%lf %lf",xopt[i],yopt[i]);
	}
}

	

