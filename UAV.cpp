#include <stdio.h>
#include <math.h>

// Table 1, paper 1
double const H = 10;
double const Vmax = 6.4;            
double const Noise = 0.000000001;
double const S = 35;   
double const Ps = 0.04;              
double const Pwpt = pow(10,6.3);            
double const w0 = 0.001;      
double const alpha = 2.3;         
double const Pb = 0.000001;
double const Pu = 0.005;     
double const DelT = 0.5;
double const muy = 0.84;
double const Nmax = 0.5;
double const E = 0.5772156649;
double const B = 4; //mhz
double const I = 0.1;

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
double const v0 = 7.2;

double const ws[2] = {2,0};   //vi tri nguon
double const wd[2] = {7,0};  //vi tri dich

double Tn = 0.5;   
int const N = 10; //So khe thoi gian
double const iV=floor(sqrt(Tn*Vmax*Tn*Vmax-1));  
double x[100],xopt[100],sum[100];
double y[100],yopt[100],F22_a[100],F22_b[100],F22_d1[100],F22_d2[100];
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
	f22d1 = Efly+   Tn * DelT * (Pb + Pu);
	f22d2 = muy * (1 - Tn) * DelT * w0 * Pwpt / pow(dsu, alpha);
}

int BScheck(int i) {     //Data tren UAV lon hon Data tai dich
	if (F22_b[i-1] + f22b + sqrt(Noise) * S < F22_a[i-1] + f22a)
	   return 0;
	   else return 1;
}

int Scheck() {          //Data toi thieu
	if (F22_a[N] < S)
	   return 0;
	   else return 1;
}

int Echeck(int i) {     // Rang buoc nang luong
	if (F22_d1[i-1] + f22d1 > F22_d2[i-1] + f22d2)
	   return 0;
	   else return 1;
}

void max()
{
Calc(10, 10, sqrt((10-x[9])*(10-x[9])+(10-y[9])*(10-y[9])));
F22_a[N]=F22_a[N-1]+f22a;
if (F22_a[N]>f22amax) {    // Cap nhat gia tri cuc dai
	if (Scheck()){
	for(int t=0;t<N;t++)
	{xopt[t]=x[t]; yopt[t]=y[t];}
	f22amax=F22_a[N];    
}
}
}
void find(int i,int j){
	i=i+1; 
	if (j-iV<0) jmin[i]=0 ; else jmin[i]=j-iV;
	if (j+iV>10) jmax[i]=10;  else jmax[i]=j+iV;
	if (i==N) max();
 	if (i<N) {for(int k=jmin[i];k<=jmax[i];k++) 
	{x[i]=i;
	y[i]=k;
	Calc(x[i], y[i], sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1])));   //  Tinh toan thong so
	if	((sqrt((x[i]-x[N])*(x[i]-x[N]) + (y[i]-y[N])*(y[i]-y[N]))<Tn*Vmax*(N-i))&&BScheck(i)&&Echeck(i)){  // Kiem tra cac rang buoc
	F22_a[i]=F22_a[i-1]+f22a; 
	F22_b[i]=F22_b[i-1]+f22b;
	F22_d1[i]=F22_d1[i-1]+f22d1;
	F22_d2[i]=F22_d2[i-1]+f22d2;
	find(x[i],y[i]);}
	};
}}

int main () {
	// init
	// vi tri dau
		x[0] = 0;  y[0] = 10;
	// vi tri cuoi
		xopt[N]=x[N] = 10; yopt[N]=y[N] = 10;
	find(0,10);
	printf("\n%lf", f22amax);
	for(int i=0;i<=N;i++) { printf("\n%lf %lf",xopt[i],yopt[i]);
	};
}
