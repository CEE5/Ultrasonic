#include <iostream>

using namespace std;

#include <math.h>

#include <fstream>

#define inSize 256

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>

#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)



    int x[inSize];
    int y[inSize];

    int g[2*inSize-1];


void four1(double data[], int nn, int isign)
{
    int n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    tempr = data[j];     data[j] = data[i];     data[i] = tempr;
	    tempr = data[j+1]; data[j+1] = data[i+1]; data[i+1] = tempr;
	}
	m = n >> 1;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }
    mmax = 2;
    while (n > mmax) {
	istep = 2*mmax;
	theta = TWOPI/(isign*mmax);
	wtemp = sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j =i + mmax;
		tempr = wr*data[j]   - wi*data[j+1];
		tempi = wr*data[j+1] + wi*data[j];
		data[j]   = data[i]   - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr)*wpr - wi*wpi + wr;
	    wi = wi*wpr + wtemp*wpi + wi;
	}
	mmax = istep;
    }
}


int main()
{

    fstream datei;
    datei.open("output.txt", ios::out);



    for(int i=0;i<inSize;i++){

        x[i]= 0xff * sin(i*PI/40);
        y[i]= 0xff * sin((i-50)*PI/40);

        datei << "x:\t" <<x[i]<<"\t \t y:\t"<<y[i]<<endl;

    }
    datei<<endl<<endl<<"Diskrete Korrelation"<<endl;




    for(int m=0;m<2*inSize-1;m++){
    g[m]=0;
        for(int n=0;n<inSize;n++){

           // g[m]+= x[n]*y[n+(m-(inSize-1))];
            g[m]+= x[n]*y[n+m];
        }
        g[m]=g[m]/inSize;

        datei<<g[m]<<endl;

    }


    datei.close();



    return 0;
}
