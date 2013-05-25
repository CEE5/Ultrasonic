#include <iostream>

using namespace std;

#include <math.h>
#define PI 3.14159265
#include <fstream>


#include <complex.h>


#define inSize 1024



void IFFT(double _Complex *data, unsigned N)
{
    unsigned butterflySize;
// size for actual butterfly calculation
    int i, j, k;
// local index variables
    double _Complex wActual; // actual rotation factor
    double _Complex wStep;
// step rotation factors
    double _Complex tmp;
// temp. value for butterfly calculation
// loop over all level of FFT
    for(butterflySize=N/2; butterflySize>0; butterflySize/=2)
    {
// evaluate angle step and set first angle
        wStep = cos(+PI/butterflySize)+ sin(+PI/butterflySize)*_Complex_I;
        wActual = 1 +0 * _Complex_I;
// loop over number of butterflys
        for(j=0; j<butterflySize; j++)
        {
// loop over number of FFTs
            for(i=j; i<N; i+=2*butterflySize)
            {
// get index of second element
                k = i+butterflySize;
// perform butterfly calculation
                tmp = data[i];
// store one element
                data[i] += data[k];
// take sum
                data[k] = tmp-data[k]; // take difference
                data[k] *= wActual;
// multiply with rotation factor
            }
// evaluate next rotation factor
            wActual *= wStep;
        }
    }
// perform bit reversal
    j = 0;
    for(i=0; i<N; i++)
    {
        if(j>i)
        {
// swap numbers
            tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
        }
        k = N/2;
        while(k>=2 && j>=k)
        {
            j -= k;
            k /= 2;
        }
        j += k;
        data[i]=data[i]/N;
    }
}
void FFT(double _Complex *data, unsigned N)
{
    unsigned butterflySize;
// size for actual butterfly calculation
    int i, j, k;
// local index variables
    double _Complex wActual; // actual rotation factor
    double _Complex wStep;
// step rotation factors
    double _Complex tmp;
// temp. value for butterfly calculation
// loop over all level of FFT
    for(butterflySize=N/2; butterflySize>0; butterflySize/=2)
    {
// evaluate angle step and set first angle
        wStep = cos(-PI/butterflySize)+ sin(-PI/butterflySize)*_Complex_I;
        wActual = 1 +0 * _Complex_I;
// loop over number of butterflys
        for(j=0; j<butterflySize; j++)
        {
// loop over number of FFTs
            for(i=j; i<N; i+=2*butterflySize)
            {
// get index of second element
                k = i+butterflySize;
// perform butterfly calculation
                tmp = data[i];
// store one element
                data[i] += data[k];
// take sum
                data[k] = tmp-data[k]; // take difference
                data[k] *= wActual;
// multiply with rotation factor
            }
// evaluate next rotation factor
            wActual *= wStep;
        }
    }
// perform bit reversal
    j = 0;
    for(i=0; i<N; i++)
    {
        if(j>i)
        {
// swap numbers
            tmp = data[i];
            data[i] = data[j];
            data[j] = tmp;
        }
        k = N/2;
        while(k>=2 && j>=k)
        {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}



int main()
{
    double _Complex x[inSize];
    double _Complex y[inSize];

    int g[2*inSize-1];

    fstream datei;
    datei.open("output.txt", ios::out);



    for(int i = 1; i<inSize; i++)
    {

        x[i]= 0xff * sin(double(i)/80*PI*2);
        y[i]= 0xff * sin(double(i)/80*PI*2);

       // datei << "x:\t" <<x[i]<<"\t \t y:\t"<<y[i]<<endl;

    }

    /*
    datei<<endl<<endl<<"Diskrete Korrelation"<<endl;




    for(int m=0; m<2*inSize-1; m++)
    {
        g[m]=0;
        for(int n=0; n<inSize; n++)
        {

            // g[m]+= x[n]*y[n+(m-(inSize-1))];
            g[m]+= x[n]*y[n+m];
        }
        g[m]=g[m]/inSize;

        datei<<g[m]<<endl;

    }*/

    FFT(x,1024);
    IFFT(x,1024);

    for(int i=0; i<inSize; i++)
    {

        datei << creal(x[i])<<endl;
    }






    datei.close();



    return 0;
}
