#include <iostream>

using namespace std;

#include <math.h>
#define PI 3.14159265
#include <fstream>

#define inSize 256

    int x[inSize];
    int y[inSize];

    int g[2*inSize-1];


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
