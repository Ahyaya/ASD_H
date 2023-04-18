#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static const double pi = 3.141592653589793238463;

int main(void){
    int64_t pf, N = 1e4;            // number of nonuniform points
    FILE *fptr_t=fopen("exampleSignals/t.dat","w");
    FILE *fptr_x=fopen("exampleSignals/X.dat","w");
    FILE *fptr_x_u2n=fopen("exampleSignals/X_u2n.dat","w");
    FILE *fptr_freq_u2n=fopen("exampleSignals/freq_u2n.dat","w");

    double t,X;

    for (pf=0; pf<N; ++pf) {
        t= 12.345678+(pf+0.1*rand()/RAND_MAX)/100.0;   //  introduce sampling noise in timestamp
        X = 1024*cos(2*pi*t*1.0+1.234) + 256*cos(2*pi*t*2.0+0.123) + 64*cos(2*pi*t*4.0) + 16*cos(2*pi*t*8.0+5.67) + 4*cos(2*pi*t*16.0-5.67);
        fprintf(fptr_t,"%.6e\n",t);
        fprintf(fptr_x,"%.6e\n",X);
    }
    for (pf=0; pf<N; ++pf) {
        t=pf/100.0;
        X = 1024*cos(2*pi*t*1.0+1.234) + 256*cos(2*pi*t*2.0+0.123) + 64*cos(2*pi*t*4.0) + 16*cos(2*pi*t*8.0+5.67) + 4*cos(2*pi*t*16.0-5.67);
        fprintf(fptr_x_u2n,"%.6e\n",X);
    }
    int64_t Nf=2e3;
    for(pf=0;pf<Nf;++pf){
        fprintf(fptr_freq_u2n,"%.6e\n",(pf+1)/100.0);
    }
    fclose(fptr_t); fclose(fptr_x); fclose(fptr_x_u2n); fclose(fptr_freq_u2n);
    fprintf(stderr,"X.dat t.dat X_u2n.dat are updated in ./exampleSignals\ntotal lines: %ld\n", N);
    return 0;
}