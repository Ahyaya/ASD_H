#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/ASD_H.h"

int main(void){

//  Test if default option is working
    fprintf(stderr,"testing simple input\n");
    ASD_H_auto("u2u_HANN.test", NULL, "exampleSignals/X_u2n.dat", NULL, NULL);

//  User define options
    fftOption_t opts;
    default_fftOption_init(&opts);

//  Test if smooth filter mode is working
    fprintf(stderr,"testing filter mode 1\n");
    opts.smoothMode=1;
    ASD_H_auto("filter_HANN.test", NULL, "exampleSignals/X_u2n.dat", NULL, &opts);

    opts.winOption.winFuncName=KAISER;
    ASD_H_auto("filter_KAISER.test", NULL, "exampleSignals/X_u2n.dat", NULL, &opts);
    ASD_H_auto("filter_KAISER_n2u.test", NULL, "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

    fprintf(stderr,"testing filter mode 2\n");
    opts.smoothMode=2;
    opts.winOption.winFuncName=HANN;
    ASD_H_auto("n2u_HANN.test", NULL, "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

//
//  start to test N2UviaFile at several frequencies
//
    fprintf(stderr,"============\nStart to test non-uniform to uniform transform.\n============\n");
//  set window function to Kaiser window which is the most widely-used
    opts.smoothMode=0;
    opts.winOption.winFuncName=KAISER;
    opts.winOption.alpha=10;
    opts.tolerance=1e-12;
    ASD_H_auto("n2u_KAISER.test", NULL, "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

//  change window function for further test BLACKMAN
//    opts.winOption.winFuncName=BLACKMAN;
//    ASD_H_auto("n2u_BLACKMAN.test", NULL, "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

//  change window function for further test HANN (default)
    opts.winOption.winFuncName=HANN;
    ASD_H_auto("default.test", NULL, "exampleSignals/X.dat", "exampleSignals/t.dat", NULL);

//  change window function for further test NONE
//    opts.winOption.winFuncName=NONE;
//    opts.debug=1;
//    ASD_H_auto("n2u_NONE.test", NULL, "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);
//
//  start to test U2NviaFile at several frequencies
//
    
    fprintf(stderr,"============\nStart to test uniform to non-uniform transform.\n============\n");
//    opts.winOption.winFuncName=BLACKMAN;
//    opts.dt=0.01;
//    opts.debug=0;
//    ASD_H_auto("u2n_BLACKMAN.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X_u2n.dat", NULL, &opts);

//  default setting has been tested, but mind the sampling time dt=0.001 by default
    opts.winOption.winFuncName=HANN;
    ASD_H_auto("u2n_HANN.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X_u2n.dat", NULL, &opts);

//    opts.winOption.winFuncName=NONE;
//    ASD_H_auto("u2n_NONE.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X_u2n.dat", NULL, &opts);

    opts.winOption.winFuncName=KAISER;
    opts.debug=0;
    ASD_H_auto("u2n_KAISER.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X_u2n.dat", NULL, &opts);

//
//  start to test N2NviaFile at several frequencies
//
    fprintf(stderr,"============\nStart to test non-uniform to non-uniform transform.\n============\n");
    opts.winOption.winFuncName=KAISER;
    opts.winOption.alpha=10;
    opts.tolerance=1e-12;
    opts.debug=0;
    ASD_H_auto("n2n_KAISER.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

//    opts.winOption.winFuncName=BLACKMAN;
//    ASD_H_auto("n2n_BLACKMAN.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

    opts.winOption.winFuncName=HANN;
    ASD_H_auto("n2n_HANN.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

//    opts.winOption.winFuncName=NONE;
//    opts.debug=1;
//    ASD_H_auto("n2n_NONE.test", "exampleSignals/freq_u2n.dat", "exampleSignals/X.dat", "exampleSignals/t.dat", &opts);

    return 0;
}