#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <sys/stat.h>
#include "../include/finufft.h"
#include "../include/ASD_H.h"

static const double pi=3.141592653589793238463;

/* PURPOSE: The modified Bessel function of the first kind I_0(x). */
double bessel0(double x){
    double ax, ans;
    double y;
    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75;
        y = y * y;
        ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }else{
        y = 3.75 / ax;
        ans = (exp(ax) / sqrt(ax))*(0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))));
    }
    return ans;
}

/* PURPOSE: The modified Bessel function of the first kind I_1(x). */
double bessel1(double x){
    double ax, ans;
    double y;
    if ((ax = fabs(x)) < 3.75) {
        y = x / 3.75, y = y * y;
        ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    }else{
        y = 3.75 / ax;
        ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

double quickmedian(double array[], int64_t length){
    double gmin=array[0], gmax=array[0], segmentLen, segmentSum[4096]={0};
    int64_t ptr, segmentCount[4096]={0}, Count;
    int sPtr;
    for(ptr=1;ptr<length;++ptr){
        gmax=array[ptr]>gmax?array[ptr]:gmax;
        gmin=array[ptr]<gmin?array[ptr]:gmin;
    }
    gmax+=0.0001*(gmax-gmin);
    gmin-=0.0001*(gmax-gmin);
    segmentLen=(gmax-gmin)/4096;
    for(ptr=0;ptr<length;++ptr){
        sPtr=(int) ((array[ptr]-gmin)/segmentLen);
        segmentSum[sPtr]+=array[ptr];
        segmentCount[sPtr]+=1;
    }
    for(sPtr=0, Count=0; sPtr<4096; ++sPtr){
        Count+=segmentCount[sPtr];
        if(Count > (length>>1)) break;
    }
    return segmentSum[sPtr]/segmentCount[sPtr];
}

int default_fftOption_init(fftOption_t *opts){
    opts->winOption.winFuncName=HANN;
    opts->winOption.alpha=10;
    opts->winLengthRatio=1;
    opts->overlapRatio=0;
    opts->smoothMode=0;
    opts->smoothNodes=256;
    opts->debug=0;
    opts->tolerance=1e-9;
    opts->dt=0.01;
    return 0;
}

int nufft_winFunc_create(double *winFunc, int64_t N, double *t, fftWinFuncOpt_t *opts){
    fftWinFuncOpt_t optParse;
    if(N<1){
        fprintf(stderr,"window function cannot be created with N<1.\n");
        return -1;
    }
    if(opts==NULL){
        optParse.winFuncName=0;
    }else{
        optParse.winFuncName=opts->winFuncName;
        optParse.alpha=opts->alpha;
    }
    int64_t pf;
    double *buff = (double *)malloc(sizeof(double)*(2*N+1));
    double *diff_t = buff+N, *win_t = buff+N;
    double sumWT=0;
    if(t==NULL){
        t=buff;
        for(pf=0;pf<N;++pf){
            t[pf]=(double)pf;
            win_t[pf]=1;
        }
    }else{
        for(pf=1;pf<N;++pf){
            diff_t[pf] = t[pf] - t[pf-1];
        }
        diff_t[0] = diff_t[1];
        diff_t[N] = diff_t[N-1];
        for(pf=0;pf<N;++pf){
            win_t[pf] = 0.5*(diff_t[pf]+diff_t[pf+1]);
            sumWT += win_t[pf];
        }
        for(pf=0;pf<N;++pf){
            win_t[pf] *= N/sumWT;
        }
    }

    double t_duration=t[N-1]-t[0], sumCG=0;
    double sigma, tmp_0, tmp_1;
    switch (optParse.winFuncName)
    {
    case -1://NONE
        for(pf=0;pf<N;++pf){
            winFunc[pf]=1;
        }
        break;
    case 0://HANN
        for(pf=0;pf<N;++pf){
            winFunc[pf]=0.5-0.5*cos(2*pi*(t[pf]-t[0])/t_duration);
        }
        break;
    case 1://HAMM
        for(pf=0;pf<N;++pf){
            winFunc[pf]=0.54-0.46*cos(2*pi*(t[pf]-t[0])/t_duration);
        }
        break;
    case 2://FLAT
        for(pf=0;pf<N;++pf){
            winFunc[pf]=0.21557895-0.41663158*cos(2*pi*(t[pf]-t[0])/t_duration)+0.277263158*cos(4*pi*(t[pf]-t[0])/t_duration)-0.083578947*cos(6*pi*(t[pf]-t[0])/t_duration)+0.006947368*cos(8*pi*(t[pf]-t[0])/t_duration);
        }
        break;
    case 3://BLACKMAN
        for(pf=0;pf<N;++pf){
            winFunc[pf]=7938.0/18608-8240.0/18608*cos(2*pi*(t[pf]-t[0])/t_duration)+1430.0/18608*cos(4*pi*(t[pf]-t[0])/t_duration);
        }
        break;
    case 4://POISSON
        for(pf=0;pf<N;++pf){
            winFunc[pf]=exp(-fabs((t[pf]-t[0])-0.5*t_duration)/t_duration);
        }
        break;
    case 5://KAISER
        for(pf=0;pf<N;++pf){
            winFunc[pf]=bessel0(optParse.alpha*pi*sqrt(1-(2*(t[pf]-t[0])/t_duration-1)*(2*(t[pf]-t[0])/t_duration-1)))/bessel0(optParse.alpha*pi);
        }
        break;
    case 6://NORM,GAUSS
        for(pf=0;pf<N;++pf){
            winFunc[pf]=exp(-0.5*((2*(t[pf]-t[0])/t_duration-1)*optParse.alpha/3)*((2*(t[pf]-t[0])/t_duration-1)*optParse.alpha/3));
        }
        break;
    case 7://WELCH
        for(pf=0;pf<N;++pf){
            winFunc[pf]=1-(2*(t[pf]-t[0])-t_duration)*(2*(t[pf]-t[0])-t_duration)/(t_duration*t_duration);
        }
        break;
    case 8://NUTTALL
        for(pf=0;pf<N;++pf){
            winFunc[pf]=0.355768-0.487396*cos(2*pi*(t[pf]-t[0])/t_duration)+0.144232*cos(4*pi*(t[pf]-t[0])/t_duration)-0.012604*cos(6*pi*(t[pf]-t[0])/t_duration);
        }
        break;
    case 9://PCOS
        for(pf=0;pf<N;++pf){
            winFunc[pf]=pow(cos(pi*(t[pf]-t[0])/t_duration-pi/2), opts->alpha);
        }
        break;
    case 10://PLANCK
        sigma=1/opts->alpha;
        sigma=sigma<0?0.3:sigma;
        sigma=sigma<1?sigma:0.3;
        sigma=sigma<0.5?sigma:1-sigma;
        tmp_0=sigma*t_duration;
        tmp_1=(1-sigma)*t_duration+1;
        winFunc[0]=0;
        winFunc[N-1]=0;
        for(pf=1;t[pf]<tmp_0;++pf){
            winFunc[pf]=1/(exp(tmp_0*(1/t[pf]+1/(t[pf]-tmp_0)))+1);
        }
        for(;t[pf]<tmp_1;++pf){
            winFunc[pf]=1;
        }
        for(;pf<N;++pf){
            winFunc[pf]=1/(exp(tmp_0*(1/(t_duration-t[pf])+1/(tmp_1-t[pf])))+1);
        }
        break;
    default:
        fprintf(stderr,"Unknown window function type, using Hanning window by default.\n");
        for(pf=0;pf<N;++pf){
            winFunc[pf]=0.5-0.5*cos(2*pi*(t[pf]-t[0])/t_duration);
        }
        break;
    }
    /*Compute coherent gain's level*/
    for(pf=0;pf<N;++pf){
        sumCG += winFunc[pf];
    }
    for(pf=0;pf<N;++pf){
        winFunc[pf] *= win_t[pf]*N/sumCG;
    }
    free(buff);
    return 0;
}

int ASD_filter(fftOutput_t *Result, int N_stack){
    N_stack=N_stack>256?N_stack:256;
    int64_t pf, mpf, Nf=Result->N, pstack, qstack, nstack, hstack;
    double *freq=Result->freq, *ASD=Result->ASD, bndfreq[N_stack];
    while(freq[0]<1e-8){
        ++freq;
        ++ASD;
        --Nf;
    }
    double expo_0=log(freq[0]), dexpo=(log(freq[Nf-1])-log(freq[0]))/N_stack, sumstackASD, sumstackFreq, meanstackASD, meanstackFreq, peakASD, peakfreq, medianASD;
    for(pf=0;pf<N_stack;++pf){
        bndfreq[pf]=exp(expo_0+pf*dexpo);
    }
    for(pf=0, mpf=0, pstack=0, nstack=0; pf<Nf; ++nstack, pf+=pstack){
        for(pstack=0,sumstackASD=0,sumstackFreq=0,peakASD=0;pf+pstack<Nf && freq[pf+pstack]<bndfreq[nstack+1];++pstack){
            sumstackASD+=ASD[pf+pstack];
            sumstackFreq+=freq[pf+pstack];
            if(ASD[pf+pstack]>peakASD){
                peakASD=ASD[pf+pstack];
                peakfreq=freq[pf+pstack];
            }
        }
        if(pstack<8){
            for(qstack=0;qstack<pstack;++qstack){
                Result->fltFreq[mpf]=freq[pf+qstack];
                Result->fltASD[mpf]=ASD[pf+qstack];
                ++mpf;
            }
            continue;
        }
        medianASD=quickmedian(&ASD[pf],pstack);
        meanstackASD=sumstackASD/pstack;
        meanstackFreq=sumstackFreq/pstack;
        if(peakASD>8*medianASD){
            hstack=pstack>>1;
            if(peakfreq>freq[pf+hstack]){
                Result->fltFreq[mpf]=freq[pf+hstack];Result->fltASD[mpf]=ASD[pf+hstack];++mpf;
                Result->fltFreq[mpf]=peakfreq;Result->fltASD[mpf]=peakASD;++mpf;
            }else{
                Result->fltFreq[mpf]=peakfreq;Result->fltASD[mpf]=peakASD;++mpf;
                Result->fltFreq[mpf]=freq[pf+hstack];Result->fltASD[mpf]=ASD[pf+hstack];++mpf;
            }
        }else{
            Result->fltFreq[mpf]=meanstackFreq;Result->fltASD[mpf]=medianASD;++mpf;
        }
        
    }
    Result->fltN=mpf;
    return 0;
}

int ASD_H_N2U(fftOutput_t *Result, int64_t N, double *X, double *t, fftOption_t *opts){
    fftOption_t optParse;
    default_fftOption_init(&optParse);
    if(opts==NULL){
        //Don't know what to do :D
    }else{
        optParse.winOption.winFuncName=opts->winOption.winFuncName;
        optParse.winOption.alpha=opts->winOption.alpha;
        optParse.winLengthRatio=opts->winLengthRatio;
        optParse.overlapRatio=opts->overlapRatio;
        optParse.smoothMode=opts->smoothMode;
        optParse.smoothNodes=opts->smoothNodes;
        optParse.debug=opts->debug;
        optParse.tolerance=opts->tolerance;
        optParse.dt=opts->dt;
    }

    int64_t pf;

    int64_t halfN=(N>>1)+1;

    double complex *Xbuff = (double complex *)malloc(sizeof(double complex)*N);
    double *winFunc = (double *)malloc(sizeof(double)*N);
    double *xmode = (double *)malloc(sizeof(double)*N);
    double complex *nufft_coef = (double complex *)malloc(sizeof(double complex)*N);

    //  create window function according to the winOption parsed by user.
    nufft_winFunc_create(winFunc, N, t, &optParse.winOption);

    double c_average=0;
    double t_duration=t[N-1]-t[0];
    for(pf=0;pf<N;++pf){
        c_average+=X[pf];
        Xbuff[pf]=X[pf];
        xmode[pf]=2*pi*((t[pf]-t[0])/t_duration-0.5);
    }
    c_average/=N;
    for(pf=0;pf<N;++pf){
        Xbuff[pf]=(Xbuff[pf]-c_average)*winFunc[pf];
    }

    finufft_opts finfftOpts;                      // make an opts struct
    finufft_default_opts(&finfftOpts);          // set default opts (must do this)
    finfftOpts.modeord = 1;
    finfftOpts.nthreads = 1;
    finfftOpts.debug=optParse.debug;
    finufft1d1(N,xmode,Xbuff,+1,optParse.tolerance,N,nufft_coef,&finfftOpts);
    for(pf=0;pf<halfN;++pf){
        Result->freq[pf]=pf/t_duration;
        Result->ASD[pf]=cabs(nufft_coef[pf])*sqrt(2*t_duration)/(N-1);
    }
    Result->N=halfN;
    if(optParse.smoothMode){
        ASD_filter(Result,optParse.smoothNodes);
    }

    free(Xbuff);free(winFunc);free(xmode);free(nufft_coef);
    
    return 0;
}

int ASD_H_U2N(fftOutput_t *Result, int64_t Nf, double *freq, int64_t Nx, double *X, fftOption_t *opts){
    fftOption_t optParse;
    default_fftOption_init(&optParse);
    if(opts==NULL){
        //Don't know what to do :D
    }else{
        optParse.winOption.winFuncName=opts->winOption.winFuncName;
        optParse.winOption.alpha=opts->winOption.alpha;
        optParse.winLengthRatio=opts->winLengthRatio;
        optParse.overlapRatio=opts->overlapRatio;
        optParse.smoothMode=opts->smoothMode;
        optParse.smoothNodes=opts->smoothNodes;
        optParse.debug=opts->debug;
        optParse.tolerance=opts->tolerance;
        optParse.dt=opts->dt;
    }

    int64_t pf, autoTick=0;

    double complex *Xbuff = (double complex *)malloc(sizeof(double complex)*Nx);
    double *winFunc = (double *)malloc(sizeof(double)*Nx);
    double *xmode = (double *)malloc(sizeof(double)*Nf);
    double complex *nufft_coef = (double complex *)malloc(sizeof(double complex)*Nf);

    //  create window function according to the winOption parsed by user.
    nufft_winFunc_create(winFunc, Nx, NULL, &optParse.winOption);

    double c_average=0;
    double t_duration=(Nx-1)*optParse.dt;

    for(pf=0;pf<Nx;++pf){
        c_average+=X[pf];
        Xbuff[pf]=X[pf];
    }
    c_average/=Nx;
    for(pf=0;pf<Nx;++pf){
        Xbuff[pf]=(Xbuff[pf]-c_average)*winFunc[pf];
    }
    for(pf=0;pf<Nf;++pf){
        xmode[pf] = 2*pi*freq[pf]*optParse.dt;
    }

    finufft_opts finfftOpts;                      // make an opts struct
    finufft_default_opts(&finfftOpts);          // set default opts (must do this)
    finfftOpts.modeord = 0;
    finfftOpts.nthreads = 1;
    finfftOpts.debug=optParse.debug;
    
    finufft1d2(Nf,xmode,nufft_coef,+1,optParse.tolerance,Nx,Xbuff,&finfftOpts);

    for(pf=0;pf<Nf;++pf){
        Result->freq[pf]=freq[pf];
        Result->ASD[pf]=cabs(nufft_coef[pf])*sqrt(2*t_duration)/(Nx-1);
    }
    Result->N=Nf;

    if(optParse.smoothMode){
        ASD_filter(Result,optParse.smoothNodes);
    }

    free(Xbuff);free(winFunc);free(xmode);free(nufft_coef);
    
    return 0;
}

int ASD_H_N2N(fftOutput_t *Result, int64_t Nf, double *freq, int64_t Nx, double *X, double *t, fftOption_t *opts){
    fftOption_t optParse;
    default_fftOption_init(&optParse);
    if(opts==NULL){
        //Don't know what to do :D
    }else{
        optParse.winOption.winFuncName=opts->winOption.winFuncName;
        optParse.winOption.alpha=opts->winOption.alpha;
        optParse.winLengthRatio=opts->winLengthRatio;
        optParse.overlapRatio=opts->overlapRatio;
        optParse.smoothMode=opts->smoothMode;
        optParse.smoothNodes=opts->smoothNodes;
        optParse.debug=opts->debug;
        optParse.tolerance=opts->tolerance;
        optParse.dt=opts->dt;
    }

    int64_t pf;

    double complex *Xbuff = (double complex *)malloc(sizeof(double complex)*Nx);
    double *xmode = (double *)malloc(sizeof(double)*Nx);
    double *winFunc = (double *)malloc(sizeof(double)*Nx);
    double *smode = (double *)malloc(sizeof(double)*Nf);
    double complex *nufft_coef = (double complex *)malloc(sizeof(double complex)*Nf);

    //  create window function according to the winOption parsed by user.
    nufft_winFunc_create(winFunc, Nx, t, &optParse.winOption);

    double c_average=0;
    double t_duration=t[Nx-1]-t[0];
    for(pf=0;pf<Nx;++pf){
        c_average+=X[pf];
        Xbuff[pf]=X[pf];
    }
    c_average/=Nx;
    for(pf=0;pf<Nx;++pf){
        Xbuff[pf]=(Xbuff[pf]-c_average)*winFunc[pf];
        xmode[pf]=t[pf]-t[0];
    }
    for(pf=0;pf<Nf;++pf){
        smode[pf] = 2*pi*freq[pf];
    }

    finufft_opts finfftOpts;                      // make an opts struct
    finufft_default_opts(&finfftOpts);          // set default opts (must do this)
    finfftOpts.modeord = 0;
    finfftOpts.nthreads = 1;
    finfftOpts.debug=optParse.debug;
    
    finufft1d3(Nx,xmode,Xbuff,+1,optParse.tolerance,Nf,smode,nufft_coef,&finfftOpts);

    for(pf=0;pf<Nf;++pf){
        Result->freq[pf]=freq[pf];
        Result->ASD[pf]=cabs(nufft_coef[pf])*sqrt(2*t_duration)/(Nx-1);
    }
    Result->N=Nf;

    if(optParse.smoothMode){
        ASD_filter(Result,optParse.smoothNodes);
    }

    free(Xbuff);free(xmode);free(winFunc);free(smode);free(nufft_coef);
    
    return 0;
}

int ASD_H_N2UviaFile(char outFilePath[], char xFilePath[], char tFilePath[], fftOption_t *opts){
    FILE *outf, *xfile, *tfile, *outfSp;
    char pathSp[256];
    if((xfile=fopen(xFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",xFilePath);
        return -1;
    }
    if((tfile=fopen(tFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",tFilePath);
        return -1;
    }
    if((outf=fopen(outFilePath,"w"))==NULL){
        fprintf(stderr,"unable to write %s\n",outFilePath);
        return -1;
    }
    struct stat statBuff;
    int64_t sizeXfile, sizeTfile, estLine, N, hN, pf;

    if(stat(xFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(xFile).\n");
        return -1;
    }
    sizeXfile=statBuff.st_size;

    if(stat(tFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(tFile).\n");
        return -1;
    }
    sizeTfile=statBuff.st_size;

    // estimate how many lines are required.
    estLine=(sizeXfile>sizeTfile?sizeXfile:sizeTfile)>>2;

    fftOutput_t Result;

    double* X = (double *)malloc(sizeof(double)*estLine);
    double* t = (double *)malloc(sizeof(double)*estLine);

    // read data from file
    for(N=0; N<estLine; ++N){
        if(fscanf(tfile,"%lf%*[^0-9+-]",&t[N])<1) break;
        if(fscanf(xfile,"%lf%*[^0-9+-]",&X[N])<1) break;
    }
    fclose(xfile); fclose(tfile);

    hN=(N>>1)+1;
    Result.ASD=(double *)malloc(sizeof(double)*hN);
    Result.freq=(double *)malloc(sizeof(double)*hN);

    int N_stack, fltFlag=0;
    if(opts){
        N_stack=opts->smoothNodes;
        fltFlag=opts->smoothMode;
    }else{
        N_stack=256;
        fltFlag=0;
    }
    Result.fltFreq=(double *)malloc(sizeof(double)*32*N_stack);
    Result.fltASD=(double *)malloc(sizeof(double)*32*N_stack);

    ASD_H_N2U(&Result, N, X, t, opts);

    if(fltFlag==1){
        for(pf=0;pf<Result.fltN;++pf){
            fprintf(outf,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
        }
    }else if(fltFlag==2){
        sprintf(pathSp,"%s_flt.csv",outFilePath);
        if((outfSp=fopen(pathSp,"w"))==NULL){
            fprintf(stderr,"unable to write %s\n",pathSp);
        }else{
            fprintf(outfSp,"freq,ASD\n");
            for(pf=0;pf<Result.fltN;++pf){
                fprintf(outfSp,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
            }
            fclose(outfSp);
        }
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }else{
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }
    
    fclose(outf);

    free(X); free(t); free(Result.freq); free(Result.ASD); free(Result.fltASD); free(Result.fltFreq);

    return 0;
}

int ASD_H_U2NviaFile(char outFilePath[], char freqFilePath[], char xFilePath[], fftOption_t *opts){
    FILE *outf, *xfile, *freqfile, *outfSp;
    char pathSp[256];
    if((xfile=fopen(xFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",xFilePath);
        return -1;
    }
    if((freqfile=fopen(freqFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",freqFilePath);
        return -1;
    }
    if((outf=fopen(outFilePath,"w"))==NULL){
        fprintf(stderr,"unable to write %s\n",outFilePath);
        return -1;
    }
    struct stat statBuff;
    int64_t sizeXfile, sizeFreqfile, estLineX, estLineFreq, Nx, Nf, pf;

    if(stat(xFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(xFile).\n");
        return -1;
    }
    sizeXfile=statBuff.st_size;

    if(stat(freqFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(tFile).\n");
        return -1;
    }
    sizeFreqfile=statBuff.st_size;

    // estimate how many lines are required.
    estLineX=sizeXfile>>2;
    estLineFreq=sizeFreqfile>>2;

    fftOutput_t Result;
    double* X = (double *)malloc(sizeof(double)*estLineX);
    double* freq = (double *)malloc(sizeof(double)*estLineFreq);

    // read X data from file
    for(Nx=0; Nx<estLineX; ++Nx){
        if(fscanf(xfile,"%lf%*[^0-9+-]",&X[Nx])<1) break;
    }
    for(Nf=0; Nf<estLineFreq; ++Nf){
        if(fscanf(freqfile,"%lf%*[^0-9+-]",&freq[Nf])<1) break;
    }
    fclose(xfile); fclose(freqfile);

    Result.ASD=(double *)malloc(sizeof(double)*Nf);
    Result.freq=(double *)malloc(sizeof(double)*Nf);

    int N_stack, fltFlag=0;
    if(opts){
        N_stack=opts->smoothNodes;
        fltFlag=opts->smoothMode;
    }else{
        N_stack=256;
        fltFlag=0;
    }
    Result.fltFreq=(double *)malloc(sizeof(double)*32*N_stack);
    Result.fltASD=(double *)malloc(sizeof(double)*32*N_stack);

    ASD_H_U2N(&Result, Nf, freq, Nx, X, opts);

    if(fltFlag==1){
        for(pf=0;pf<Result.fltN;++pf){
            fprintf(outf,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
        }
    }else if(fltFlag==2){
        sprintf(pathSp,"%s_flt.csv",outFilePath);
        if((outfSp=fopen(pathSp,"w"))==NULL){
            fprintf(stderr,"unable to write %s\n",pathSp);
        }else{
            fprintf(outfSp,"freq,ASD\n");
            for(pf=0;pf<Result.fltN;++pf){
                fprintf(outfSp,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
            }
            fclose(outfSp);
        }
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }else{
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }

    fclose(outf);
    
    free(X); free(freq); free(Result.freq); free(Result.ASD); free(Result.fltFreq); free(Result.fltASD);

    return 0;
}

int ASD_H_N2NviaFile(char outFilePath[], char freqFilePath[], char xFilePath[], char tFilePath[], fftOption_t *opts){
    FILE *outf, *xfile, *tfile, *freqfile, *outfSp;
    char pathSp[256];
    if((xfile=fopen(xFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",xFilePath);
        return -1;
    }
    if((tfile=fopen(tFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",tFilePath);
        return -1;
    }
    if((freqfile=fopen(freqFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",freqFilePath);
        return -1;
    }
    if((outf=fopen(outFilePath,"w"))==NULL){
        fprintf(stderr,"unable to write %s\n",outFilePath);
        return -1;
    }
    struct stat statBuff;
    int64_t sizeXfile, sizeTfile, sizeFreqfile, estLineX, estLineFreq, Nx, Nf, pf;

    if(stat(xFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(xFile).\n");
        return -1;
    }
    sizeXfile=statBuff.st_size;

    if(stat(tFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(tFile).\n");
        return -1;
    }
    sizeTfile=statBuff.st_size;

    if(stat(freqFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(tFile).\n");
        return -1;
    }
    sizeFreqfile=statBuff.st_size;

    // estimate how many lines are required.
    estLineX=(sizeXfile>sizeTfile?sizeXfile:sizeTfile)>>2;
    estLineFreq=sizeFreqfile>>2;

    fftOutput_t Result;
    
    double* X = (double *)malloc(sizeof(double)*estLineX);
    double* t = (double *)malloc(sizeof(double)*estLineX);
    double* freq = (double *)malloc(sizeof(double)*estLineFreq);

    // read X data from file
    for(Nx=0; Nx<estLineX; ++Nx){
        if(fscanf(xfile,"%lf%*[^0-9+-]",&X[Nx])<1) break;
        if(fscanf(tfile,"%lf%*[^0-9+-]",&t[Nx])<1) break;
    }
    for(Nf=0; Nf<estLineFreq; ++Nf){
        if(fscanf(freqfile,"%lf%*[^0-9+-]",&freq[Nf])<1) break;
    }
    fclose(xfile); fclose(tfile); fclose(freqfile);

    Result.ASD=(double *)malloc(sizeof(double)*Nf);
    Result.freq=(double *)malloc(sizeof(double)*Nf);

    int N_stack, fltFlag=0;
    if(opts){
        N_stack=opts->smoothNodes;
        fltFlag=opts->smoothMode;
    }else{
        N_stack=256;
        fltFlag=0;
    }
    Result.fltFreq=(double *)malloc(sizeof(double)*32*N_stack);
    Result.fltASD=(double *)malloc(sizeof(double)*32*N_stack);

    ASD_H_N2N(&Result, Nf, freq, Nx, X, t, opts);

    if(fltFlag==1){
        for(pf=0;pf<Result.fltN;++pf){
            fprintf(outf,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
        }
    }else if(fltFlag==2){
        sprintf(pathSp,"%s_flt.csv",outFilePath);
        if((outfSp=fopen(pathSp,"w"))==NULL){
            fprintf(stderr,"unable to write %s\n",pathSp);
        }else{
            fprintf(outfSp,"freq,ASD\n");
            for(pf=0;pf<Result.fltN;++pf){
                fprintf(outfSp,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
            }
            fclose(outfSp);
        }
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }else{
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }

    fclose(outf);
    
    free(X); free(t); free(freq); free(Result.ASD); free(Result.freq); free(Result.fltFreq); free(Result.fltASD);

    return 0;
}

int ASD_H_U2UviaFile(char outFilePath[], char xFilePath[], fftOption_t *opts){
    FILE *outf, *xfile, *outfSp;
    char pathSp[256];
    if((xfile=fopen(xFilePath,"r"))==NULL){
        fprintf(stderr,"unable to read %s\n",xFilePath);
        return -1;
    }
    if((outf=fopen(outFilePath,"w"))==NULL){
        fprintf(stderr,"unable to write %s\n",outFilePath);
        return -1;
    }
    struct stat statBuff;
    int64_t sizeXfile, estLineX, Nx, pf;

    if(stat(xFilePath,&statBuff)<0){
        fprintf(stderr,"error when stat(xFile).\n");
        return -1;
    }
    sizeXfile=statBuff.st_size;

    // estimate how many lines are required.
    estLineX=sizeXfile>>2;

    fftOutput_t Result;
    double* X = (double *)malloc(sizeof(double)*estLineX);
    // read X data from file
    for(Nx=0; Nx<estLineX; ++Nx){
        if(fscanf(xfile,"%lf%*[^0-9+-]",&X[Nx])<1) break;
    }
    fclose(xfile);

    int64_t Nf=(Nx>>1);
    double t_duration=opts==NULL?(Nx-1)*0.01:(Nx-1)*(opts->dt);
    double *freq=(double *)malloc(sizeof(double)*Nf);
    for(pf=0;pf<Nf;++pf){
        freq[pf]=pf/t_duration;
    }

    Result.ASD=(double *)malloc(sizeof(double)*Nf);
    Result.freq=(double *)malloc(sizeof(double)*Nf);

    int N_stack, fltFlag=0;
    if(opts){
        N_stack=opts->smoothNodes;
        fltFlag=opts->smoothMode;
    }else{
        N_stack=256;
        fltFlag=0;
    }
    Result.fltFreq=(double *)malloc(sizeof(double)*32*N_stack);
    Result.fltASD=(double *)malloc(sizeof(double)*32*N_stack);

    ASD_H_U2N(&Result, Nf, freq, Nx, X, opts);

    if(fltFlag==1){
        for(pf=0;pf<Result.fltN;++pf){
            fprintf(outf,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
        }
    }else if(fltFlag==2){
        sprintf(pathSp,"%s_flt.csv",outFilePath);
        if((outfSp=fopen(pathSp,"w"))==NULL){
            fprintf(stderr,"unable to write %s\n",pathSp);
        }else{
            fprintf(outfSp,"freq,ASD\n");
            for(pf=0;pf<Result.fltN;++pf){
                fprintf(outfSp,"%e,%e\n",Result.fltFreq[pf],Result.fltASD[pf]);
            }
            fclose(outfSp);
        }
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }else{
        for(pf=1;pf<Result.N;++pf){
            fprintf(outf,"%e,%e\n",Result.freq[pf],Result.ASD[pf]);
        }
    }
    
    fclose(outf);
    
    free(X); free(Result.ASD); free(freq); free(Result.freq); free(Result.fltASD); free(Result.fltFreq);

    return 0;
}

int ASD_H_auto(char outFilePath[], char freqFilePath[], char xFilePath[], char tFilePath[], fftOption_t *opts){
    int status;
    fftOption_t optsDefault;
    default_fftOption_init(&optsDefault);
    if(opts==NULL){
        opts=&optsDefault;
    }
    if(freqFilePath==NULL){
        if(tFilePath==NULL){
            status=ASD_H_U2UviaFile(outFilePath, xFilePath, opts);
            goto autoExit;
        }else{
            status=ASD_H_N2UviaFile(outFilePath, xFilePath, tFilePath, opts);
            goto autoExit;
        }
    }else{
        if(tFilePath==NULL){
            status=ASD_H_U2NviaFile(outFilePath, freqFilePath, xFilePath, opts);
            goto autoExit;
        }else{
            status=ASD_H_N2NviaFile(outFilePath, freqFilePath, xFilePath, tFilePath, opts);
            goto autoExit;
        }
    }

    autoExit:
    return status;
}