#ifndef ASD_H
#define ASD_H
#endif

#define NONE -1
#define HANN 0
#define HAMM 1
#define FLAT 2
#define BLACKMAN 3
#define POISSON 4
#define KAISER 5
#define NORM 6
#define GAUSS 6
#define WELCH 7
#define NUTTALL 8
#define PCOS 9
#define PLANCK 10

struct fftWinFunc_option_type
{
    int winFuncName;
    double alpha;
};

typedef struct fftWinFunc_option_type fftWinFuncOpt_t;

struct fftOutput_type
{
    int64_t N;
    double *freq;
    double *ASD;
    double *PSD;
    double *AS;
    double *PS;
    double *phase;
    double *fltFreq;
    double *fltASD;
    int fltN;
};

typedef struct fftOutput_type fftOutput_t;

struct fftOption_type
{
    fftWinFuncOpt_t winOption;
    double dt;
    double winLengthRatio;
    double overlapRatio;
    int smoothMode;
    int smoothNodes;
    int debug;
    double tolerance;
};

typedef struct fftOption_type fftOption_t;

double bessel0(double x);
double bessel1(double x);
double quickmedian(double array[], int64_t length);

int default_fftOption_init(fftOption_t *opts);
int nufft_winFunc_create(double *winFunc, int64_t N, double *t, fftWinFuncOpt_t *opts);
int ASD_filter(fftOutput_t *Result, int N_stack);
int ASD_H_N2U(fftOutput_t *Result, int64_t N, double *X, double *t, fftOption_t *opts);
int ASD_H_U2N(fftOutput_t *Result, int64_t Nf, double *freq, int64_t Nx, double *X, fftOption_t *opts);
int ASD_H_N2N(fftOutput_t *Result, int64_t Nf, double *freq, int64_t Nx, double *X, double *t, fftOption_t *opts);
int ASD_H_N2UviaFile(char outFilePath[], char xFilePath[], char tFilePath[], fftOption_t *opts);
int ASD_H_U2NviaFile(char outFilePath[], char freqFilePath[], char xFilePath[], fftOption_t *opts);
int ASD_H_N2NviaFile(char outFilePath[], char freqFilePath[], char xFilePath[], char tFilePath[], fftOption_t *opts);
int ASD_H_U2UviaFile(char outFilePath[], char xFilePath[], fftOption_t *opts);
int ASD_H_auto(char outFilePath[], char freqFilePath[], char xFilePath[], char tFilePath[], fftOption_t *opts);