#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* data;
double pi = 3.14159265358979323846;
double pi2 = 6.28318530717958647692;

#define REAL 0
#define IMAG 1
//data[point][real or imag]

//returns answer in data[]
void cmult(double* data, double* data2, double* output){
    output[REAL] = (data[REAL] * data2[REAL]) - (data[IMAG] * data2[IMAG]);
    output[IMAG] = (data[REAL] * data2[IMAG]) + (data[IMAG] * data2[REAL]);
}

void DTFT(double** x, double** X, int N){
    for(int k = 0; k < N; k++){
        X[k][REAL] = 0;
        X[k][IMAG] = 0;
        for(int n = 0; n < N; n++){
            double val[3][2];
            val[0][REAL] = x[n][0];
            val[0][IMAG] = x[n][1];
            val[1][REAL] = cos(pi2*k*n/N);
            val[1][IMAG] = -1*sin(pi2*k*n/N);
            cmult(val[0], val[1], val[2]);
            X[k][REAL]= X[k][REAL] + val[2][REAL];
            X[k][IMAG]= X[k][IMAG] + val[2][IMAG];
        }
    }
}

void cooley_tukey(double** xData, double** X, int N, int columns, int rows){
    if((columns*rows) != N){
        fprintf(stderr, "Invalid matrix dimensions (%d x %d != %d).\n", columns, rows, N);
        exit(1);
    }
    /* Reshape */
    double*** x2 = (double***)malloc(sizeof(double**) * rows);
    double*** x2_0 = (double***)malloc(sizeof(double**) * rows);
    for(int x = 0; x < rows; x++){
        x2[x] = (double**)malloc(sizeof(double*) * columns);
        x2_0[x] = (double**)malloc(sizeof(double*) * columns);
        for(int y = 0; y < columns; y++){
            x2[x][y] = (double*)malloc(sizeof(double) * 2);
            x2_0[x][y] = (double*)malloc(sizeof(double) * 2);
            x2[x][y][REAL] = xData[rows*y+x][REAL];
            x2[x][y][IMAG] = xData[rows*y+x][IMAG];
        }
    }
    double*** x2_1 = (double***)malloc(sizeof(double**) * columns);
    double*** x2_2 = (double***)malloc(sizeof(double**) * columns);
    for(int y = 0; y < columns; y++){
        x2_1[y] = (double**)malloc(sizeof(double*) * rows);
        x2_2[y] = (double**)malloc(sizeof(double*) * rows);
        for(int x = 0; x < rows; x++){
            x2_1[y][x] = (double*)malloc(sizeof(double) * 2);
            x2_2[y][x] = (double*)malloc(sizeof(double) * 2);
        }
    }
    /* DTFT rows */
    for(int x = 0; x < rows; x++){
        DTFT(x2[x], x2_0[x], columns);
    }
    /* Twiddle Factors and rearrange for DTFT by row*/
    for(int x = 0; x < rows; x++){
        for(int y = 0; y < columns; y++){
            double val[2];
            val[0] = cos(pi2*x*y/N);
            val[1] = -1*sin(pi2*x*y/N);
            cmult(x2_0[x][y], val, x2_1[y][x]);
        }
    }
    /* DTFT columns */
    for(int y = 0; y < columns; y++){
        DTFT(x2_1[y], x2_2[y], rows);
    }
    /* Reshape into 1 dimensional array */
    for(int pos = 0; pos < N; pos++){
        int x = pos%rows;
        int y = (pos-x)/rows;
        X[pos][0] = x2_2[y][x][0];
        X[pos][1] = x2_2[y][x][1];
    }

    /* Free Memory */
    for(int x = 0; x < rows; x++){
        for(int y = 0; y < columns; y++){
            free(x2[x][y]);
            free(x2_0[x][y]);
        }
        free(x2[x]);
        free(x2_0[x]);
    }
    for(int y = 0; y < columns; y++){
        for(int x = 0; x < rows; x++){
            free(x2_1[y][x]);
            free(x2_2[y][x]);
        }
        free(x2_1[y]);
        free(x2_2[y]);
    }
    free(x2_1);
    free(x2_2);
    free(x2);
    free(x2_0);
}

void absolute(double** data, double* absoluteData, int N){
    for(int i = 0; i < N; i++){
        absoluteData[i] = sqrt(pow(data[i][0],2)+pow(data[i][1],2));
    }
}

void real2complex(double* real, double** complex, int N){
    for(int i = 0; i < N; i++){
        complex[i][0] = real[i];
        complex[i][1] = 0.0;
    }
}

int main(int argc, char** argv){
    int N = 35; // number of points
    double Fs = 1000; // sampling frequency
    int Fsignal = 90; // signalnal frequency
    double** X_DTFT = (double**)malloc(sizeof(double*) * N);
    double** X_FFT = (double**)malloc(sizeof(double*) * N);
    double** sig = (double**)malloc(sizeof(double*) * N);
    double* X = (double*)malloc(sizeof(double) * N);
    for(int i = 0; i < N; i++){
        X_DTFT[i] = (double*)malloc(sizeof(double) * 2);
        X_FFT[i] = (double*)malloc(sizeof(double) * 2);
        sig[i] = (double*)malloc(sizeof(double) * 2);
    }
    fprintf(stdout, "Signal Data:\n");
    for(int i = 0; i < N; i++){
        sig[i][0] = sin(2*pi*(double)Fsignal*i/(double)Fs);
        sig[i][1] = 0.0;
        fprintf(stdout, "%f\n", sig[i][0]);
    }
    DTFT(sig, X_DTFT, N);
    absolute(X_DTFT, X, N);
    fprintf(stdout, "\nDTFT:\n");
    for(int i = 0; i < N; i++){
        fprintf(stdout, "%f\n", X[i]);
    }
    printf("\n Cooley-Tukey\n\n");
    for(int i = 0; i < N; i++){
        sig[i][0] = sin(2*pi*(double)Fsignal*i/(double)Fs);
        sig[i][1] = 0.0;
    }
    cooley_tukey(sig, X_FFT, 35, 7, 5);
    absolute(X_FFT, X, N);
    for(int i = 0; i < N; i++){
        fprintf(stdout, "%f\n", X[i]);
    }
    return 0;
}
