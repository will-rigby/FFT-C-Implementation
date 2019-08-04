#include <stdlib.h>
#include <math.h>

float* data;
float pi = 3.14159265358979323846;
float pi2 = 6.28318530717958647692;

#define REAL 0
#define IMAG 1
//data[point][real or imag]

//returns answer in data[]
void cmult(float* data, float* data2, float* output){
    output[REAL] = (data[REAL] * data2[REAL]) - (data[IMAG] * data2[IMAG]);
    output[IMAG] = (data[REAL] * data2[IMAG]) + (data[IMAG] * data2[REAL]);
}

void DTFT(float** x, float** X, int N){
    for(int k = 0; k < N; k++){
        X[k][REAL] = 0;
        X[k][IMAG] = 0;
        for(int n = 0; n < N; n++){
            float val[3][2];
            val[0][REAL] = x[n][0];
            val[0][IMAG] = x[n][1];
            val[1][REAL] = cosf(pi2*k*n/N);
            val[1][IMAG] = -1*sinf(pi2*k*n/N);
            cmult(val[0], val[1], val[2]);
            X[k][REAL]= X[k][REAL] + val[2][REAL];
            X[k][IMAG]= X[k][IMAG] + val[2][IMAG];
        }
    }
}

void cooley_tukey(float** xData, float** X, int N, int columns, int rows){
    /* Reshape */
    float*** x2 = (float***)malloc(sizeof(float**) * rows);
    float*** x2_0 = (float***)malloc(sizeof(float**) * rows);
    for(int x = 0; x < rows; x++){
        x2[x] = (float**)malloc(sizeof(float*) * columns);
        x2_0[x] = (float**)malloc(sizeof(float*) * columns);
        for(int y = 0; y < columns; y++){
            x2[x][y] = (float*)malloc(sizeof(float) * 2);
            x2_0[x][y] = (float*)malloc(sizeof(float) * 2);
            x2[x][y][REAL] = xData[rows*y+x][REAL];
            x2[x][y][IMAG] = xData[rows*y+x][IMAG];
        }
    }
    float*** x2_1 = (float***)malloc(sizeof(float**) * columns);
    float*** x2_2 = (float***)malloc(sizeof(float**) * columns);
    for(int y = 0; y < columns; y++){
        x2_1[y] = (float**)malloc(sizeof(float*) * rows);
        x2_2[y] = (float**)malloc(sizeof(float*) * rows);
        for(int x = 0; x < rows; x++){
            x2_1[y][x] = (float*)malloc(sizeof(float) * 2);
            x2_2[y][x] = (float*)malloc(sizeof(float) * 2);
        }
    }
    /* DTFT rows */
    for(int x = 0; x < rows; x++){
        DTFT(x2[x], x2_0[x], columns);
    }
    /* Twiddle Factors and rearrange for DTFT by row*/
    for(int x = 0; x < rows; x++){
        for(int y = 0; y < columns; y++){
            float val[2];
            val[0] = cosf(pi2*x*y/N);
            val[1] = -1*sinf(pi2*x*y/N);
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

void absolute(float** data, float* absoluteData, int N){
    for(int i = 0; i < N; i++){
        absoluteData[i] = sqrtf(powf(data[i][0],2)+powf(data[i][1],2));
    }
}

void real2complex(float* real, float** complex, int N){
    for(int i = 0; i < N; i++){
        complex[i][0] = real[i];
        complex[i][1] = 0.0;
    }
}

int main(int argc, char** argv){
    int N = 35; // number of points
    float Fs = 1000; // sampling frequency
    int Fsignal = 90; // signalnal frequency
    float** X_DTFT = (float**)malloc(sizeof(float*) * N);
    float** sig = (float**)malloc(sizeof(float*) * N);
    float* X = (float*)malloc(sizeof(float) * N);
    for(int i = 0; i < N; i++){
        X_DTFT[i] = (float*)malloc(sizeof(float) * 2);
        sig[i] = (float*)malloc(sizeof(float) * 2);
    }
    for(int i = 0; i < N; i++){
        sig[i][0] = sinf(2*pi*(float)Fsignal*i/(float)Fs);
        sig[i][1] = 0.0;    }
    DTFT(sig, X_DTFT, N);
    absolute(X_DTFT, X, N);
    for(int i = 0; i < N; i++){
        free(X_DTFT[i]);
        free(sig[i]);
    }
    free(X_DTFT);
    free(sig);
    free(X);
    return 0;
}
