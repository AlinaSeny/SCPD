#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#define  Max(a,b) ((a)>(b)?(a):(b))
#define  Min(a,b) ((a)<(b)?(a):(b))


float   maxeps = 0.1e-7;
int itmax = 100;
int i,j,k;

float eps, sum;
float A [N][N];

void relax();
void init();
void verify();

int wrank, wsize;
int block, startrow, lastrow;

int main(int argc, char **argv)
{
    double avg_time = 0.;

    //printf("%d\n", argc);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int rep_num = 0; rep_num < REPETITIONS; rep_num++) {
        if (!wrank) {
#ifdef DEBUG
            printf("Started repetition: %d\n", rep_num);
#endif
        }
        block = (N - 2) / wsize;
        startrow = block * wrank + 1;
        if ((wrank == (wsize - 1)) && (N - 2) % wsize != 0) {
            lastrow = startrow + block + (N - 2) % wsize - 1;
        } else {
            lastrow = startrow + block - 1;
        }
        int it;
        init();

        double start, finish;

        if (!wrank) {
            start = MPI_Wtime();
#ifdef DEBUG
            printf("wsize: %d\n", wsize);
#endif
        }

        for (it = 1; it <= itmax; it++) {
            eps = 0.;
            relax();

            if (!wrank) {
#ifdef DEBUG
                if (!(it % 1000)) {
                    printf("it=%4i   eps=%.6f\n", it, eps);
                }
#endif
            }

            if (eps < maxeps) break;
        }

        if (!wrank) {
            finish = MPI_Wtime();
            avg_time += finish - start;
#ifdef DEBUG
            printf("Iterations: %d\n", it);
            printf("Time: %fs\n", finish - start);
#endif
        }

        verify();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if (!wrank) {
        avg_time /= REPETITIONS;
        printf("%lf", avg_time);
    }

    MPI_Finalize();

    return 0;
}

void init()
{
    for(i=0; i<=N-1; i++)
        for(j=0; j<=N-1; j++)
        {
            if(i==0 || i==N-1 || j==0 || j==N-1) A[i][j]= 0.;
            else A[i][j]= ( 1. + i + j ) ;
        }
}

void relax()
{
    for(i=2; i<=N-3; i++)
        for(j=startrow; j<=lastrow; j++)
        {
            A[i][j] = (A[i-1][j]+A[i+1][j]+A[i-2][j]+A[i+2][j])/4.;
        }
    
    for (int w = 0; w < wsize; w++) {
        for (i = block * w + 1; i <= block * (w + 1); i++) {
            MPI_Gather(&A[i][startrow], block, MPI_FLOAT, &A[i][1], block, MPI_FLOAT, w, MPI_COMM_WORLD);
        }
    }

    float local_eps = eps;

    for(i=startrow; i<=lastrow; i++)
        for(j=2; j<=N-3; j++)
        {
            float e=A[i][j];
            A[i][j] = (A[i][j - 1] + A[i][j + 1] + A[i][j - 2] + A[i][j + 2]) / 4.;
            local_eps=Max(local_eps, fabs(e-A[i][j]));
        }

    MPI_Allreduce(&local_eps, &eps, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

    for (int w = 0; w < wsize; w++) {
        for (i = block * w + 1; i <= block * (w + 1); i++) {
            MPI_Scatter(&A[i][1], block, MPI_FLOAT, &A[i][startrow], block, MPI_FLOAT, w, MPI_COMM_WORLD);
        }
    }
}

void verify()
{
    float s=0.;
    for(i=0; i<=N-1; i++)
        for(j=startrow; j<=lastrow; j++)
        {
            s=s+A[i][j]*(i+1)*(j+1)/(N*N);
        }

    MPI_Allreduce(&s, &sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    if (!wrank) {
#ifdef DEBUG
        printf("\tS = %f\n", sum);
#endif
    }

}
