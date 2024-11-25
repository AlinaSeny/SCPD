#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#define  Max(a,b) ((a)>(b)?(a):(b))
#define  Min(a,b) ((a)<(b)?(a):(b))

float maxeps = 0.1e-7;
int itmax = 100;
int i,j,k;
int num_of_threads;

float eps;
float A [N][N];

void relax();
void init();
void verify();

int main(int argc, char **argv)
{
    if (argc < 2) {
        return 1;
    }

    num_of_threads = (int)strtol(argv[1], NULL, 10);
    omp_set_num_threads(num_of_threads);
    double avg_time = 0.;

    for (int rep_num = 0; rep_num < REPETITIONS; rep_num++) {
#ifdef DEBUG
        printf("Started repetition: %d\n", rep_num);
#endif

        int it;
        init();

        double start = omp_get_wtime();

        for (it = 1; it <= itmax; it++) {
            eps = 0;
            relax();
#ifdef DEBUG
            if (!(it % 1000)) {
                printf("it=%4i   eps=%6f\n", it, eps);
            }
#endif
            if (eps < maxeps) break;
        }

        double finish = omp_get_wtime();

#ifdef DEBUG
        printf("Iterations: %d\n", it);
        printf("Time: %fs\n", finish - start);
#endif

        avg_time += finish - start;
        verify();
    }
    avg_time /= REPETITIONS;
    printf("%lf\n", avg_time);
    return 0;
}

void init()
{
    for(i=0; i<=N-1; i++) {
        for (j = 0; j <= N - 1; j++) {
            if (i == 0 || i == N - 1 || j == 0 || j == N - 1) A[i][j] = 0.;
            else A[i][j] = (1. + i + j);
        }
    }
}

void relax() {
#pragma omp parallel shared(A, eps, num_of_threads) private(i, j, k) default(none)
    {
        float e;
#pragma omp single
        {
            int per_task = (N - 2) / num_of_threads;
            int num_of_tasks = (N - 2) / per_task;
            if ((N - 2) % per_task != 0) {
                num_of_tasks += 1;
            }
            for (i = 2; i <= N - 3; i++) {
                for (k = 0; k < num_of_tasks; k++) {
#pragma omp task
                    {
                        for (j = k * per_task + 1; j <= Min((k + 1) * per_task, N - 2); j++) {
                            A[i][j] = (A[i - 1][j] + A[i + 1][j] + A[i - 2][j] + A[i + 2][j]) / 4.;
                        }
                    }
                }
#pragma omp taskwait
            }

            for (k = 0; k < num_of_tasks; k++) {
#pragma omp task private(e)
                {
                for (i = k * per_task + 1; i <= Min((k + 1) * per_task, N - 2); i++) {
                    for (j = 2; j <= N - 3; j++) {
                        float tmp = A[i][j];
                        A[i][j] = (A[i][j - 1] + A[i][j + 1] + A[i][j - 2] + A[i][j + 2]) / 4.;
                        e = Max(e, fabs(tmp - A[i][j]));
                    }
                }
#pragma omp critical
                    eps = Max(eps, e);
                }
            }
        }
        }
    }

void verify()
{
    float s;

    s=0.;
    for(i=0; i<=N-1; i++) {
        for (j = 0; j <= N - 1; j++) {
            s = s + A[i][j] * (i + 1) * (j + 1) / (N * N);
        }
    }
#ifdef DEBUG
    printf("\tS = %f\n", s);
#endif
}
