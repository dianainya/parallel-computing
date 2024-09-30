#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

#define SCHEDULE

int main(int argc, char *argv[]) {
    int i, N, THREADS;
    struct timeval T1, T2;

    N = atoi(argv[1]);
    THREADS = atoi(argv[2]);

    #ifdef _OPENMP
    omp_set_num_threads(THREADS);
    #endif

    gettimeofday(&T1, NULL);
    int N2 = N / 2; // N2 = N/2
    double *restrict M1 = (double *) malloc(N * sizeof(double));
    double *restrict M2 = (double *) malloc(N2 * sizeof(double));
    double *restrict M2_old = (double *) malloc(N2 * sizeof(double));

    double A = 455.0;
    double min = 1.0;
    double max = A;
    double max_2 = max * 10;
    double key = 0.0, X = 0.0;  // Initialize key and X
    long delta_ms;
    int j, k;

    unsigned seeds[THREADS];
    for (i = 0; i < THREADS; i++) {
        seeds[i] = i;
    }

    for (i = 0; i < 100; ++i) {
        X = 0.0;
        key = 0.0;

        #pragma omp parallel default(none) shared(N, N2, M1, M2, M2_old, seeds, min, max, max_2, key, X) private(j, k)
        {
            // GENERATE
            #pragma omp for nowait
            for (j = 0; j < N; ++j) {
                M1[j] = ((double) rand_r(seeds + omp_get_thread_num()) / (RAND_MAX)) * (max - min) + min;
            }
            #pragma omp for
            for (k = 0; k < N2; ++k) {
                M2[k] = ((double) rand_r(seeds + omp_get_thread_num()) / (RAND_MAX)) * (max_2 - max) + max;
            }

            // MAP
            #if defined(SCHEDULE)
            #pragma omp for schedule(runtime) nowait
            #else
            #pragma omp for nowait
            #endif
            for (j = 0; j < N; ++j) {
                M1[j] = pow(M1[j] / M_PI, 3);
            }

            #if defined(SCHEDULE)
            #pragma omp for schedule(runtime) nowait
            #else
            #pragma omp for
            #endif
            for (k = 0; k < N2; ++k) {
                M2_old[k] = M2[k];
            }

            #if defined(SCHEDULE)
            #pragma omp for schedule(runtime)
            #else
            #pragma omp for
            #endif
            for (k = 1; k < N2; ++k) {
                M2[k] = M2[k] + M2_old[k - 1];
            }

            #if defined(SCHEDULE)
            #pragma omp for schedule(runtime)
            #else
            #pragma omp for
            #endif
            for (k = 0; k < N2; ++k) {
                M2[k] = fabs(sin(M2[k]));
            }

            // MERGE
            #if defined(SCHEDULE)
            #pragma omp for schedule(runtime)
            #else
            #pragma omp for
            #endif
            for (k = 0; k < N2; ++k) {
                M2[k] = M1[k] * M2[k];
            }

            // REDUCE
            #pragma omp for schedule(static) reduction(min:key)
            for (k = 1; k < N2; ++k) {
                if (M2[k] != 0 && (key == 0 || M2[k] < key)) {
                    key = M2[k];
                }
            }

            #if defined(SCHEDULE)
            #pragma omp for reduction(+ : X) schedule(runtime)
            #else
            #pragma omp for reduction(+ : X)
            #endif
            for (k = 0; k < N2; ++k) {
                if (((int) (M2[k] / key) % 2) == 0) {
                    X += sin(M2[k]);
                }
            }
        }
    }

//    printf("%f\n", X);

    gettimeofday(&T2, NULL);
    delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
    printf("%ld\n", delta_ms);

    free(M1);
    free(M2);
    free(M2_old);

    return 0;
}
