/** plu.c: rewrote lu.c as parallel plu.c program **/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

#define N 8

double A[N][N], L[N][N], U[N][N];
double B[N], b[N];
int    P[N];
double Y[N], X[N];

pthread_barrier_t barrier;

int print(char c, double x[N][N])
{
	int i, j;
	printf("------------- %c -----------------------\n", c);
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
			printf("%6.2f  ", x[i][j]);
		printf("\n");
	}
}

int printV(char c, double x[N])
{
	int i;
	printf("--------- %c vector-----------\n", c);
	for (i = 0; i < N; i++)
		printf("%6.2f ", x[i]);
	printf("\n");
}

int printP()
{
	int i;
	printf("--------- P vector-----------\n");
	for (i = 0; i < N; i++)
		printf("%d ", P[i]);
	printf("\n");
}

// LU decomposition function
int *plu(void* arg)
{
	int i, j, k, m, itemp, prow;
	int myid = (int)arg;
	double max;
	double temp;
	double si;

	// LU decomposition loop
	for (k = 0; k < N; k++) {
		if (k == myid) {
			max = 0;
			printf("partial pivoting by thread %d on row %d: ", myid, k);
			for (i = k; i < N; i++) {
				if (max < fabs(A[i][k])) { // partial pivoting
					max = fabs(A[i][k]);
					j = i;
				}
				prow = i;
			}
			printf("pivot_row= %d pivot=%6.2f\n", prow, A[prow][k]);
			if (max == 0) {
				printf("zero pivot: singular A matrix\n");
				exit(1);
			}

			// swap P[k] and P[j];
			itemp = P[k]; P[k] = P[j]; P[j] = itemp;


			// swap row A[k] and row A[j]
			for (m = 0; m < N; m++) {
				temp = A[k][m]; A[k][m] = A[j][m]; A[j][m] = temp;
			}

			//swap L[k][0,k-2] and L[j][0,k-2]
			for (m = 0; m < k - 2; m++) {
				temp = L[k][m]; L[k][m] = L[j][m]; L[j][m] = temp;
			}
		}

		pthread_barrier_wait(&barrier);
		// compute L U entries
		U[k][k] = A[k][k];

		for (i = k + 1; i < N; i++) {
			L[i][k] = A[i][k] / U[k][k];
			U[k][i] = A[k][i];
		}

		// row reductions on
		for (i = k + 1; i < N; i++) {
			if (i == myid) {
				printf("thread %d do row %d\n", myid, i);
				for (m = k + 1; m < N; m++) {
					A[i][m] -= L[i][k] * U[k][m];
				}
			}
		}

		pthread_barrier_wait(&barrier);
		if (k == myid) {
			print('A', A); print('L', L); print('U', U); printP();
			getchar();
		}
	}
}

int main(int argc, char* argv[])
{
	int i, j, k;
	double si;
	pthread_t threads[N];

	printf("main: initialize matrix A[N][N], B[N], L, U and P\n");
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			A[i][j] = 1.0;

	for (i = 0; i < N; i++)
		A[i][N - 1 - i] = 1.0 * N;

	for (i = 0; i < N; i++) {
		B[i] = (N) * (N + 1) / 2 + (N - i) * (N - 1);
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			U[i][j] = 0.0;
			L[i][j] = 0.0;
			if (i == j)
				L[i][j] = 1.0;
		}
	}

	for (i = 0; i < N; i++) {
		P[i] = i;
	}

	print('A', A);  print('L', L);  print('U', U);
	printV('B', B);	printP();

	pthread_barrier_init(&barrier, NULL, N); //set up barrier

	printf("main: create N=%d working threads\n", N);
	for (i = 0; i < N; i++) {
		pthread_create(&threads[i], NULL, plu, (void*)i);
	}

	printf("main: wait for all %d working threads to join\n", N);
	for (i = 0; i < N; i++) {
		pthread_join(threads[i], NULL);
	}

	// P L U are all computed; solve P*U*L*X = P*B

	printf("main: back substitution : ");
	printf("The solution is :\n");
	// apply P to B to get b[ ]
	printV('B', B);
	for (i = 0; i < N; i++) {
		b[i] = B[P[i]];
	}
	printV('b', b);

	// solve L*Y = PB = b
	for (i = 0; i < N; i++) {  // forwar substitution
		Y[i] = b[i];
		for (j = 0; j < i; j++) {
			Y[i] -= L[i][j] * Y[j];
		}
	}
	printV('Y', Y);

	// solve U*X=Y
	for (i = N - 1; i >= 0; i--) { // backward substitution
		si = 0.0;
		for (j = i + 1; j < N; j++)
			si += U[i][j] * X[j];
		X[i] = (Y[i] - si) / U[i][i];
	}
	printV('X', X);
}