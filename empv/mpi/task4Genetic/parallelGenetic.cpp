#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <mpi.h>

using namespace std;

double frand() // вещественное случайное число в диапазоне [0,1)
{
    return double(rand()) / RAND_MAX;
}

double eval(double* a, int n)
{
    double sum = 0;
    // // sphere
    // for( int i = 0; i < n; i++ )
    //     sum += a[i] * a[i];

    // Rosenbrok
    for( int i = 0; i < n - 1; i++ )
        sum += 100 * (a[i] * a[i] - a[i + 1]) * (a[i] * a[i] - a[i + 1]) + (a[i] - 1) * (a[i] - 1);

    // // Rastring
    // for( int i = 0; i < n; i++ )
    //     sum += a[i] * a[i] - 10 * cos(2 * 3.14159265 * a[i]) + 10;

    return sum;
}

void init(double* P, int m, int n)
{
    for( int k = 0; k < m; k++ )
        for( int i = 0; i < n; i++ )
            P[k * n + i] = frand() * 200 - 100;
}

void shuffle(double* P, int m, int n)
{
    for( int k = 0; k < m; k++ )
    {
        int l = rand() % m;
        for( int i = 0; i < n; i++ )
            swap(P[k * n + i], P[l * n + i]);
    }
}

void select(double* P, int m, int n)
{
    double pwin = 0.75;
    shuffle(P, m, n);
    for( int k = 0; k < m / 2; k++ )
    {
        int a = 2 * k;
        int b = 2 * k + 1;
        double fa = eval(P + a * n, n);
        double fb = eval(P + b * n, n);
        double p = frand();
        if ( ((fa < fb) && (p <= pwin)) || ((fa > fb) && (p >= pwin)) )
            for( int i = 0; i < n; i++ )
                P[b * n + i] = P[a * n + i];
        else
            for( int i = 0; i < n; i++ )
                P[a * n + i] = P[b * n + i];
    }
}

void crossover(double* P, int m, int n)
{
    shuffle(P, m, n);

    // // 1 point
    // for( int k = 0; k < m / 2; k++ )
    // {
    //     int a = 2 * k;
    //     int b = 2 * k + 1;
    //     int j = rand() % n;
    //     for( int i = j; i < n; i++ )
    //         swap(P[a * n + i], P[b * n + i]);
    // }

    // // 2 point
    // for( int k = 0; k < m / 2; k++ )
    // {
    //     int a = 2 * k;
    //     int b = 2 * k + 1;
    //     int j = rand() % n;
    //     int jj = rand() % n;
    //     int tmp = min(j, jj);
    //     jj = max(j, jj);
    //     j = tmp;
    //     for( int i = j; i < jj; i++ )
    //         swap(P[a * n + i], P[b * n + i]);
    // }

    // uniform
    double pcross = 0.5;
    for( int k = 0; k < m / 2; k++ )
    {
        int a = 2 * k;
        int b = 2 * k + 1;
        for( int i = 0; i < n; i++ )
            if (frand() <= pcross)
                swap(P[a * n + i], P[b * n + i]);
    }
}

void mutate(double* P, int m, int n)
{
    double pmut = 0.1;
    for( int k = 0; k < m; k++ )
        for( int i = 0; i < n; i++ )
        {
            double rnd = frand();
            if( rnd <= pmut )
                P[k * n + i] += rnd * 20 - 1;
        }
}

void printthebest(double* P, int m, int n, int myRank, int nProc, int iter)
{
    int k0 = 0;
    double f0 = eval(P, n);
    double m0 = f0;
    for( int k = 1; k < m; k++)
    {
        double f = eval(P + k * n, n);
        m0 += f;
        if( f < f0 )
        {
            f0 = f;
            k0 = k;
        }
    }

    int kg = 0;
    double fg = 0.0;
    double mg = 0.0;

    MPI_Reduce(&m0, &mg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&f0, &fg, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    if (myRank == 0)
        cout << iter << ":   " << "Best error = " << fg << "   " << "Mean = " << mg / (m * nProc) << endl;
}

void migrate(double *P, int m, int n, int s, int left, int right)
{
    shuffle(P, m, n);

    double *tmp = new double [s * n];

    // MPI_Sendrecv(P, s * n, MPI_DOUBLE, left, 0, P, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(P, s * n, MPI_DOUBLE, left, 0, tmp, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    memcpy(P, tmp, s * n * sizeof(double));

    delete [] tmp;
    return;
}

void runGA(int n, int m, int T, int s, int k, int myRank, int nProc)
{
    double* P = new double[n * m];
    init(P, m, n);

    int left = myRank - 1;
    int right = myRank + 1;

    if (left == -1)
        left = nProc - 1;

    if (right == nProc)
        right = 0;

    long double time = MPI_Wtime();
    for( int t = 1; t < T + 1; t++ )
    {
        select(P, m, n);
        crossover(P, m, n);
        mutate(P, m, n);

        if (t % k == 0)
            migrate(P, m, n, s, left, right);

        printthebest(P, m, n, myRank, nProc, t);
    }
    time = MPI_Wtime() - time;
    long double timeR = 0.0;
    MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myRank == 0)
        cout << "> Final time of computation = " << timeR << endl;

    delete[] P;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int n = atoi(argv[1]); // gene size
    int m = atoi(argv[2]); // genome size
    int T = atoi(argv[3]); // iter num
    int s = atoi(argv[4]); // migrate size
    int k = atoi(argv[5]); // migrate iter

    srand(time(NULL) + myRank);

    runGA(n, m, T, s, k, myRank, nProc);

    MPI_Finalize();
    return 0;
}