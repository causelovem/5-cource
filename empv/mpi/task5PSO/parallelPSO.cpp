#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <mpi.h>

using namespace std;

double a = -100.0; // initial position limits 
double b = 100.0;  // initial position limits
double tau = 1.0;  // delta time
double eta = 0.1;
double v0 = 10.0;   // initial velocity range
double K = 0.2;

double pi = 3.1415926;

double frand() // вещественное случайное число в диапазоне [0,1)
{
    return double(rand()) / RAND_MAX * 0.9999;
}

double F(double* x, int n) // сферическая функция
{
    double s = 0.0;

    // // sphere
    // for( int i = 0; i < n; i++ )
    //     s += x[i] * x[i];

    // // Rosenbrok
    // for( int i = 0; i < n - 1; i++ )
    //     s += 100 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (x[i] - 1) * (x[i] - 1);

    // Rastring
    for( int i = 0; i < n; i++ )
        s += x[i] * x[i] - 10 * cos(2 * 3.14159265 * x[i]) + 10;

    return s;
}

double create_particle(double* x, double* v, int n)
{
    double d = 0.0;
    for( int i = 0; i < n; i++ )
    {
        x[i] = a + (b - a) * frand();
        v[i] = v0 * (2 * frand() - 1);
    }

    return F(x, n);
}

void copy_data(double* x, int n, double* y)
{
    for( int i = 0; i < n; i++ )
        y[i] = x[i];
}

void update_particle(double* x, double* v, double* p, double* g, int n)
{
    double alpha = frand();
    double beta = frand(); 
    for( int i = 0; i < n; i++ )
    {
        v[i] = eta * v[i] + alpha * (p[i] - x[i]) + beta * (g[i] - x[i]);
        x[i] = x[i] + tau * v[i];
    }
    
    if( F(x,n) < F(p, n) )
        copy_data(x, n, p);
}

void print_best(double* g, double fbest, int n, int iter)
{
    cout << iter << ":   " << "Best error = " << fbest << "   " << "Solution = ";
    for( int i = 0; i < n; i++ )
        cout << g[i] << " ";
    cout << endl;
}

void print_all(ostream& f, double* x, int n, int m)
{
    for( int i = 0; i < n * m; i++ )
        f << x[i] << " ";
    f << endl;
}

double synchronizePSO(double* g, int n, int myRank, int nProc)
{
    struct
    {
        double func;
        int rank;
    } mybest, allbest;

    mybest.func = F(g, n);
    mybest.rank = myRank;

    MPI_Allreduce(&mybest, &allbest, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    MPI_Bcast(g, n, MPI_DOUBLE, allbest.rank, MPI_COMM_WORLD);

    return F(g, n);
}

void shuffle(double* x, double* v, double* p, int m, int n)
{
    for( int k = 0; k < m; k++ )
    {
        int l = rand() % m;
        for( int i = 0; i < n; i++ )
        {
            swap(x[k * n + i], x[l * n + i]);
            swap(v[k * n + i], v[l * n + i]);
            swap(p[k * n + i], p[l * n + i]);
        }
    }
}

void migrate(double* x, double* v, double* p, int m, int n, int s, int left, int right)
{
    shuffle(x, v, p, m, n);

    double *tmp = new double [s * n];

    // MPI_Sendrecv(x, s * n, MPI_DOUBLE, left, 0, x, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // MPI_Sendrecv(v, s * n, MPI_DOUBLE, left, 0, v, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // MPI_Sendrecv(p, s * n, MPI_DOUBLE, left, 0, p, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Sendrecv(x, s * n, MPI_DOUBLE, left, 0, tmp, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(x, tmp, s * n * sizeof(double));
    delete [] tmp;

    tmp = new double [s * n];
    MPI_Sendrecv(v, s * n, MPI_DOUBLE, left, 0, tmp, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(v, tmp, s * n * sizeof(double));
    delete [] tmp;

    tmp = new double [s * n];
    MPI_Sendrecv(p, s * n, MPI_DOUBLE, left, 0, tmp, s * n, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(p, tmp, s * n * sizeof(double));
    delete [] tmp;

    return;
}

void run_psoa(int n, int m, int T, int s, int kk, int myRank, int nProc)
{
    double* x = new double[n * m];
    double* v = new double[n * m];
    double* p = new double[n * m];
    double* g = new double[n];

    int left = myRank - 1;
    int right = myRank + 1;

    if (left == -1)
        left = nProc - 1;

    if (right == nProc)
        right = 0;

    // initializing
    double fbest = 1e30;
    int k0 = 0;
    for( int k = 0; k < m; k++ )
    { 
        double f = create_particle(x + n * k, v + n * k, n);
        copy_data(x + n * k, n, p + n * k);
        if( f < fbest )
        {
            fbest = f;
            k0 = k;
        }
    }
    copy_data(x + n * k0, n, g);

    int dt = T / 100;
    if( dt == 0 )
        dt = 1;

    ofstream fout("trek.dat");
    fout << (T / dt) + 1 << " " << n * m << endl;
    print_all(fout, x, n, m);

    // main loop
    // fbest = synchronizePSO(g, n, myRank, nProc);
    long double time = MPI_Wtime();
    k0 = 0;
    for( int t = 1; t <= T; t++ )
    {
        for( int k = 0; k < m; k++ )
            update_particle(x + n * k, v + n * k, p + n * k, g, n);

        for( int k = 0; k < m; k++ )
        {
            double f = F(x + n * k, n); 

            if( f < fbest )
            {
                fbest = f;
                k0 = k;
            }
        }
        copy_data(x + n * k0, n, g);

        // fbest = synchronizePSO(g, n, myRank, nProc);
        if (t % kk == 0)
            migrate(x, v, p, m, n, s, left, right);

        double fbestall = 0;
        MPI_Reduce(&fbest, &fbestall, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (myRank == 0)
            if( t % dt == 0 )
            {
                // print_best(g, fbest, n, t);
                print_best(g, fbestall, n, t); 
                // print_all(fout, x, n, m);
            }
    }
    time = MPI_Wtime() - time;
    long double timeR = 0.0;
    MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    fout.close();

    if (myRank == 0)
    {
        print_best(g, fbest, n, -1);
        cout << "> Final time of computation = " << timeR << endl;
    }

    delete[] x;
    delete[] v;
    delete[] p;
    delete[] g;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int n = atoi(argv[1]); // task dimension
    int m = atoi(argv[2]); // swarm size
    int T = atoi(argv[3]); // duration
    int s = atoi(argv[4]); // migrate size
    int k = atoi(argv[5]); // migrate iter

    srand(time(NULL) + myRank);

    run_psoa(n, m, T, s, k, myRank, nProc);

    MPI_Finalize();
    return 0;
}