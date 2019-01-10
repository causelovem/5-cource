#include <fstream>
#include <cmath>
#include <mpi.h>

using namespace std;

int f(int* data, int i, int j, int n)
{
    int state = data[i * (n + 2) + j];
    int s = -state;
    for( int ii = i - 1; ii <= i + 1; ii++ )
        for( int jj = j - 1; jj <= j + 1; jj++ )
            s += data[ii * (n + 2) + jj];
    if( state == 0 && s == 3 )
        return 1;
    if( state == 1 && (s < 2 || s > 3) ) 
        return 0;
    return state;
}

void update_data(int n, int* data, int* temp)
{
    for( int i = 1; i <= n; i++ )
        for( int j = 1; j <= n; j++ )
            temp[i * (n + 2) + j] = f(data, i, j, n);
}

void init(int n, int* data, int* temp)
{
    for( int i = 0; i < (n + 2) * (n + 2); i++ )
        data[i] = temp[i] = 0;
    int n0 = 1 + n / 2;
    int m0 = 1 + n / 2;
    data[(n0 - 1) * (n + 2) + m0] = 1;
    data[n0 * (n + 2) + m0 + 1] = 1;
    for( int i = 0; i < 3; i++ )
        data[(n0 + 1) * (n + 2) + m0 + i - 1] = 1;
}

void setup_boundaries(int n, int* data)
{
    for( int i = 0; i < n + 2; i++ )
    {
        data[i * (n + 2) + 0] = data[i * (n + 2) + n];
        data[i * (n + 2) + n + 1] = data[i * (n + 2) + 1];
    }
    for( int j = 0; j < n + 2; j++ )
    {
        data[0 * (n + 2) + j] = data[n * (n + 2) + j];
        data[(n + 1) * (n + 2) + j] = data[1 * (n + 2) + j];
    }
}

void distribute_data(int* data, int n, int p, int myRank)
{
    MPI_Datatype batch;
    MPI_Type_vector(n, n, n * p, MPI_INT, &batch);
    MPI_Type_commit(&batch);
    MPI_Status status;

    if (myRank == 0)
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                MPI_Send(&data[i * p * n * n + j * n], 1, batch, i * p + j + 1, 0, MPI_COMM_WORLD);
    else
        MPI_Recv(&data, 1, batch, 0, 0, MPI_COMM_WORLD, &status);

    return;
}

void setup_boundaries_MPI(int* data, int n, int p, int myRank)
{
    int i = (myRank - 1) / p;
    int j = (myRank - 1) % p;

    int left = myRank - 1;
    if (j == 0)
        left += p;
    int right = myRank + 1;
    if (j == p - 1)
        right -= p;
    int above = myRank - p;
    if (i == 0)
        above += p * (p - 1);
    int below = myRank + p;
    if (i == p - 1)
        below -= p * (p - 1);
}

void run_life(int n, int T, int myRank, int dim)
{
    int N = dim * n;
    int* data;
    int* temp;

    if (myRank == 0)
    {
        data = new int[(N + 2) * (N + 2)];
        // temp = new int[(N + 2) * (N + 2)];
        init(N, data, temp);
        setup_boundaries(N, data);
    }
    else
    {
        data = new int[(n + 2) * (n + 2)];
        temp = new int[(n + 2) * (n + 2)];
    }
    distribute_data(data, n, dim, myRank);
    
    for( int t = 0; t < T; t++ )
    {
        // setup_boundaries(n, data);
        update_data(n, data, temp);
        swap(data, temp);
    }

    ofstream f("output.dat");
    for( int i = 1; i <= n; i++ )
    {
        for( int j = 1; j <= n; j++ )
            f << data[i * (n + 2) + j];
        f << endl;
    }
    f.close();

    delete[] data;
    if (myRank != 0)
        delete[] temp;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int n = atoi(argv[1]);
    int T = atoi(argv[2]);

    int dim = sqrt(nProc - 1);

    if (dim * dim != nProc - 1)
        if (myRank == 0)
            cout << ">Can not parallel program, because your proc number is not a 2-rd power of some N." << endl;


    run_life(n, T, myRank, dim);


    MPI_Finalize();
    return 0;
}