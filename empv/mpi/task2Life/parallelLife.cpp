#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
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

// void init(int n, int* data, int* temp)
void init(int n, int* data)
{
    for( int i = 0; i < (n + 2) * (n + 2); i++ )
        // data[i] = temp[i] = 0;
        data[i] = 0;
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
    MPI_Type_vector(n + 2, n + 2, n * p + 2, MPI_INT, &batch);
    MPI_Type_commit(&batch);

    if (myRank == 0)
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                MPI_Send(&(data[i * (n * p + 2) * n + j * n]), 1, batch, i * p + j + 1, 0, MPI_COMM_WORLD);
    else
        MPI_Recv(data, (n + 2) * (n + 2), MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


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
        above += p * p;
    int below = myRank + p;
    if (i == p - 1)
        below -= p * p;

    MPI_Datatype batchHor;
    MPI_Type_vector(1, n + 2, n + 2, MPI_INT, &batchHor);
    MPI_Type_commit(&batchHor);

    MPI_Datatype batchVert;
    MPI_Type_vector(n + 2, 1, n + 2, MPI_INT, &batchVert);
    MPI_Type_commit(&batchVert);


    MPI_Sendrecv(&(data[1]), 1, batchVert, left, 0, &(data[n + 1]), 1, batchVert, right, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&(data[n]), 1, batchVert, right, 0, &(data[0]), 1, batchVert, left, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&(data[n + 2]), 1, batchHor, above, 0, &(data[(n + 1) * (n + 2)]), 1, batchHor, below, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&(data[(n + 2) * n]), 1, batchHor, below, 0, &(data[0]), 1, batchHor, above, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Type_free(&batchHor);
    MPI_Type_free(&batchVert);
    return;
}

void collect_data(int* data, int n, int p, int myRank)
{
    MPI_Datatype batch;
    MPI_Type_vector(n + 2, n + 2, n * p + 2, MPI_INT, &batch);
    MPI_Type_commit(&batch);

    if (myRank == 0)
        for (int i = 0; i < p; i++)
            for (int j = 0; j < p; j++)
                MPI_Recv(&(data[i * (n * p + 2) * n + j * n]), 1, batch, i * p + j + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else
        MPI_Send(data, (n + 2) * (n + 2), MPI_INT, 0, 0, MPI_COMM_WORLD);

    MPI_Type_free(&batch);
    return;
}

void run_life(int n, int T, int myRank, int dim, int nProc)
{
    int N = dim * n;
    int* data;
    int* temp;


    if (myRank == 0)
    {
        data = new int[(N + 2) * (N + 2)];
        init(N, data);
        setup_boundaries(N, data);
    }
    else
    {
        data = new int[(n + 2) * (n + 2)];
        temp = new int[(n + 2) * (n + 2)];
    }
    distribute_data(data, n, dim, myRank);
    
    long double time = MPI_Wtime();
    for( int t = 0; t < T; t++ )
    {
        if (myRank != 0)
        {
            update_data(n, data, temp);
            swap(data, temp);
            setup_boundaries_MPI(data, n, dim, myRank);
        }
    }
    time = MPI_Wtime() - time;
    long double timeR = 0.0;
    MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


    collect_data(data, n, dim, myRank);


    if (myRank == 0)
    {
        ofstream fo("output.dat");
        for( int i = 1; i <= N; i++ )
        {
            for( int j = 1; j <= N; j++ )
                fo << data[i * (N + 2) + j];
            fo << endl;
        }
        fo.close();

        ofstream fs("stat.dat", fstream::app);
        fs << timeR << " " << nProc << " " << n << " " << T << endl;
        fs.close();
        cout << "> Final time of computation = " << timeR << endl;
    }

    delete[] data;
    if (myRank != 0)
        delete[] temp;

    return;
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

    if ((dim * dim != nProc - 1) || (nProc == 1))
    {
        if (myRank == 0)
            cout << ">Can not parallel program, because your proc number is not like p ^ 2 + 1." << endl;
        MPI_Finalize();
        return 0;
    }

    // long double time = MPI_Wtime();

    run_life(n, T, myRank, dim, nProc);
    
    // time = MPI_Wtime() - time;
    // if (myRank == 0)
    //     time = 0.0;
    // long double timeR = 0.0;
    // MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // if (myRank == 0)
    // {
    //     ofstream fs("stat.dat", fstream::app);
    //     fs << timeR << " " << nProc << " " << n << " " << T << endl;
    //     fs.close();
    //     cout << "> Final time of computation = " << timeR << endl;
    // }


    MPI_Finalize();
    return 0;
}