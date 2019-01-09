#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <mpi.h>

using namespace std;

double frand(double a, double b)
{
	return a + (b - a) * (rand() / double(RAND_MAX));
}

int do_walk(int a, int b, int x, double p, double& t)
{
	int step = 0;
	while( x > a && x < b )
	{
		if( frand(0,1) < p )
			x += 1;
		else
			x -= 1;
		t += 1.0;
		step += 1;
	}
	return x;
}

void run_mc(int a, int b, int x, double p, int N, double &wRes, double &tRes)
{
	double t = 0.0; 
	double w = 0.0; 

	for( int i = 0; i < N; i++ )
	{
		int out = do_walk(a, b, x, p, t);
		if( out == b )
			w += 1;
	}

	wRes = w;
	tRes = t;
}

int main(int argc, char** argv)
{
	srand(time(0));
	MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	int x = atoi(argv[3]);
	double p = atof(argv[4]);
	int N = atoi(argv[5]);

	int pN = N / nProc;
	
	if (myRank == 0)
		pN += N % nProc;

	double wRes = 0.0;
	double tRes = 0.0;

	long double time = MPI_Wtime();
	run_mc(a, b, x, p, pN, wRes, tRes);
    time = MPI_Wtime() - time;


	double wResR = 0.0;
	double tResR = 0.0;
	MPI_Reduce(&wRes, &wResR, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&tRes, &tResR, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	
    long double timeR = 0.0;
	MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (myRank == 0)
	{
	    ofstream fo("output.dat", fstream::app);
		fo << wResR / N << " " << tResR / N << endl;
		fo.close();

		ofstream fs("stat.dat", fstream::app);
		fs << timeR << " " << nProc << " " << a << " " << b << " " << x << " " << p << " " << N << endl;
		fs.close();
        // cout << "> Final time of computation = " << timeR << endl;
    }


	MPI_Finalize();
	return 0;
}
