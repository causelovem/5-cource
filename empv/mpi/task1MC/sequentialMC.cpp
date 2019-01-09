#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

double frand(double a, double b)
{
	return a+(b-a)*(rand()/double(RAND_MAX));
}

int do_walk(int a, int b, int x, double p, double& t)
{
	int step = 0;
	while( x>a && x<b )
	{
		if( frand(0,1)<p )
			x += 1;
		else
			x -= 1;
		t += 1.0;
		step += 1;
	}
	return x;
}

void run_mc(int a, int b, int x, double p, int N)
{
	// srand(time(0));
	double t = 0.0; 
	double w = 0.0; 

	for( int i=0; i<N; i++ )
	{
		int out = do_walk(a, b, x, p, t);
		if( out == b )
			w += 1;
	}

	ofstream f("output.dat");
	f << w/N << " " << t/N << endl;
	f.close();
}

int main(int argc, char** argv)
{
	int a = atoi(argv[1]);
	int b = atoi(argv[2]);
	int x = atoi(argv[3]);
	double p = atof(argv[4]);
	int N = atoi(argv[5]);

	run_mc(a, b, x, p, N);

	return 0;
}
