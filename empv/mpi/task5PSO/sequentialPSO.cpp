#include <iostream>
#include <fstream>
#include <cmath>

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
	return double(rand())/RAND_MAX*0.9999;
}

double F(double* x, int n) // сферическая функция
{
	double s = 0.0;
	for( int i=0; i<n; i++ )
	{
		s += x[i]*x[i];
	}

	return s;
}

double create_particle(double* x, double* v, int n)
{
	double d = 0.0;
	for( int i=0; i<n; i++ )
	{
		x[i] = a+(b-a)*frand();
		v[i] = v0*(2*frand()-1);
	}
	return F(x, n);
}

void copy_data(double* x, int n, double* y)
{
	for( int i=0; i<n; i++ )
		y[i] = x[i];
}

void update_particle(double* x, double* v, double* p, double* g, int n)
{
	double alpha = frand();
	double beta = frand(); 
	for( int i=0; i<n; i++ )
	{
		v[i] = eta*v[i] + alpha*(p[i]-x[i]) + beta*(g[i]-x[i]);
		x[i] = x[i] + tau*v[i];
	}
	if( F(x,n) < F(p, n) )
		copy_data(x, n, p);
}

void print_best(double* g, double fbest, int n)
{
	cout << fbest << " : ";
	for( int i=0; i<n; i++ )
		cout << g[i] << " ";
	cout << endl;
}

void print_all(ostream& f, double* x, int n, int m)
{
	for( int i=0; i<n*m; i++ )
		f << x[i] << " ";
	f << endl;
}

void run_psoa(int n, int m, int T)
{
	double* x = new double[n*m];
	double* v = new double[n*m];
	double* p = new double[n*m];
	double* g = new double[n];

	// initializing
	double fbest = 1e30;
	for( int k=0; k<m; k++ )
	{ 
		double f = create_particle(x+n*k,v+n*k,n);
		copy_data(x+n*k,n,p+n*k);
		if( f < fbest )
		{
			fbest = f;
			copy_data(x+n*k,n,g);
		}
	}

	int dt = T/100;
	if( dt==0 )
		dt = 1;

	ofstream fout("trek.dat");
	fout << (T/dt)+1 << " " << n*m << endl;
	print_all(fout, x, n, m);

	// main loop
	for( int t = 1; t<=T; t++ )
	{
		for( int k=0; k<m; k++ )
			update_particle(x+n*k, v+n*k, p+n*k, g, n);
		for( int k=0; k<m; k++ )
		{
			double f = F(x+n*k, n); 
			//cout << f << " ";
			if( f < fbest )
			{
				fbest = f;
				copy_data(x+n*k,n,g);
			}
		}
		if( t%dt==0 )
		{
			print_best(g, fbest, n); 
			print_all(fout, x, n, m);
		}
	}

	fout.close();

	print_best(g, fbest, n);

	delete[] x;
	delete[] v;
	delete[] p;
	delete[] g;
}

int main(int argc, char** argv)
{
	int n = atoi(argv[1]); // task dimension
	int m = atoi(argv[2]); // swarm size
	int T = atoi(argv[3]); // duration

	run_psoa(n, m, T);

	return 0;
}