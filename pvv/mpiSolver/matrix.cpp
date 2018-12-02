#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <omp.h>
#include <mpi.h>

using namespace std;

double globtime = 0.0;

class Vector;

class Matrix
{
    private:
        double *A = NULL;
        int *JA = NULL;
        int *IA = NULL;
        // int *rows = NULL;
        int sizeIA = 0;
        int sizeA = 0;

    public:
        // (startI, startJ, startK, Nxp, Nyp, Nzp)
        // Matrix(int Nx, int Ny, int Nz, int Px, int Py, int Pz)
        Matrix(int Px, int Py, int Pz, int Nxp, int Nyp, int Nzp, int Nx, int Ny, int Nz, int *rows)
        {
            sizeIA = Nxp * Nyp * Nzp + 1;

            int Ppx = Nxp + Px,
                Ppy = Nyp + Py,
                Ppz = Nzp + Pz;

            sizeA = 10 * sizeIA;
            // cout << Nx << ' ' << Ny << ' ' << Nz << ' ' << Ppx << ' ' << Ppy << ' ' << Ppz << ' ' << Px << ' ' << Py << ' ' << Pz <<  endl;


            IA = new int [sizeIA];
            JA = new int [sizeA];
            A = new double [sizeA];
            // rows = new int [sizeIA];
            IA[0] = 0;

            int ia = 0, ija = 0, iia = 1;
            int cnt = 0;

            for (int k = Pz; k < Ppz; k++)
                for (int j = Py; j < Ppy; j++)
                    for (int i = Px; i < Ppx; i++)
                    {
                        int num = i + Nx * j + Nx * Ny * k;
                        rows[iia - 1] = num;
                        double sum = 0.0;
                        // double s = sin(i + j + 1);

                        if (k > 0)
                        {
                            // A[ia] = s;
                            A[ia] = sin(2 * num - Nx * Ny + 1);
                            sum += fabs(A[ia++]);
                            JA[ija++] = num - Nx * Ny;
                            cnt++;
                        }

                        if (j > 0)
                        {
                            // A[ia] = s;
                            A[ia] = sin(2 * num - Nx + 1);
                            sum += fabs(A[ia++]);
                            JA[ija++] = num - Nx;
                            cnt++;
                        }

                        if (i > 0)
                        {
                            // A[ia] = s;
                            A[ia] = sin(2 * num - 1 + 1);
                            sum += fabs(A[ia++]);
                            JA[ija++] = num - 1;
                            cnt++;
                        }

                        int tia = ia++, tija = ija++;

                        if (i < Nx - 1)
                        {
                            // A[ia] = s;
                            A[ia] = sin(2 * num + 1 + 1);
                            sum += fabs(A[ia++]);
                            JA[ija++] = num + 1;
                            cnt++;
                        }

                        if (j < Ny - 1)
                        {
                            // A[ia] = s;
                            A[ia] = sin(2 * num + Nx + 1);
                            sum += fabs(A[ia++]);
                            JA[ija++] = num + Nx;
                            cnt++;
                        }

                        if (k < Nz - 1)
                        {
                            // A[ia] = s;
                            A[ia] = sin(2 * num + Nx * Ny + 1);
                            sum += fabs(A[ia++]);
                            JA[ija++] = num + Nx * Ny;
                            cnt++;
                        }

                        A[tia] = fabs(sum) * 1.1;
                        JA[tija] = num;
                        IA[iia++] = ++cnt;
                    }

            sizeA = ija;
        }

        Matrix(const Matrix &mat)
        {
            sizeIA = mat.sizeIA;
            sizeA = mat.sizeA;

            IA = new int [sizeIA];
            JA = new int [sizeA];
            A = new double [sizeA];  

            // memcpy(JA, mat.JA, sizeA * sizeof(int));
            // memcpy(A, mat.A, sizeA * sizeof(double));
            // memcpy(IA, mat.IA, sizeIA * sizeof(int));
            
            #pragma omp parallel for
            for (int i = 0; i < sizeA; i++)
            {
                JA[i] = mat.JA[i];
                A[i] = mat.A[i];
            }   

            #pragma omp parallel for
            for (int i = 0; i < sizeIA; i++)
                IA[i] = mat.IA[i];
        }

        ~Matrix()
        {
            delete [] A;
            delete [] JA;
            delete [] IA;
            // delete [] rows;

            sizeIA = 0;
            sizeA = 0;
        }

        Matrix(const Matrix &mat, int k)
        {
            for (int i = 0; i < mat.sizeIA - 1; i++)
                if (mat.get(i, i) > 0.000001)
                    sizeA++;

            sizeIA = mat.sizeIA;

            IA = new int [sizeIA];
            JA = new int [sizeA];
            A = new double [sizeA];
            IA[0] = 0;
            int cnt = 0;

            int ia = 0, iia = 1;

            for (int i = 0; i < sizeIA - 1; i++)
            {
                double diag = mat.get(i, i);

                if (diag > 0.000001)
                {
                    A[ia] = 1.0 / diag;
                    // A[ia] = diag;
                    JA[ia++] = i;
                    cnt++;
                }

                IA[iia++] = cnt;
            }
        }

        double get(int i, int j) const
        {
            for (int k = IA[i]; k < IA[i + 1]; k++)
                if (j == JA[k])
                    return A[k];
            return 0.0;
        }

        void printCSR()
        {
            cout << "A" << endl;
            for (int i = 0; i < sizeA; i++)
                cout << A[i] << ' ';
            cout << endl;

            cout << "JA" << endl;
            for (int i = 0; i < sizeA; i++)
                cout << JA[i] << ' ';
            cout << endl;

            cout << "IA" << endl;
            for (int i = 0; i < sizeIA; i++)
                cout << IA[i] << ' ';
            cout << endl;

            // cout << "rows" << endl;
            // for (int i = 0; i < sizeIA; i++)
            //     cout << rows[i] << ' ';
            // cout << endl;

            return;
        }

        // doesnt work in parallel
        void print(int flg = 0)
        {
            if (flg == 1)
                printCSR();

            cout << "Matrix" << endl;
            for (int i = 0; i < sizeIA - 1; i++)
            {
                for (int j = 0; j < sizeIA - 1; j++)
                {
                    // cout.width(6);
                    // cout.setf(ios::left);
                    cout << get(i, j) << ' ';
                }
                cout << endl;
            }

            return;
        }

        friend int SpMV(const Matrix &mat, const Vector &vec, Vector &res);
};

class Vector
{
    private:
        double *A = NULL;
        int size = 0;
    public:
        // Vector(int s)
        // {
        //     size = s;
        //     A = new double [s];

        //     #pragma omp parallel for
        //     for (int i = 0; i < size; i++)
        //         A[i] = sin(i);
        //         // A[i] = rand() % 10;
        // }

        Vector(int s, int start)
        {
            size = s;
            A = new double [size];

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                A[i] = sin(start + i);
                // A[i] = rand() % 10;
        }

        Vector(const Vector &vec)
        {
            size = vec.size;
            A = new double [size];

            // memcpy(A, vec.A, size * sizeof(double));

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                A[i] = vec.A[i];
        }

        Vector(int s, double c)
        {
            size = s;
            A = new double [size];

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                A[i] = c;
        }

        Vector(int s, double *vec)
        {
            size = s;
            A = new double [size];

            // memcpy(A, vec, size * sizeof(double));

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                A[i] = vec[i];
        }

        ~Vector()
        {
            delete [] A;
            size = 0;
        }

        Vector operator = (const Vector &vec)
        {
            if (A != NULL)
                delete [] A;

            size = vec.size;
            A = new double [size];

            const int ss = (size / 4) * 4;
            // float time = omp_get_wtime();

            // memcpy(A, vec.A, size * sizeof(double));

            #pragma omp parallel for
            for (int i = 0; i < ss; i += 4)
            {
                A[i] = vec.A[i];
                A[i + 1] = vec.A[i + 1];
                A[i + 2] = vec.A[i + 2];
                A[i + 3] = vec.A[i + 3];
            }

            // #pragma omp parallel for
            for (int i = ss; i < size; i++)
                A[i] = vec.A[i];

            // #pragma omp parallel for
            // for (int i = 0; i < size; i++)
            //     A[i] = vec.A[i];

            // time = omp_get_wtime() - time;
            // globtime += time;
            // cout << "> Time of computation = " << (time) << endl;

            return *this;
        }

        Vector operator = (const double c)
        {
            const int ss = (size / 4) * 4;

            // float time = omp_get_wtime();
            #pragma omp parallel for
            for (int i = 0; i < ss; i += 4)
            {
                A[i] = c;
                A[i + 1] = c;
                A[i + 2] = c;
                A[i + 3] = c;
            }

            // #pragma omp parallel for
            for (int i = ss; i < size; i++)
                A[i] = c;

            // #pragma omp parallel for
            // for (int i = 0; i < size; i++)
            //     A[i] = c;

            // time = omp_get_wtime() - time;
            // globtime += time;
            // cout << "> Time of computation = " << (time) << endl;

            return *this;
        }

        void print()
        {
            cout << "Vector" << endl;
            for (int i = 0; i < size; i++)
            {
                if (fabs(A[i]) < 0.000001)
                    cout << 0 << ' ';
                else
                    cout << A[i] << ' ';
            }
            cout << endl;

            return;
        }
    
        friend double dot(const Vector &vec1, const Vector &vec2)
        {
            if (vec1.size != vec2.size)
                cout << "Can't dot: different lenghts!" << endl;

            double res = 0.0;

            const int ss = (vec1.size / 4) * 4;
            // float time = omp_get_wtime();

            #pragma omp parallel for reduction(+:res)
            for (int i = 0; i < ss; i += 4)
            {
                res += vec1.A[i] * vec2.A[i];
                res += vec1.A[i + 1] * vec2.A[i + 1];
                res += vec1.A[i + 2] * vec2.A[i + 2];
                res += vec1.A[i + 3] * vec2.A[i + 3];
            }

            // #pragma omp parallel for reduction(+:res)
            for (int i = ss; i < vec1.size; i++)
                res += vec1.A[i] * vec2.A[i];

            // #pragma omp parallel for reduction(+:res)
            // for (int i = 0; i < vec1.size; i++)
            //     res += vec1.A[i] * vec2.A[i];

            // time = omp_get_wtime() - time;
            // globtime += time;
            // cout << "> Time of computation = " << (time) << endl;

            return res;
        }

        friend int axpby(Vector &vec1, const Vector &vec2, double a, double b)
        {
            if (vec1.size != vec2.size)
            {
                cout << "Can't axpby: different lenghts!" << endl;
                return -1;
            }

            const int ss = (vec1.size / 4) * 4;
            // float time = omp_get_wtime();

            #pragma omp parallel for
            for (int i = 0; i < ss; i += 4)
            {
                vec1.A[i] = a * vec1.A[i] + b * vec2.A[i];
                vec1.A[i + 1] = a * vec1.A[i + 1] + b * vec2.A[i + 1];
                vec1.A[i + 2] = a * vec1.A[i + 2] + b * vec2.A[i + 2];
                vec1.A[i + 3] = a * vec1.A[i + 3] + b * vec2.A[i + 3];
            }

            // #pragma omp parallel for
            for (int i = ss; i < vec1.size; i++)
                vec1.A[i] = a * vec1.A[i] + b * vec2.A[i];

            // #pragma omp parallel for
            // for (int i = 0; i < vec1.size; i++)
            //     vec1.A[i] = a * vec1.A[i] + b * vec2.A[i];

            // time = omp_get_wtime() - time;
            // globtime += time;
            // cout << "> Time of computation = " << (time) << endl;

            return 0;
        }

        friend int SpMV(const Matrix &mat, const Vector &vec, Vector &res)
        {
            if (mat.sizeIA - 1 != vec.size)
            {
                cout << "Can't SpMV: different lenghts!" << endl;
                return -1;
            }

            // for (int i = 0; i < vec.size; i++)
            //     for (int j = 0; j < vec.size; j++)
            //         res.A[i] += mat.get(i, j) * vec.A[j];

            res = 0.0;

            // much faster
            #pragma omp parallel for
            for (int i = 0; i < mat.sizeIA - 1; i++)
            {
                res.A[i] = 0.0;
                for (int j = mat.IA[i]; j < mat.IA[i + 1]; j++)
                    res.A[i] += mat.A[j] * vec.A[mat.JA[j]];
            }

            return 0;
        }
};


int solve(int N, Matrix &A, Vector &BB, double tol, int maxit, int debug)
{
    int I = 0;
    Vector XX(N, 0.0);
    Matrix DD(A, 1);

    Vector PP(N, 0.0), PP2(N, 0.0), RR(N, 0.0), RR2(N, 0.0), TT(N, 0.0), VV(N, 0.0), SS(N, 0.0), SS2(N, 0.0);
    double initres = 0.0, res = 0.0, mineps = 1e-15, eps = 0.0;
    double Rhoi_1 = 1.0, alphai = 1.0, wi = 1.0, betai_1 = 1.0, Rhoi_2 = 1.0, alphai_1 = 1.0, wi_1 = 1.0, RhoMin = 1e-60;

    RR = BB;
    RR2 = BB;
    initres = sqrt(dot(RR, RR));
    eps = max(mineps, tol * initres);
    res = initres;


    for (I = 0; I < maxit; I++)
    {
        if (debug != 0)
            cout << "It = " << I << "; res = " << res << "; tol = " << res / initres << endl;

        if (res < eps)
            break;
        if (res > initres / mineps)
            return -1;

        if (I == 0)
            Rhoi_1 = initres * initres;
        else
            Rhoi_1 = dot(RR2, RR);
        if (fabs(Rhoi_1) < RhoMin)
            return -1;

        if (I == 0)
            PP = RR;
        else
        {
            betai_1 = (Rhoi_1 * alphai_1) / (Rhoi_2  * wi_1);
            axpby(PP, RR, betai_1, 1.0);
            axpby(PP, VV, 1.0, -wi_1 * betai_1);
        }

        SpMV(DD, PP, PP2);
        SpMV(A, PP2, VV);

        alphai = dot(RR2, VV);

        if (fabs(alphai) < RhoMin)
            return -3;

        alphai = Rhoi_1 / alphai;

        SS = RR;
        axpby(SS, VV, 1.0, -alphai);

        SpMV(DD, SS, SS2);
        SpMV(A, SS2, TT);

        wi = dot(TT, TT);
        if (fabs(wi) < RhoMin)
            return -4;
        wi = dot(TT, SS) / wi;
        if (fabs(wi) < RhoMin)
            return -5;

        axpby(XX, PP2, 1.0, alphai);
        axpby(XX, SS2, 1.0, wi);

        RR = SS;

        axpby(RR, TT, 1.0, -wi);

        alphai_1 = alphai;

        Rhoi_2 = Rhoi_1;
        wi_1 = wi;

        res = sqrt(dot(RR, RR));
    }
    cout << "> Final discrepancy = " << res << endl;

    // if (debug != 0)
    // {
        // Vector check(N, 0.0);

        // SpMV(A, XX, check);

        // cout << "> Solution" << endl;
        // XX.print();
        // cout << "> Ax" << endl;
        // check.print();
        // cout << "> BB" << endl;
        // BB.print();
    // }

    return I;
}

int testFunc(int Nx, int Ny, int Nz, int N)
{
    int thread[6] = {1, 2, 4, 8, 10, 16};
    int xx[5] = {10, 10, 10, 100, 100};
    int yy[5] = {10, 10, 100, 100, 100};
    int zz[5] = {10, 100, 100, 100, 1000};
    for (int k = 0; k < 5; k++)
    {
        Nx = xx[k];
        Ny = yy[k];
        Nz = zz[k];
        N = Nx * Ny * Nz;
        cout << Nx << ' ' << Ny << ' ' << Nz << ' ' << N << endl;

        long double seqTimeDot = 0.0;
        long double seqTimeAxpby = 0.0;
        long double seqTimeSpmv = 0.0;

        for (int i = 0; i < 6; i++)
        {
            Vector testVec1(N, 0), testVec2(N, 0);
            Vector testRes(N, 0.0);
            Matrix testMat(Nx, Ny, Nz, 1, 1, 1, 1, 1, 1, NULL);

            int corner = 4 * 8;
            int edge = 5 * 4 * ((Nx - 2) + (Ny - 2) + (Nz - 2));
            int face = 6 * 2 * (((Nx - 2) * (Ny - 2)) + ((Nx - 2) * (Nz - 2)) + ((Ny - 2) * (Nz - 2)));
            int inner = 7 * (Nx - 2) * (Ny - 2) * (Nz - 2);
            int sizeA = corner + edge + face + inner;
            // int num = 1 << i;
            int num = thread[i];

            omp_set_num_threads(num);
            cout << "> Threads = " << num << endl << endl;

            long double operations = 1e-9 * 2 * N;
            long double testTime = omp_get_wtime();
            cout << "> DOT = " << dot(testVec1, testVec2) << endl;
            testTime = omp_get_wtime() - testTime;
            cout << "> Time of DOT = " << testTime << endl;
            cout << "> GFLOPS = " << operations / testTime << endl;
            if (i == 0)
                seqTimeDot = testTime;
            cout << "> SpeedUp = " << seqTimeDot / testTime << endl;
            cout << endl;


            operations = 1e-9 * 3 * N;
            testTime = omp_get_wtime();
            axpby(testVec1, testVec2, 1.0, 2.0);
            testTime = omp_get_wtime() - testTime;
            cout << "> AXPBY L2 norm = " << sqrt(dot(testVec1, testVec1)) << endl;
            cout << "> Time of AXPBY = " << testTime << endl;
            cout << "> GFLOPS = " << operations / testTime << endl;
            if (i == 0)
                seqTimeAxpby = testTime;
            cout << "> SpeedUp = " << seqTimeAxpby / testTime << endl;
            cout << endl;


            operations = 1e-9 * 2 * N * 7;
            testTime = omp_get_wtime();
            SpMV(testMat, testVec2, testVec1);
            testTime = omp_get_wtime() - testTime;
            cout << "> SpMV L2 norm = " << sqrt(dot(testVec1, testVec1)) << endl;
            cout << "> Time of SpMV = " << testTime << endl;
            cout << "> GFLOPS = " << operations / testTime << endl;
            if (i == 0)
                seqTimeSpmv = testTime;
            cout << "> SpeedUp = " << seqTimeSpmv / testTime << endl;
            cout << endl;


            if ((k == 4) && (i == 5))
            {
                Matrix A(Nx, Ny, Nz, 1, 1, 1, 1, 1, 1, NULL);
                Vector BB(N, 0);

                cout << "> Solver test..." << endl << endl;
                long double seqTimeSolve = 0.0;
                for (int j = 0; j < 6; j++)
                {
                    omp_set_num_threads(thread[j]);
                    cout << "> Number of threads = " << thread[j] << endl;
                    // 5 dot 6 axpby 4 spmv N diag
                    long double operations = 1e-9 * (5 * (2 * N) + 6 * (3 * N) + 4 * (2 * N * 7) + N);
                    long double time = omp_get_wtime();
                    int res = solve(N, A, BB, 0.000000001, 1000, 1);
                    time = omp_get_wtime() - time;
                    cout << "> Numder of iters = " << res << endl;
                    cout << "> Final time of computation = " << time << endl;
                    cout << "> GFLOPS = " << operations / time << endl;
                    if (j == 0)
                        seqTimeSolve = time;
                    cout << "> SpeedUp = " << seqTimeSolve / time << endl;
                    cout << endl;
                }
            }
        }
    }

    return 0;
}

int main (int argc, char **argv)
{
    srand(time(0));

    if (argc != 11)
    {
        cout << "> Wrong number of in params, check your command" << endl;
        cout << "<Nx> <Ny> <Nz> <tol> <maxit> <omp> <Px> <Py> <Pz> <debug>" << endl;
        return -1;
    }

    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    omp_set_num_threads(atoi(argv[6]));

    int Nx = atoi(argv[1]), Ny = atoi(argv[2]), Nz = atoi(argv[3]);
    int N = Nx * Ny * Nz, maxit = atoi(argv[5]);
    double tol = atof(argv[4]);
    int Px = atoi(argv[7]), Py = atoi(argv[8]), Pz = atoi(argv[9]);
    int debug = atoi(argv[10]);

    if (Px * Py * Pz != nProc)
    {
        cout << "> Px * Py * Pz != nProc, check your command" << endl;
        return -1;
    }
    // cout << Nx << ' ' << Ny << ' ' << Nz << ' ' << tol << ' ' << maxit << endl;

    int kk = myRank / (Px * Py),
        jj = (myRank % (Px * Py)) / Pz,
        ii = myRank - kk * (Px * Py) - jj * Pz;

    // int dim = Nx * Ny * Nz;

    int Nxp = Nx / Px,
        Nyp = Ny / Py,
        Nzp = Nz / Pz;

    int startI = ii * Nxp,
        startJ = jj * Nyp,
        startK = kk * Nzp;

    int Np = N / nProc;

    if (myRank == nProc - 1)
    {
        Nxp += Nx % Px;
        Nyp += Ny % Py;
        Nzp += Nz % Pz;

        Np += N % nProc;
    }

    int rowsSize = Nxp * Nyp * Nzp;


    rowsSize = 10 * rowsSize;

    int *rows = new int [rowsSize];
    int *global2loc = new int [N];


    for (int i = 0; i < rowsSize; i++)
        rows[i] = -1;

    Matrix A(startI, startJ, startK, Nxp, Nyp, Nzp, Nx, Ny, Nz, rows);
    for (int i = 0; i < N; i++)
        global2loc[i] = -1;

    for (int i = 0; i < rowsSize; i++)
        if (rows[i] != -1)
            global2loc[rows[i]] = i;

    // if (myRank == 2)
    // {
    //     for (int i = 0; i < rowsSize; i++)
    //         cout << rows[i] << ' ';
    //     cout << endl << endl;
    //     for (int i = 0; i < N; i++)
    //         cout << global2loc[i] << ' ';
    //     cout << endl;
    // }




    delete [] rows;
    delete [] global2loc;
    MPI_Finalize();
    return 0;


    Vector BB(Nxp * Nyp * Nzp, Np);

    // if (myRank == 0)
    // {
    //     // cout << Nxp << ' ' << Nyp << ' ' << Nzp << ' ' << startI << ' ' << startJ << ' ' << startK << endl;
    //     A.print(1);
    //     // cout << A.get(0, 100) << endl;
    // }




    if (debug != 0)
        testFunc(Nx, Ny, Nz, N);
    else
    {
        omp_set_num_threads(atoi(argv[6]));
        // 5 dot 6 axpby 4 spmv N diag
        long double operations = 1e-9 * (5 * (2 * N) + 6 * (3 * N) + 4 * (2 * N * 7) + N);
        cout << "> Number of threads = " << atoi(argv[6]) << endl;
        long double time = omp_get_wtime();
        int res = solve(N, A, BB, tol, maxit, debug);
        time = omp_get_wtime() - time;
        cout << "> Number of iters = " << res << endl;
        cout << "> Final time of computation = " << time << endl;
        cout << "> GFLOPS = " << operations / time << endl;
    }
    // cout << "> Operation time = " << globtime << endl;

    // delete [] rows;
    // delete [] global2loc;

    MPI_Finalize();
    return 0;
}
