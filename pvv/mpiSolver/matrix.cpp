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
        int sizeIA = 0;
        int sizeA = 0;

    public:
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

        friend void makeHalo(const Matrix &mat, int *global2loc, int *part, int *rows, int ** &halo, int ** &sHalo, int haloSize)
        {
            for (int i = 0; i < mat.sizeIA - 1; i++)
            {
                for (int k = mat.IA[i]; k < mat.IA[i + 1]; k++)
                    if (global2loc[mat.JA[k]] ==  -1)
                    {
                        for (int j = 0; j < haloSize; j++)
                            if (halo[part[mat.JA[k]]][j] == -1)
                            {
                                halo[part[mat.JA[k]]][j] = mat.JA[k];
                                break;
                            }

                        for (int j = 0; j < haloSize; j++)
                            if (sHalo[part[mat.JA[k]]][j] == -1)
                            {
                                sHalo[part[mat.JA[k]]][j] = rows[i];
                                break;
                            }
                    }
            }

            return;
        }

        friend void changeJA(Matrix &mat, int *global2loc)
        {
            for (int i = 0; i < mat.sizeA; i++)
                mat.JA[i] = global2loc[mat.JA[i]];
            return;
        }
};

class Vector
{
    private:
        double *A = NULL;
        int size = 0;
        int locSize = 0;
    public:
        Vector(int s)
        {
            size = s;
            locSize = s;
            A = new double [s];

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                A[i] = sin(i);
                // A[i] = rand() % 10;
        }

        Vector(int selfSize, int haloSize, int *rows)
        {
            size = selfSize;
            locSize = selfSize + haloSize;
            A = new double [locSize];

            #pragma omp parallel for
            for (int i = 0; i < selfSize; i++)
                A[i] = sin(rows[i]);
        }

        Vector(const Vector &vec)
        {
            size = vec.size;
            locSize = vec.locSize;
            A = new double [locSize];

            // memcpy(A, vec.A, size * sizeof(double));

            #pragma omp parallel for
            for (int i = 0; i < locSize; i++)
                A[i] = vec.A[i];
        }

        Vector(int selfSize, int haloSize, double c)
        {
            size = selfSize;
            locSize = selfSize + haloSize;
            A = new double [locSize];

            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                A[i] = c;
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
            locSize = vec.locSize;
            A = new double [locSize];

            const int ss = (locSize / 4) * 4;
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

            for (int i = ss; i < locSize; i++)
                A[i] = vec.A[i];

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

            for (int i = ss; i < vec1.size; i++)
                res += vec1.A[i] * vec2.A[i];

            MPI_Allreduce(&res, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

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

            for (int i = ss; i < vec1.size; i++)
                vec1.A[i] = a * vec1.A[i] + b * vec2.A[i];

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

            res = 0.0;

            #pragma omp parallel for
            for (int i = 0; i < mat.sizeIA - 1; i++)
            {
                res.A[i] = 0.0;
                for (int j = mat.IA[i]; j < mat.IA[i + 1]; j++)
                    res.A[i] += mat.A[j] * vec.A[mat.JA[j]];
            }

            return 0;
        }

        friend void sync(Vector &vec, int ** &halo, int ** &sHalo, int haloOptSize)
        {
            MPI_Request *ss = new MPI_Request [haloOptSize];
            MPI_Request *rr = new MPI_Request [haloOptSize];

            double **recv = new double* [haloOptSize]; // halo
            double **send = new double* [haloOptSize]; // sHalo

            for (int i = 0; i < haloOptSize; i++)
            {
                recv[i] = new double [halo[i][1]];
                send[i] = new double [sHalo[i][1]];

                for (int j = 0; j < sHalo[i][1]; j++)
                    send[i][j] = vec.A[sHalo[i][j + 2]];

                MPI_Isend(send[i], sHalo[i][1], MPI_DOUBLE, sHalo[i][0], 0, MPI_COMM_WORLD, &(ss[i]));
                MPI_Irecv(recv[i], halo[i][1], MPI_DOUBLE, halo[i][0], MPI_ANY_TAG, MPI_COMM_WORLD, &(rr[i]));
            }

            MPI_Waitall(haloOptSize, ss, MPI_STATUS_IGNORE);
            MPI_Waitall(haloOptSize, rr, MPI_STATUS_IGNORE);

            for (int i = 0; i < haloOptSize; i++)
                for (int j = 0; j < halo[i][1]; j++)
                    vec.A[halo[i][j + 2]] = recv[i][j];


            for (int i = 0; i < haloOptSize; i++)
            {
                delete [] recv[i];
                delete [] send[i];
            }
            delete [] recv;
            delete [] send;

            delete [] ss;
            delete [] rr;
            return;
        }
};


int solve(int N, Matrix &A, Vector &BB, double tol, int maxit, int debug, int ** &halo, int ** &sHalo, int haloOptSize,
        int rowsSize, int haloRedSize, int *rows)
{
    int I = 0;
    Vector XX(rowsSize, haloRedSize, rows);
    Matrix DD(A, 1);

    Vector PP(rowsSize, haloRedSize, rows), PP2(rowsSize, haloRedSize, rows),
    RR(rowsSize, haloRedSize, rows), RR2(rowsSize, haloRedSize, rows), TT(rowsSize, haloRedSize, rows),
    VV(rowsSize, haloRedSize, rows), SS(rowsSize, haloRedSize, rows), SS2(rowsSize, haloRedSize, rows);
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

        sync(PP, halo, sHalo, haloOptSize);
        SpMV(DD, PP, PP2);
        sync(PP2, halo, sHalo, haloOptSize);
        SpMV(A, PP2, VV);

        alphai = dot(RR2, VV);

        if (fabs(alphai) < RhoMin)
            return -3;

        alphai = Rhoi_1 / alphai;

        SS = RR;
        axpby(SS, VV, 1.0, -alphai);

        sync(SS, halo, sHalo, haloOptSize);
        SpMV(DD, SS, SS2);
        sync(SS2, halo, sHalo, haloOptSize);
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

int testFunc(int Px, int Py, int Pz, int myRank, int nProc)
{
    int thread[6] = {1, 2, 4, 8, 10, 16};
    int xx[5] = {10, 10, 10, 100, 100};
    int yy[5] = {10, 10, 100, 100, 100};
    int zz[5] = {10, 100, 100, 100, 1000};
    for (int k = 0; k < 5; k++)
    {
        int Nx = xx[k];
        int Ny = yy[k];
        int Nz = zz[k];
        int N = Nx * Ny * Nz;

        int kk = myRank / (Px * Py),
        jj = (myRank % (Px * Py)) / Px,
        ii = myRank - kk * (Px * Py) - jj * Px;

        int Nxp = int(ceil(double(Nx) / double(Px))),
            Nyp = int(ceil(double(Ny) / double(Py))),
            Nzp = int(ceil(double(Nz) / double(Pz)));

        int startI = ii * Nxp,
            startJ = jj * Nyp,
            startK = kk * Nzp;

        if (ii == Px - 1)
            Nxp = Nx - ii * Nxp;
        if (jj == Py - 1)
            Nyp = Ny - jj * Nyp;
        if (kk == Pz - 1)
            Nzp = Nz - kk * Nzp;

        int rowsSize = Nxp * Nyp * Nzp;

        int *rows = new int [rowsSize];
        int *global2loc = new int [N];
        int *part = new int [N];

        int **halo = NULL, **sHalo = NULL;

        for (int i = 0; i < rowsSize; i++)
            rows[i] = -1;

        // part
        for (int i = 0; i < nProc; i++)
        {
            int kk = i / (Px * Py),
                jj = (i % (Px * Py)) / Px,
                ii = i - kk * (Px * Py) - jj * Px;

            int Nxp = int(ceil(double(Nx) / double(Px))),
                Nyp = int(ceil(double(Ny) / double(Py))),
                Nzp = int(ceil(double(Nz) / double(Pz)));

            int startI = ii * Nxp,
                startJ = jj * Nyp,
                startK = kk * Nzp;

            if (ii == Px - 1)
                Nxp = Nx - ii * Nxp;
            if (jj == Py - 1)
                Nyp = Ny - jj * Nyp;
            if (kk == Pz - 1)
                Nzp = Nz - kk * Nzp;

            int Ppx = Nxp + startI,
                Ppy = Nyp + startJ,
                Ppz = Nzp + startK;

            for (int k = startK; k < Ppz; k++)
                for (int j = startJ; j < Ppy; j++)
                    for (int l = startI; l < Ppx; l++)
                        part[l + Nx * j + Nx * Ny * k] = i;
        }


        Matrix A(startI, startJ, startK, Nxp, Nyp, Nzp, Nx, Ny, Nz, rows);

        // g2l
        for (int i = 0; i < N; i++)
            global2loc[i] = -1;
        for (int i = 0; i < rowsSize; i++)
            if (rows[i] != -1)
                global2loc[rows[i]] = i;

        halo = new int* [nProc];
        sHalo = new int* [nProc];

        // halo shalo
        int haloSize = 2 * (Nxp * Nyp + Nxp * Nzp + Nyp * Nzp);
        for (int i = 0; i < nProc; i++)
        {
            halo[i] = new int [haloSize];
            sHalo[i] = new int [haloSize];

            for (int j = 0; j < haloSize; j++)
            {
                halo[i][j] = -1;
                sHalo[i][j] = -1;
            }
        }
        makeHalo(A, global2loc, part, rows, halo, sHalo, haloSize);

        int haloRedSize = 6 * haloSize;
        int *haloRed = new int [haloRedSize];

        int curHR = 0;
        for (int i = 0; i < nProc; i++)
            if (halo[i][0] != -1)
                for (int j = 0; j < haloSize; j++)
                    if (halo[i][j] != -1)
                        haloRed[curHR++] = halo[i][j];
        haloRedSize = curHR;

        // update g2l
        for (int i = 0; i < haloRedSize; i++)
            global2loc[haloRed[i]] = i + rowsSize;
        changeJA(A, global2loc);

        // optimized halo
        int *tmpH = new int [nProc];
        int *tmpSH = new int [nProc];
        int haloOptSize = 0;

        for (int i = 0; i < nProc; i++)
        {
            tmpH[i] = 0;
            tmpSH[i] = 0;
        }

        for (int i = 0; i < nProc; i++)
        {
            if (halo[i][0] != -1)
            {
                haloOptSize++;
                for (int j = 0; j < haloSize; j++)
                    if (halo[i][j] != -1)
                        tmpH[i]++;
            }
            if (sHalo[i][0] != -1)
            {
                for (int j = 0; j < haloSize; j++)
                    if (sHalo[i][j] != -1)
                        tmpSH[i]++;
            }
        }


        int **haloOpt = NULL, **sHaloOpt = NULL;
        haloOpt = new int* [haloOptSize];
        sHaloOpt = new int* [haloOptSize];


        int curH = 0, curSH = 0;
        for (int i = 0; i < nProc; i++)
        {
            if (halo[i][0] != -1)
            {
                haloOpt[curH] = new int [tmpH[i] + 2];
                haloOpt[curH][0] = i;
                haloOpt[curH][1] = tmpH[i];
                for (int j = 0; j < tmpH[i]; j++)
                    haloOpt[curH][j + 2] = halo[i][j];
                curH++;
            }

            if (sHalo[i][0] != -1)
            {
                sHaloOpt[curSH] = new int [tmpSH[i] + 2];
                sHaloOpt[curSH][0] = i;
                sHaloOpt[curSH][1] = tmpSH[i];
                for (int j = 0; j < tmpSH[i]; j++)
                    sHaloOpt[curSH][j + 2] = sHalo[i][j];
                curSH++;
            }
        }

        // update opt halo by g2l
        for (int i = 0; i < haloOptSize; i++)
        {
            for (int j = 0; j < haloOpt[i][1]; j++)
                haloOpt[i][j + 2] = global2loc[haloOpt[i][j + 2]];
            for (int j = 0; j < sHaloOpt[i][1]; j++)
                sHaloOpt[i][j + 2] = global2loc[sHaloOpt[i][j + 2]];
        }

        // cout << Nx << ' ' << Ny << ' ' << Nz << ' ' << N << endl;

        long double seqTimeDot = 0.0;
        long double seqTimeAxpby = 0.0;
        long double seqTimeSpmv = 0.0;

        for (int i = 0; i < 6; i++)
        {
            int *rows = new int [rowsSize];
            Matrix testMat(startI, startJ, startK, Nxp, Nyp, Nzp, Nx, Ny, Nz, rows);
            Vector testVec1(rowsSize, haloRedSize, rows), testVec2(rowsSize, haloRedSize, rows);
            Vector testRes(rowsSize, haloRedSize, 0.0);

            int num = thread[i];

            omp_set_num_threads(num);
            if (myRank == 0)
                cout << "> Threads = " << num << endl << endl;

            long double operations = 1e-9 * 2 * N;
            long double testTimeR = 0.0;
            long double testTime = omp_get_wtime();
            long double dotRes = dot(testVec1, testVec2);
            testTime = omp_get_wtime() - testTime;
            MPI_Reduce(&testTime, &testTimeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myRank == 0)
            {
                cout << "> DOT = " << dotRes << endl;
                cout << "> Time of DOT = " << testTimeR << endl;
                cout << "> GFLOPS = " << operations / testTimeR << endl;
                if (i == 0)
                    seqTimeDot = testTimeR;
                cout << "> SpeedUp = " << seqTimeDot / testTimeR << endl;
                cout << endl;
            }


            operations = 1e-9 * 3 * N;
            testTime = omp_get_wtime();
            axpby(testVec1, testVec2, 1.0, 2.0);
            testTime = omp_get_wtime() - testTime;
            MPI_Reduce(&testTime, &testTimeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myRank == 0)
            {
                cout << "> AXPBY L2 norm = " << sqrt(dot(testVec1, testVec1)) << endl;
                cout << "> Time of AXPBY = " << testTimeR << endl;
                cout << "> GFLOPS = " << operations / testTimeR << endl;
                if (i == 0)
                    seqTimeAxpby = testTimeR;
                cout << "> SpeedUp = " << seqTimeAxpby / testTimeR << endl;
                cout << endl;
            }


            operations = 1e-9 * 2 * N * 7;
            testTime = omp_get_wtime();
            sync(testVec2, halo, sHalo, haloOptSize);
            SpMV(testMat, testVec2, testVec1);
            testTime = omp_get_wtime() - testTime;
            MPI_Reduce(&testTime, &testTimeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myRank == 0)
            {
                cout << "> SpMV L2 norm = " << sqrt(dot(testVec1, testVec1)) << endl;
                cout << "> Time of SpMV = " << testTimeR << endl;
                cout << "> GFLOPS = " << operations / testTimeR << endl;
                if (i == 0)
                    seqTimeSpmv = testTimeR;
                cout << "> SpeedUp = " << seqTimeSpmv / testTimeR << endl;
                cout << endl;
            }


            if ((k == 4) && (i == 5))
            {
                Matrix A(Nx, Ny, Nz, 1, 1, 1, 1, 1, 1, NULL);
                Vector BB(N);

                if (myRank == 0)
                    cout << "> Solver test..." << endl << endl;
                long double seqTimeSolve = 0.0;
                for (int j = 0; j < 6; j++)
                {
                    omp_set_num_threads(thread[j]);
                    if (myRank == 0)
                        cout << "> Number of threads = " << thread[j] << endl;
                    // 5 dot 6 axpby 4 spmv N diag
                    long double operations = 1e-9 * (5 * (2 * N) + 6 * (3 * N) + 4 * (2 * N * 7) + N);
                    long double time = omp_get_wtime();
                    int res = solve(N, A, BB, 0.000000001, 1000, 1, halo, sHalo, haloOptSize, rowsSize, haloRedSize, rows);
                    time = omp_get_wtime() - time;
                    MPI_Reduce(&time, &testTimeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    if (myRank == 0)
                    {
                        cout << "> Numder of iters = " << res << endl;
                        cout << "> Final time of computation = " << testTimeR << endl;
                        cout << "> GFLOPS = " << operations / testTimeR << endl;
                        if (j == 0)
                            seqTimeSolve = testTimeR;
                        cout << "> SpeedUp = " << seqTimeSolve / testTimeR << endl;
                        cout << endl;
                    }
                }
            }
            delete [] rows;
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
        MPI_Finalize();
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
        MPI_Finalize();
        return -1;
    }
    // cout << Nx << ' ' << Ny << ' ' << Nz << ' ' << tol << ' ' << maxit << endl;

    // sizes count
    int kk = myRank / (Px * Py),
        jj = (myRank % (Px * Py)) / Px,
        ii = myRank - kk * (Px * Py) - jj * Px;

    int Nxp = int(ceil(double(Nx) / double(Px))),
        Nyp = int(ceil(double(Ny) / double(Py))),
        Nzp = int(ceil(double(Nz) / double(Pz)));

    int startI = ii * Nxp,
        startJ = jj * Nyp,
        startK = kk * Nzp;

    if (ii == Px - 1)
        Nxp = Nx - ii * Nxp;
    if (jj == Py - 1)
        Nyp = Ny - jj * Nyp;
    if (kk == Pz - 1)
        Nzp = Nz - kk * Nzp;

    int rowsSize = Nxp * Nyp * Nzp;

    int *rows = new int [rowsSize];
    int *global2loc = new int [N];
    int *part = new int [N];

    int **halo = NULL, **sHalo = NULL;

    for (int i = 0; i < rowsSize; i++)
        rows[i] = -1;

    // part
    for (int i = 0; i < nProc; i++)
    {
        int kk = i / (Px * Py),
            jj = (i % (Px * Py)) / Px,
            ii = i - kk * (Px * Py) - jj * Px;

        int Nxp = int(ceil(double(Nx) / double(Px))),
            Nyp = int(ceil(double(Ny) / double(Py))),
            Nzp = int(ceil(double(Nz) / double(Pz)));

        int startI = ii * Nxp,
            startJ = jj * Nyp,
            startK = kk * Nzp;

        if (ii == Px - 1)
            Nxp = Nx - ii * Nxp;
        if (jj == Py - 1)
            Nyp = Ny - jj * Nyp;
        if (kk == Pz - 1)
            Nzp = Nz - kk * Nzp;

        int Ppx = Nxp + startI,
            Ppy = Nyp + startJ,
            Ppz = Nzp + startK;

        for (int k = startK; k < Ppz; k++)
            for (int j = startJ; j < Ppy; j++)
                for (int l = startI; l < Ppx; l++)
                    part[l + Nx * j + Nx * Ny * k] = i;
    }


    Matrix A(startI, startJ, startK, Nxp, Nyp, Nzp, Nx, Ny, Nz, rows);

    // g2l
    for (int i = 0; i < N; i++)
        global2loc[i] = -1;
    for (int i = 0; i < rowsSize; i++)
        if (rows[i] != -1)
            global2loc[rows[i]] = i;

    halo = new int* [nProc];
    sHalo = new int* [nProc];

    // halo shalo
    int haloSize = 2 * (Nxp * Nyp + Nxp * Nzp + Nyp * Nzp);
    for (int i = 0; i < nProc; i++)
    {
        halo[i] = new int [haloSize];
        sHalo[i] = new int [haloSize];

        for (int j = 0; j < haloSize; j++)
        {
            halo[i][j] = -1;
            sHalo[i][j] = -1;
        }
    }
    makeHalo(A, global2loc, part, rows, halo, sHalo, haloSize);

    int haloRedSize = 6 * haloSize;
    int *haloRed = new int [haloRedSize];

    int curHR = 0;
    for (int i = 0; i < nProc; i++)
        if (halo[i][0] != -1)
            for (int j = 0; j < haloSize; j++)
                if (halo[i][j] != -1)
                    haloRed[curHR++] = halo[i][j];
    haloRedSize = curHR;

    // update g2l
    for (int i = 0; i < haloRedSize; i++)
        global2loc[haloRed[i]] = i + rowsSize;
    changeJA(A, global2loc);

    Vector BB(rowsSize, haloRedSize, rows);


    // optimized halo
    int *tmpH = new int [nProc];
    int *tmpSH = new int [nProc];
    int haloOptSize = 0;

    for (int i = 0; i < nProc; i++)
    {
        tmpH[i] = 0;
        tmpSH[i] = 0;
    }

    for (int i = 0; i < nProc; i++)
    {
        if (halo[i][0] != -1)
        {
            haloOptSize++;
            for (int j = 0; j < haloSize; j++)
                if (halo[i][j] != -1)
                    tmpH[i]++;
        }
        if (sHalo[i][0] != -1)
        {
            for (int j = 0; j < haloSize; j++)
                if (sHalo[i][j] != -1)
                    tmpSH[i]++;
        }
    }


    int **haloOpt = NULL, **sHaloOpt = NULL;
    haloOpt = new int* [haloOptSize];
    sHaloOpt = new int* [haloOptSize];


    int curH = 0, curSH = 0;
    for (int i = 0; i < nProc; i++)
    {
        if (halo[i][0] != -1)
        {
            haloOpt[curH] = new int [tmpH[i] + 2];
            haloOpt[curH][0] = i;
            haloOpt[curH][1] = tmpH[i];
            for (int j = 0; j < tmpH[i]; j++)
                haloOpt[curH][j + 2] = halo[i][j];
            curH++;
        }

        if (sHalo[i][0] != -1)
        {
            sHaloOpt[curSH] = new int [tmpSH[i] + 2];
            sHaloOpt[curSH][0] = i;
            sHaloOpt[curSH][1] = tmpSH[i];
            for (int j = 0; j < tmpSH[i]; j++)
                sHaloOpt[curSH][j + 2] = sHalo[i][j];
            curSH++;
        }
    }

    // update opt halo by g2l
    for (int i = 0; i < haloOptSize; i++)
    {
        for (int j = 0; j < haloOpt[i][1]; j++)
            haloOpt[i][j + 2] = global2loc[haloOpt[i][j + 2]];
        for (int j = 0; j < sHaloOpt[i][1]; j++)
            sHaloOpt[i][j + 2] = global2loc[sHaloOpt[i][j + 2]];
    }


    // print
    // if (myRank == 0)
    {
        // cout << "ROWS" << endl;
        // for (int i = 0; i < rowsSize; i++)
        //     cout << rows[i] << ' ';
        // cout << endl << endl;
        // cout << "G2L" << endl;
        // for (int i = 0; i < N; i++)
        //     cout << global2loc[i] << ' ';
        // cout << endl << endl;
        // cout << "PART" << endl;
        // for (int i = 0; i < N; i++)
        //     cout << part[i] << ' ';
        // cout << endl;

        // A.printCSR();
        // cout << endl;
        // BB.print();

        // cout << "HALO" << endl;
        // for (int i = 0; i < nProc; i++)
        // {
        //     for (int j = 0; j < haloSize; j++)
        //         if (halo[i][j] != -1)
        //         {
        //             cout << i << "    ";
        //             for (int j = 0; j < haloSize; j++)
        //                  if (halo[i][j] != -1)
        //                     cout << halo[i][j] << ' ';
        //             cout << endl;
        //             break;
        //         }
        // }

        // cout << "SHALO" << endl;
        // for (int i = 0; i < nProc; i++)
        // {
        //     for (int j = 0; j < haloSize; j++)
        //         if (sHalo[i][j] != -1)
        //         {
        //             cout << i << "    ";
        //             for (int j = 0; j < haloSize; j++)
        //                  if (halo[i][j] != -1)
        //                     cout << sHalo[i][j] << ' ';
        //             cout << endl;
        //             break;
        //         }
        // }
        // cout << endl;

        // cout << "HALOOPT" << endl;
        // for (int i = 0; i < haloOptSize; i++)
        // {
        //     cout << haloOpt[i][0] << "    ";
        //     for (int j = 0; j < haloOpt[i][1]; j++)
        //         cout << haloOpt[i][j + 2] << ' ';
        //     cout << endl;
        // }
        // cout << endl;
        // cout << "SHALOOPT" << endl;
        // for (int i = 0; i < haloOptSize; i++)
        // {
        //     cout << sHaloOpt[i][0] << "    ";
        //     for (int j = 0; j < sHaloOpt[i][1]; j++)
        //         cout << sHaloOpt[i][j + 2] << ' ';
        //     cout << endl;
        // }

        // cout << "ALLHALO" << endl;
        // for (int i = 0; i < haloRedSize; i++)
        //     cout << haloRed[i] << " ";
        // cout << endl;

        // cout << haloSize << endl;
    }


    if (debug != 0)
        testFunc(Px, Py, Pz, myRank, nProc);
    else
    {
        omp_set_num_threads(atoi(argv[6]));
        // 5 dot 6 axpby 4 spmv N diag
        long double operations = 1e-9 * (5 * (2 * N) + 6 * (3 * N) + 4 * (2 * N * 7) + N);
        if (myRank == 0)
            cout << "> Number of threads = " << atoi(argv[6]) << endl;
        long double time = omp_get_wtime();
        int res = solve(N, A, BB, tol, maxit, debug, haloOpt, sHaloOpt, haloOptSize, rowsSize, haloRedSize, rows);
        time = omp_get_wtime() - time;
        long double timeR = 0.0;
        MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (myRank == 0)
        {
            cout << "> Number of iters = " << res << endl;
            cout << "> Final time of computation = " << timeR << endl;
            cout << "> GFLOPS = " << operations / timeR << endl;
        }
    }
    if (myRank == 0)
        cout << "> Operation time = " << globtime << endl;

    delete [] rows;
    delete [] global2loc;

    for (int i = 0; i < nProc; i++)
    {
        delete [] halo[i];
        delete [] sHalo[i];
    }
    delete [] part;
    delete [] halo;
    delete [] sHalo;
    delete [] haloRed;
    for (int i = 0; i < haloOptSize; i++)
    {
        delete [] haloOpt[i];
        delete [] sHaloOpt[i];
    }
    delete [] haloOpt;
    delete [] sHaloOpt;
    delete [] tmpH;
    delete [] tmpSH;
    MPI_Finalize();
    return 0;
}
