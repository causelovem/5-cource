#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

class Vector;

class Matrix
{
    private:
        double *A;
        int *JA;
        int *IA;
        int sizeIA;
        int sizeA;

    public:
        Matrix(int Nx, int Ny, int Nz)
        {
            sizeIA = Nx * Ny * Nz;
            int corner = 4 * 8; // 3 soseda
            int edge = 5 * 4 * ((Nx - 2) + (Ny - 2) + (Nz - 2)); // 4 soseda
            int face = 6 * 2 * (((Nx - 2) * (Ny - 2)) + ((Nx - 2) * (Nz - 2)) + ((Ny - 2) * (Nz - 2))); // Xy Xz Yz
            int inner = 7 * (Nx - 2) * (Ny - 2) * (Nz - 2); // 5 sosedey
            sizeA = corner + edge + face + inner;

            IA = new int [sizeIA];
            JA = new int [sizeA];
            A = new double [sizeA];
            IA[0] = 0;

            int ia = 0, ija = 0, iia = 1;
            int cnt = 0;

            for (int k = 0; k < Nz; k++)
                for (int j = 0; j < Ny; j++)
                    for (int i = 0; i < Nx; i++)
                    {
                        int num = i + Nx * j + Nx * Ny * k;
                        double sum = 0;
                        double s = sin(i + j + 1);

                        if (k > 0)
                        {
                            A[ia] = s;
                            sum += A[ia++];
                            JA[ija++] = num - Nx * Ny;
                            cnt++;
                        }

                        if (j > 0)
                        {
                            A[ia] = s;
                            sum += A[ia++];
                            JA[ija++] = num - Nx;
                            cnt++;
                        }

                        if (i > 0)
                        {
                            A[ia] = s;
                            sum += A[ia++];
                            JA[ija++] = num - 1;
                            cnt++;
                        }

                        int tia = ia++, tija = ija++;

                        if (i < Nx - 1)
                        {
                            A[ia] = s;
                            sum += A[ia++];
                            JA[ija++] = num + 1;
                            cnt++;
                        }

                        if (j < Ny - 1)
                        {
                            A[ia] = s;
                            sum += A[ia++];
                            JA[ija++] = num + Nx;
                            cnt++;
                        }

                        if (k < Nz - 1)
                        {
                            A[ia] = s;
                            sum += A[ia++];
                            JA[ija++] = num + Nx * Ny;
                            cnt++;
                        }


                        A[tia] = sum * 1.1;
                        JA[tija] = num;
                        IA[iia++] = ++cnt;
                    }
        }

        ~Matrix()
        {
            delete [] A;
            delete [] JA;
            delete [] IA;
        }

        double get(int i, int j) const
        {
            for (int k = IA[i]; k < IA[i + 1]; k++)
                if (j == JA[k])
                    return A[k];
            return 0;
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
        }

        void print(int flg = 0)
        {
            if (flg == 0)
                printCSR();

            cout << "Matrix" << endl;
            for (int i = 0; i < sizeIA; i++)
            {
                for (int j = 0; j < sizeIA; j++)
                    cout << get(i, j) << ' ';
                cout << endl;
            }
        }

        friend void SpMV(const Matrix &mat, const Vector &vec, Vector &res);
};

class Vector
{
    private:
        double *A;
        int size;
    public:
        Vector(int s)
        {
            size = s;
            A = new double [s];

            for (int i = 0; i < size; i++)
                A[i] = rand() % 100;
        }

        Vector(int s, double c)
        {
            size = s;
            A = new double [s];

            for (int i = 0; i < size; i++)
                A[i] = c;
        }

        ~Vector()
        {
            delete [] A;
        }

        Vector operator = (const Vector &vec)
        {
            if (A != NULL)
                delete [] A;

            size = vec.size;
            A = new double [size];

            for (int i = 0; i < size; i++)
                A[i] = vec.A[i];

            return *this;
        }

        Vector operator = (const double c)
        {
            for (int i = 0; i < size; i++)
                A[i] = c;

            return *this;
        }

        void print()
        {
            cout << "Vector" << endl;
            for (int i = 0; i < size; i++)
                cout << A[i] << ' ';
            cout << endl;

            return;
        }
    
        friend double dot(const Vector &vec1, const Vector &vec2)
        {
            if (vec1.size != vec2.size)
                cout << "Can't dot two vectors: different lenghts!" << endl;

            double res = 0;

            for (int i = 0; i < vec1.size; i++)
                res += vec1.A[i] * vec2.A[i];

            return res;
        }

        friend void axpby(Vector &vec1, const Vector &vec2, double a, double b)
        {
            if (vec1.size != vec2.size)
                cout << "Can't axpby two vectors: different lenghts!" << endl;

            for (int i = 0; i < vec1.size; i++)
                vec1.A[i] = a * vec1.A[i] + b * vec2.A[i];

            return;
        }

        friend void SpMV(const Matrix &mat, const Vector &vec, Vector &res)
        {
            if (mat.sizeIA != vec.size)
                cout << "Can't SpMV: different lenghts!" << endl;

            res = 0;
            for (int i = 0; i < vec.size; i++)
                for (int j = 0; j < vec.size; j++)
                    res.A[i] += mat.get(i, j) * vec.A[j];

            return;
        }
};

int main ()
{
    srand(time(0));

    // Matrix A(2, 2, 2);
    // // A.print();

    // Vector B(10);
    // Vector C(10, 1);

    // B.print();
    // C.print();

    // cout << dot(B, C) << endl;

    // axpby(B, C, 1, 2);
    // B.print();

    // Vector R(8, 0);
    // Vector V(8);

    // SpMV(A, V, R);

    // R.print();

    Matrix A(2, 2, 2);
    Vector BB(8);

    

    return 0;
}