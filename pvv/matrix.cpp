#include <iostream>
#include <cmath>

using namespace std;

class Matrix
{
private:
    double *A;
    int *JA;
    int *IA;
    int sizeIA;
    int sizeA;

public:
    Matrix (int Nx, int Ny, int Nz, int rank)
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

    ~Matrix ()
    {
        delete [] A;
        delete [] JA;
        delete [] IA;
    }

    double get(int i, int j)
    {
        for (int k = IA[i]; k < IA[i + 1]; k++)
            if (j == JA[k])
                return A[k];
        return 0;
    }

    void print()
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

        cout << 'Matrix' << endl;
        for (int i = 0; i < sizeIA; i++)
        {
            for (int j = 0; j < sizeIA; j++)
                cout << get(i, j) << ' ';
            cout << endl;
        }
    }
};

int main ()
{
    Matrix A(3, 2, 2, 1);
    A.print();
}