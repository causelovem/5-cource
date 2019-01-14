#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <mpi.h>

using namespace std;

string update_data(string data, map<char, pair<string, double> >& RR)
{
    string buf = "";

    for( unsigned int i = 0; i < data.length(); i++ )
        if (((double)rand() / RAND_MAX) <= RR[data[i]].second)
            buf += RR[data[i]].first;
        else
            buf += data[i];

    return buf;
}

int get_count(int myL, int neibL)
{
    if (myL > neibL)
        return int((myL - neibL) / 2);
    return 0;
}

string align_load(string w, int myRank, int nProc)
{
    int prev = myRank - 1;
    int next = myRank + 1;

    if (prev == -1)
        prev = MPI_PROC_NULL;

    if (next == nProc)
        next = MPI_PROC_NULL;


    int myL = w.length(), nextL = w.length(), prevL = 0;
    int prevC = 0, nextC = 0;

    MPI_Sendrecv(&myL, 1, MPI_INT, prev, 0, &nextL, 1, MPI_INT, next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int c = get_count(myL, nextL);
    MPI_Sendrecv(&c, 1, MPI_INT, next, 0, &prevC, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    char *prevW = new char [prevC];
    MPI_Sendrecv((void *)w.substr(myL - c, myL).c_str(), c, MPI_CHAR, next, 0, prevW, prevC, MPI_CHAR, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    w = string(prevW, prevC) + w.substr(0, myL - c);

    myL = prevL = w.length();

    MPI_Sendrecv(&myL, 1, MPI_INT, next, 0, &prevL, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    c = get_count(myL, prevL);
    MPI_Sendrecv(&c, 1, MPI_INT, prev, 0, &nextC, 1, MPI_INT, next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    char *nextW = new char [nextC];
    MPI_Sendrecv((void *)w.substr(0, c).c_str(), c, MPI_CHAR, prev, 0, nextW, nextC, MPI_CHAR, next, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    w = w.substr(c, myL) + string(nextW, nextC);

    delete[] prevW;
    delete[] nextW;

    return w;
}

void run_lsystem(int T, int k, int myRank, int nProc)
{
    int lsize = T / k;
    int *load = new int [lsize];
    lsize = 0;
    string data;

    map< char, pair<string, double> > RR;

    // RR['a'] = make_pair("b", 1);
    // RR['b'] = make_pair("ab", 1);

    // RR['a'] = make_pair("ab", 1);
    // RR['b'] = make_pair("bc", 1);
    // RR['c'] = make_pair("c", 1);

    // RR['a'] = make_pair("aa", 0.001);

    RR['a'] = make_pair("ab", 0.01);
    RR['b'] = make_pair("a", 0.01);

    int center = nProc / 2;
    if (myRank == center)
        data = "a";
    else
        data = "";

    long double time = MPI_Wtime();
    for( int t = 1; t < T + 1; t++ )
    {
        data = update_data(data, RR);

        if (t % k == 0)
        {
            load[lsize++] = data.length();
            data = align_load(data, myRank, nProc);
        }
    }
    time = MPI_Wtime() - time;
    long double timeR = 0.0;
    MPI_Reduce(&time, &timeR, 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myRank == 0)
    {
        lsize = T / k;
        int **allLoad = new int* [nProc];

        for (int i = 0; i < nProc; i++)
            allLoad[i] = new int [lsize];
        memcpy(allLoad[0], load, lsize * sizeof(int));

        for (int i = 1; i < nProc; i++)
        {
            int len = 0;
            MPI_Recv(&len, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            char *neibW = new char [len];
            MPI_Recv(neibW, len, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            data += string(neibW, len);
            delete[] neibW;

            len = T / k;
            MPI_Recv(allLoad[i], lsize, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // ofstream fo("output.dat");
        // fo << data << endl;
        // fo.close();

        ofstream fs("stat.dat", fstream::app);

        fs << timeR << " " << nProc << " " << T << " " << k << endl;
        for (map< char, pair<string, double> >::const_iterator it = RR.begin(); it != RR.end(); it++)
            fs << "> " << it->first << "  " << it->second.first << "  " << it->second.second << endl;

        int sum = 0;
        for (int i = 0; i < lsize; i++)
        {
            sum = 0;
            for (int j = 0; j < nProc; j++)
                sum += allLoad[j][i];

            fs << (i + 1) * k << "   ";
            for (int j = 0; j < nProc; j++)
                fs << double(allLoad[j][i]) / double(sum) << ' ';
            fs << endl;
        }

        fs << data.length() << endl;

        for (int i = 0; i < nProc; i++)
            delete[] allLoad[i];
        delete[] allLoad;
        fs << endl << endl;
        fs.close();
    }
    else
    {
        int len = data.length();
        MPI_Send(&len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send((void *)data.c_str(), data.length(), MPI_CHAR, 0, 0, MPI_COMM_WORLD);

        MPI_Send(load, lsize, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    delete[] load;

    return;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int nProc = 0, myRank = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    int T = atoi(argv[1]);
    int k = atoi(argv[2]);

    srand(time(NULL) + myRank);

    run_lsystem(T, k, myRank, nProc);

    MPI_Finalize();
    return 0;
}
