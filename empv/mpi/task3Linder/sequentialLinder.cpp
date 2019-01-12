#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>

using namespace std;

string update_data(string data, map<char, string>& R)
{
    string buf = "";

    for( unsigned int i=0; i<data.length(); i++ )
        buf += R[data[i]];

    return buf;
}

void run_lsystem(int T)
{
    string data = "a";

    map<char, string> R;
    R['a'] = "b";
    R['b'] = "ab";

    for( int t=0; t<T; t++ )
        data = update_data(data,R);

    ofstream f("output.dat");
    f << data << endl;
    f.close();
}

int main(int argc, char** argv)
{
    int T = atoi(argv[1]);

    run_lsystem(T);

    return 0;
}