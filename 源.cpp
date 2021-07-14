//#include "eigen3/Eigen/Dense"
#include <iostream>
#include "include/Pipe.h"
static const int N = 100;

int main() {
	Pipe p[100];
	ifstream file;
	for (int i = 0; i < N; i++)
	{
		//file.open("pipe"+to_string(i)+".txt");
		file.open("pipe.txt");
		p[i].input(file);
		file.close();
		ofstream ofile;
		ofile.open("pipe" + to_string(i) + ".out");
		p[i].output(ofile);
	}
	
	Tdouble P[N] = {};
	Tdouble W = 2.3;
	for(int i = 0;i<N;i++)
	{	
		P[i] = 100.0 * i + 10000.0;
		p[i].boundary.Pin = &P[i];
		if (i != N-1)
			p[i].boundary.Pout = &P[i + 1];
		else
			p[i].boundary.Pout = &P[0];
		p[i].boundary.W = &W;
	}
	Tdouble err = 1.0;
	Tdouble iter = P[0];
	while (abs(err) > 0.1) {
		for (int i = 0; i < N; i++)
		{
			p[i].steady();
		}
		err = P[0] - iter;
		iter = P[0];
		cout << err <<" "<<iter<< endl;
	}


	for (int i = 0; i < N; i++)
	{
		cout << *(p[i].boundary.Pout) << endl;
	}

	return 0;
}
