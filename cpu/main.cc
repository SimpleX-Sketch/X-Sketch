#include <bits/stdc++.h>
#include <Benchmark.h>
using namespace std;

int main(int argc, char** argv) {
	init_matrix();
	// string PATH = "/share/datasets/webpage/webdocs_form01.dat"; 
	string PATH = argv[1];
	string RES = argv[2]; 
	CAIDABenchmark benchmark(PATH,RES);
	// double var_t = atoi(argv[1]);
	benchmark.Run();
	return 0;
}