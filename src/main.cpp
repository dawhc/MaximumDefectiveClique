#include "defclique/defclique.h"
#include "kdbb/kdbb.h"
#include "utils/bigraph.hpp"
#include "utils/log.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <vector>


#define ALGO_DEFAULT 0
#define ALGO_KDBB 1

void usage() {
	std::cout << "Usage: bin/run <dataset> <k> [RussianDoll/Reduction/KDBB]" << std::endl;
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		usage();
		return 0;
	}

	int mode = REDUCTION_SEARCH;
	int algo = ALGO_DEFAULT;
	int k = atoi(argv[2]);

	if (argc == 3) mode = REDUCTION_SEARCH;
	else if (!strcmp("RussianDoll", argv[3])) mode = RUSSIANDOLL_SEARCH;
	else if (!strcmp("Reduction", argv[3])) mode = REDUCTION_SEARCH;
	else if (!strcmp("KDBB", argv[3])) {
		algo = ALGO_KDBB;
	}
	else {
		usage();
		return 0;
	}
	

	auto startTimePoint = std::chrono::steady_clock::now();

	if (algo == ALGO_DEFAULT)
		defclique::run(argv[1], k, mode);
	else if (algo == ALGO_KDBB)
		kdbb::run(argv[1], k);

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);

	log("Total time spent: %ld ms\n", duration.count());

	return 0;
}
