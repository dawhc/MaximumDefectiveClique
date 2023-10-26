#include "defclique/defclique.h"
#include "utils/bigraph.hpp"
#include "utils/log.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <vector>

void usage() {
	std::cout << "Usage: bin/run <dataset> <k> [RussianDoll/Reduction]" << std::endl;
}

int main(int argc, char* argv[]) {

	if (argc < 3) {
		usage();
		return 0;
	}

	int mode = REDUCTION_SEARCH;

	if (argc >= 4 && !strcmp("RussianDoll", argv[3]))
		mode = RUSSIANDOLL_SEARCH;
	else if (argc >= 4 && !strcmp("Reduction", argv[3]) || argc == 3)
		mode = REDUCTION_SEARCH;
	else {
		usage();
		return 0;
	}

	log("Reading graph: %s ...", strrchr(argv[1], '/')+1);

	auto startTimePoint = std::chrono::steady_clock::now();

	Graph G(argv[1]);

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);

	log("Reading graph done! Time spent: %ld ms", duration.count());
	log("Graph info: n=%d, m=%d, maxdeg=%d", G.n, G.m, G.maxDeg);

	startTimePoint = std::chrono::steady_clock::now();

	defclique::run(G, atoi(argv[2]), mode);

	duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);

	log("Total time spent: %ld ms\n", duration.count());

	return 0;
}
