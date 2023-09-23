#include "defclique/defclique.h"
#include "utils/bigraph.hpp"
#include "utils/log.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <vector>

void usage() {
	std::cout << "Usage: bin/run <dataset> <k>" << std::endl;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		usage();
		return 0;
	}

	Graph G(argv[1]);

	auto startTimePoint = std::chrono::steady_clock::now();

	defclique::russianDollSearch(G, atoi(argv[2]));

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);

	log("Total time: %ld ms\n", duration.count());

	return 0;
}
