#include "defclique/defclique.h"
#include "kdbb/kdbb.h"
#include "utils/bigraph.hpp"
#include "utils/log.hpp"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "utils/cmdline.hpp"


int main(int argc, char* argv[]) {

	cmdline::parser args;

	args.add<std::string>("data", 'd', "dataset path", true, "");
	args.add<int>("key", 'k', "value of k", true, 1);
	args.add<std::string>("algo", 'a', "algorithm", false, "MDC", cmdline::oneof<std::string>("MDC", "RussianDoll", "KDBB", "PMC"));

	args.parse_check(argc, argv);

	auto dataPath = args.get<std::string>("data");
	auto k = args.get<int>("key");
	auto algo = args.get<std::string>("algo");

	auto startTimePoint = std::chrono::steady_clock::now();

	if (algo == "MDC") defclique::run(dataPath, k, REDUCTION_SEARCH);
	else if (algo == "RussianDoll") defclique::run(dataPath, k, RUSSIANDOLL_SEARCH);
	else if (algo == "KDBB") kdbb::run(dataPath, k);
	else if (algo == "PMC") kdbb::fastLB(dataPath);

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);

	log("Total time spent: %ld ms", duration.count());
	

	return 0;
}
