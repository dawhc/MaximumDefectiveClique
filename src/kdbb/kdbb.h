#ifndef KDBB_H
#define KDBB_H

#pragma once

#include "../utils/graph.hpp"
#include "../utils/vertexset.hpp"
#include <string>
#include <vector>


namespace kdbb {
	int upperbound();
	int fastLB(std::string filename);
	Graph preprocessing(Graph &G, int k, int lb);
	Graph coreReduction(Graph &G, int k);
	Graph edgeReduction(Graph &G, int k);
	int KDBB(std::string filename, int k);
	int branch(int v, int lb);
}



#endif