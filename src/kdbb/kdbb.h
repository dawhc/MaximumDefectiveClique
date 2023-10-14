#ifndef KDBB_H
#define KDBB_H

#pragma once

#include "../utils/graph.hpp"
#include "../utils/vertexset.hpp"


namespace kdbb {
	int upperbound();
	VertexSet fastLB();
	Graph preprocessing(Graph &G, int k, int lb);
	Graph coreReduction(Graph &G, int k);
	Graph edgeReduction(Graph &G, int k);
}



#endif