#ifndef DEFCLIQUE_H
#define DEFCLIQUE_H

#include "../utils/graph.hpp"
#include "../utils/vertexset.hpp"
#include "../utils/ordering.hpp"
#include <cstdint>
#include <string>

namespace defclique {
	void logSet(VertexSet &V, const std::string &name);
	Graph coreReduction(Graph &G, int k);
	Graph edgeReduction(Graph &G, int k);
	Graph subGraph(Graph &G, VertexSet V);
	VertexSet heuristic(Graph &G, int k);
	int upperbound();
	void russianDollSearch(Graph G, int k);
	void reductionSearch(Graph G, int k);
	void moveCToS(int v);
	void moveSToC(int v);
	int updateC(int v);
	void restoreC(int pos);
	bool branch(int dep);
	void add(Graph &G, VertexSet &V, std::vector<int> &degV, int v);
	void sub(Graph &G, VertexSet &V, std::vector<int> &degV, int v);
}

#endif // DEFCLIQUE_H