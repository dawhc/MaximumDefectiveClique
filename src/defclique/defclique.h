#ifndef DEFCLIQUE_H
#define DEFCLIQUE_H

#include "../utils/graph.hpp"
#include "../utils/vertexset.hpp"
#include "../utils/ordering.hpp"
#include <cstdint>
#include <string>

#define RUSSIANDOLL_SEARCH 0
#define REDUCTION_SEARCH 1
#define ONE_HOP 0
#define TWO_HOP 1

namespace defclique {
	void logSet(VertexSet &V, const std::string &name);
	Graph coreReduction(Graph &G, int k);
	Graph edgeReduction(Graph &G, int k);
	void preprocessing(Graph &G, Ordering &o, int u, int mode=TWO_HOP);
	VertexSet heuristic(Graph &G, int k);
	int upperbound();
	void run(Graph &G, int k, int mode=REDUCTION_SEARCH);
	void moveCToS(int v);
	void moveSToC(int v);
	int updateC(int v);
	void restoreC(int pos);
	int update(int v);
	void restore(int v, int posC);
	bool branch(int dep);
	void add(Graph &G, VertexSet &V, std::vector<int> &degV, int v);
	void sub(Graph &G, VertexSet &V, std::vector<int> &degV, int v);
	void addC(int v);
	void subC(int v);
}

#endif // DEFCLIQUE_H