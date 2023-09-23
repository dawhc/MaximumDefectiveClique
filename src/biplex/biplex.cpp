#include "biplex.h"
#include "linearheap.hpp"
#include "graph.hpp"
#include "vertexset.hpp"
#include <inttypes.h>
#include <chrono>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>


// #define DEBUG
#define UPPERBOUND_PRUNE
#define CALC_DEGREE
#define SMALL_K 2
// #define DESTRUCT_SETS // accelerate dataset "dblp"
// #define PRUNE_C
#define PRUNE_ONCE
// #define PRUNE_WITH_QUEUE
#define WITH_PIVOTING
// #define ENUM_LARGER_DEGREE_SIDE
// #define NO_ORDERING
#define CHOOSE_FIRST_VERTEX

#define STRINGIFY(name) #name

#define PRINT_MACRO_INFO(name) do { \
	if (#name [0] != STRINGIFY(name) [0]) \
		printf("%s: on\n", #name); \
	else \
		printf("%s: off\n", #name); \
} while (0)

namespace biplex {
	BiGraph G;
	VertexSet S[2], C[2], X[2], Q[2];
	std::vector<std::vector<uint32_t>> cand(100000);
	std::vector<std::vector<uint32_t>> nonNbrS[2];
	std::vector<uint32_t> nonNbrC;
	std::vector<uint32_t> degC[2];
	std::vector<uint32_t> ordered[2], order[2];
	std::vector<uint32_t> sup, ptrB;
	std::vector<std::vector<uint32_t>> B;
	uint32_t q, k, lb;
	Result result;
	uint64_t resultNumThres;
	std::chrono::time_point<std::chrono::steady_clock> startTime, currTime;
	uint64_t pivotingTime, bipartiteTime, calcTime;

	void setTime() {
		currTime = std::chrono::steady_clock::now();
	}

	void getTime(uint64_t &timeCount) {
		timeCount += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - currTime).count();
	}
}

/* Heuristic algorithm */
// uint32_t biplex::biplexLowerBound(const BiGraph& G) 
// {

// 	LinearHeap<uint32_t> vHeap[2] = {
// 		LinearHeap<uint32_t>(G.degree[0]), 
// 		LinearHeap<uint32_t>(G.degree[1])
// 	};

// 	while(!vHeap[0].empty() && !vHeap[1].empty()) {
// 		uint32_t side = -1u, u = -1u;
// 		for (uint32_t i = 0; i <= 1; ++i) {
// 			uint32_t v = vHeap[i].top();
// 			if (vHeap[i][v] >= vHeap[i^1].size() - k) continue;
// 			if (side == -1u || vHeap[i][v] < vHeap[side][u]) {
// 				u = v; side = i;
// 			}
// 		}
// 		if (side == -1u) break;
// 		vHeap[side].pop();
// 		for (uint32_t v : G.nbr[side][u])
// 			if (vHeap[side^1].inside(v))
// 				vHeap[side^1].dec(v);
// 	}

// 	return std::min(vHeap[0].size(), vHeap[1].size());
// }



/**
 * Core ordering, considering bigraph as a normal graph 
 */
void biplex::coreOrdering(const BiGraph& G) {
	Graph H(G.n[0] + G.n[1]);
	for (uint32_t u = 0; u < G.n[0]; ++u) 
		for (uint32_t v : G.nbr[0][u])
			H.addEdge(u, v+G.n[0]);

	LinearHeap<uint32_t> vHeap(H.degree);

	for (uint32_t side = 0; side <= 1; ++side) {
		order[side].resize(G.n[side]);
		ordered[side].resize(G.n[side]);
	}

	uint32_t nOrdered[2] = {0};

	while (!vHeap.empty()) {
		uint32_t u = vHeap.top(); vHeap.pop();
		for (uint32_t v : H.nbr[u]) {
			if (vHeap[v] <= vHeap[u]) continue;
			if (vHeap.inside(v)) vHeap.dec(v);
		}

		uint32_t side = 0;
		if (u >= G.n[0]) {
			side ^= 1;
			u -= G.n[0];
		} 

		ordered[side][nOrdered[side]] = u;
		order[side][u] = nOrdered[side]++;
	}

}

/**
 * Generate (alpha, beta)-core of a bigraph
 */
BiGraph biplex::coreReduction(const BiGraph& G, uint32_t alpha, uint32_t beta) {
	std::queue<uint32_t> q[2];
	std::vector<uint32_t> deg[2] = {G.degree[0], G.degree[1]};
	std::vector<bool> vis[2] = {
		std::vector<bool>(G.n[0]), 
		std::vector<bool>(G.n[1])
	};
	uint32_t coreNum[2] = {alpha, beta};
	uint32_t n[2] = {G.n[0], G.n[1]};
	for (uint32_t side = 0; side <= 1; ++side)
		for (uint32_t u = 0; u < G.n[side]; ++u)
			if (deg[side][u] < coreNum[side^1]) {
				vis[side][u] = true;
				q[side].push(u);
			}


	while (!q[0].empty() || !q[1].empty()) {
		for (uint32_t side = 0; side <= 1; ++side) {
			while (!q[side].empty()) {
				uint32_t u = q[side].front();
				q[side].pop();
				// Remove u from side
				--n[side];
				for (uint32_t v : G.nbr[side][u]) {
					if (vis[side^1][v]) continue;
					if (--deg[side^1][v] < coreNum[side]) {
						vis[side^1][v] = true;
						q[side^1].push(v);
					}
				}
			}
		}
	}

	BiGraph C;

	for (uint32_t u = 0; u < G.n[0]; ++u) {
		if (vis[0][u]) continue;
		for (uint32_t v : G.nbr[0][u]) {
			if (vis[1][v]) continue;
			C.addEdgeWithLabel(G.vLabel[0][u], G.vLabel[1][v]);
		}
	}

	return C;
}

/**
 * Butterfly-based reduction
 */
BiGraph biplex::butterflyReduction(const BiGraph& G, uint32_t threshold) {
	std::vector<std::pair<uint32_t, uint32_t>> edges;
	std::vector<std::unordered_map<uint32_t, uint32_t>> eid(G.n[0]);

	for (uint32_t u = 0; u < G.n[0]; ++u)
		for (uint32_t v : G.nbr[0][u]) {
			edges.push_back(std::make_pair(u, v));
			eid[u][v] = edges.size() - 1;
		}

	std::vector<uint32_t>(edges.size()).swap(sup);

	for (uint32_t i = 0; i < edges.size(); ++i) {
		uint32_t u = edges[i].first, v = edges[i].second;
		for (uint32_t w : G.nbr[1][v]) {
			if (w == u) continue;
			for (uint32_t x : G.nbr[0][w]) {
				if (x == v || !G.connect(0, u, x)) continue;
				++sup[i];
			}
		}
	}

	LinearHeap<uint32_t> supHeap(sup);

	while (!supHeap.empty()) {
		uint32_t minSupEid = supHeap.top();
		if (supHeap[minSupEid] >= threshold) break;
		supHeap.pop();

		uint32_t u = edges[minSupEid].first, v = edges[minSupEid].second;

		for (uint32_t i = 0; i < G.nbr[1][v].size(); ++i) {
			uint32_t w = G.nbr[1][v][i], ewv = eid[w][v];
			if (u == w || !supHeap.inside(ewv)) continue;

			for (uint32_t j = 0; j < G.nbr[0][w].size(); ++j) {

				uint32_t x = G.nbr[0][w][j];
				if (v == x || !G.connect(0, u, x)) continue;
				uint32_t eux = eid[u][x], ewx = eid[w][x];
				if (!supHeap.inside(eux) || !supHeap.inside(ewx)) continue;

				supHeap.dec(eux);
				supHeap.dec(ewv);
				supHeap.dec(ewx);
			}
		}
	}

	BiGraph newG;

	while (!supHeap.empty()) {
		uint32_t t = supHeap.top(); supHeap.pop();	
		newG.addEdgeWithLabel(
			G.vLabel[0][edges[t].first], 
			G.vLabel[1][edges[t].second]
		);
	}

	return newG;
}

bool biplex::pruneCX(uint32_t u, uint32_t uSide) {

	for (uint32_t side = 0; side <= 1; ++side) 
		Q[side].clear();

	for (uint32_t v : C[uSide^1]) 
		if (degC[uSide^1][v] + (uint32_t)G.connect(uSide, u, v) + k < lb)
			Q[uSide^1].pushBack(v);

	for (uint32_t w : C[uSide])
		if (degC[uSide][w] + 2*k < lb) 
			Q[uSide].pushBack(w);

#ifdef PRUNE_WITH_QUEUE
	// Prune C
	while (Q[0].size() != 0 || Q[1].size() != 0) {
		// printf("Q[0].size(): %d, Q[1].size(): %d\n", Q[0].size(), Q[1].size());
		while (Q[uSide^1].size() != 0) {
			uint32_t v = Q[uSide^1][0]; Q[uSide^1].popBack(v);
			C[uSide^1].popBack(v);
			for (uint32_t w : G.nbr[uSide^1][v]) {
				--degC[uSide][w];
				if (C[uSide].inside(w) && degC[uSide][w] + 2*k < lb) {
					Q[uSide].pushBack(w);
				} 
				else if (w == u && degC[uSide][u] + k < lb) {
					return false;
				}
			}
		}
		while (Q[uSide].size() != 0) {
			uint32_t w = Q[uSide][0]; Q[uSide].popBack(w);
			C[uSide].popBack(w);
			for (uint32_t v : G.nbr[uSide][w]) {
				--degC[uSide^1][v];
				if (C[uSide^1].inside(v) && degC[uSide^1][v] + (uint32_t)G.connect(uSide, u, v) + k < lb) {
					Q[uSide^1].pushBack(v);
				}
			}
		}
	}

	// Prune X
	for (uint32_t i = X[uSide].frontPos(); i < X[uSide].backPos(); ++i) {
		uint32_t v = X[uSide][i];
		if (degC[uSide][v] + 2*k < lb) {
			X[uSide].popBack(v); --i;
		}
	}
#else

	while (true) {
		bool flag = false;
		for (uint32_t i = C[uSide^1].frontPos(); i < C[uSide^1].backPos(); ++i) {
			uint32_t v = C[uSide^1][i];
			bool isConn = G.connect(uSide, u, v);
			if (degC[uSide^1][v] + (int32_t)isConn < q - k) {
				flag  = true;
				C[uSide^1].popBack(v); --i;
				for (uint32_t w : G.nbr[uSide^1][v]) --degC[uSide][w];
				if (degC[uSide][u] < q - k) return false;  
			}
		}
		if (flag) continue;
		for (uint32_t i = C[uSide].frontPos(); i < C[uSide].backPos(); ++i) {
			uint32_t v = C[uSide][i];
			if (degC[uSide][v] < q - 2*k) {
				flag = true;
				C[uSide].popBack(v); --i;
				for (uint32_t w : G.nbr[uSide][v]) --degC[uSide^1][w];
			}
		}

		for (uint32_t i = X[uSide].frontPos(); i < X[uSide].backPos(); ++i) {
			uint32_t v = X[uSide][i];
			if (degC[uSide][v] < q - 2*k) {
				flag = true;
				X[uSide].popBack(v); --i;
			}
		}

		if (!flag) break;
		
	}		
#endif
	return true;
}

bool biplex::constructSets(uint32_t u, uint32_t uSide) 
{
	// Clear

#ifndef DESTRUCT_SETS
	for (uint32_t side = 0; side <= 1; ++side) {
		S[side].clear();
		C[side].clear();
		X[side].clear();
		std::vector<std::vector<uint32_t>>(G.n[side]).swap(nonNbrS[side]);
		std::vector<uint32_t>(G.n[side]).swap(degC[side]);
	}
#endif


	// construct S, C and X
	S[uSide].pushBack(u);
	
	for (uint32_t v : G.nbr[uSide][u]) {
		C[uSide^1].pushBack(v);
		for (uint32_t w : G.nbr[uSide^1][v]) {
			if (w == u) continue;
			if (order[uSide][w] > order[uSide][u])
				C[uSide].pushBack(w);
			else
				X[uSide].pushBack(w);
		}
	}	

	// construct degC

	for (uint32_t side = 0; side <= 1; ++side) {
		for (uint32_t v : C[side]) {
			for (uint32_t w : G.nbr[side][v])
				++degC[side^1][w];
		}
	}


	// Prune C & X with lemmas
	if (!pruneCX(u, uSide)) return false;

	// Construct 3-hop
	for (uint32_t v : C[uSide])
		for (uint32_t w : G.nbr[uSide][v]) {
			if (C[uSide^1].inside(w)) continue;
			if (degC[uSide^1][w] + k >= lb) {
				C[uSide^1].pushBack(w);
				for (uint32_t x : G.nbr[uSide^1][w])
					++degC[uSide][x];
			}
		}

	if (!pruneCX(u, uSide)) return false;

	// Check correctness of degC
	// for (uint32_t side = 0; side <= 1; ++side) {
	// 	for (uint32_t u = 0; u < G.n[side]; ++u) {
	// 		if (!S[side].inside(u) && !C[side].inside(u) && !X[side].inside(u)) continue;
	// 		uint32_t degCu = 0;
	// 		for (uint32_t v : G.nbr[side][u])
	// 			if (C[side^1].inside(v))
	// 				++degCu;
	// 		assert(degCu == degC[side][u]);
	// 	}
	// }


	// Update non-neighbors in S of C & X
	for (uint32_t v : C[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].push_back(u);

	for (uint32_t v : X[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].push_back(u);

	return true;
}

void biplex::destructSets(uint32_t u, uint32_t uSide) {
	for (uint32_t v : C[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].pop_back();

	for (uint32_t v : X[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].pop_back();

	for (uint32_t side = 0; side <= 1; ++side) {
		for (uint32_t v : C[side]) {
			for (uint32_t w : G.nbr[side][v])
				--degC[side^1][w];
		}
		S[side].clear();
		C[side].clear();
		X[side].clear();
	}
}

/**
 * Enumeration algorithm
 */
biplex::Result biplex::enumeration(const BiGraph& Graph, uint32_t q, uint32_t k, uint64_t resultNum) 
{
	biplex::q = q;
	biplex::k = k;
	biplex::lb = q;// std::max(biplexLowerBound(Graph), q);
	biplex::resultNumThres = resultNum;

	pivotingTime = bipartiteTime = calcTime = 0;

	printf("\nMacro Info:\n");
	PRINT_MACRO_INFO(DEBUG);
	PRINT_MACRO_INFO(CALC_DEGREE);
	PRINT_MACRO_INFO(WITH_PIVOTING);
	PRINT_MACRO_INFO(ENUM_LARGER_DEGREE_SIDE);
	PRINT_MACRO_INFO(UPPERBOUND_PRUNE);
	PRINT_MACRO_INFO(DESTRUCT_SETS);
	PRINT_MACRO_INFO(NO_ORDERING);
	PRINT_MACRO_INFO(PRUNE_C);
	PRINT_MACRO_INFO(CHOOSE_FIRST_VERTEX);


	// Get the (lb-k, lb-k)-core from G
	if (lb > k)	
		G = coreReduction(Graph, lb-k, lb-k);
	else 
		G = std::move(Graph);

	// Resize sets
	std::vector<std::vector<uint32_t>>(k+1, std::vector<uint32_t>(G.n[1])).swap(B);
	std::vector<uint32_t>(std::max(G.n[0], G.n[1])).swap(sup);	

	for (uint32_t side = 0; side <= 1; ++side) {
		S[side].reserve(G.n[side]);
		C[side].reserve(G.n[side]);
		X[side].reserve(G.n[side]);
		Q[side].reserve(G.n[side]);
		std::vector<std::vector<uint32_t>>(G.n[side]).swap(nonNbrS[side]);
		std::vector<uint32_t>(G.n[side]).swap(degC[side]);	
	}

	result.numBranches = result.numBiPlexes = 0;

	// Butterfly-based reduction
	if (lb > 2 * k) {
		printf("\nStart butterfly reduction\n");
		printf("Before: m = %d, nL = %d, nR = %d\n", G.m, G.n[0], G.n[1]);
		G = butterflyReduction(G, (lb-k-1)*(lb-2*k-1));
		printf("After: m = %d, nL = %d, nR = %d\n", G.m, G.n[0], G.n[1]);
	}

#ifdef NO_ORDERING
	if (true) {
#else
	if (lb <= 2 * k) {
#endif
		for (uint32_t side = 0; side <= 1; ++side) {

			std::vector<std::vector<uint32_t>>(G.n[side]).swap(nonNbrS[side]);
			std::vector<uint32_t>(G.degree[side]).swap(degC[side]);

			for (uint32_t u = 0; u < G.n[side]; ++u) {
				nonNbrS[side][u].reserve(k);
				C[side].pushBack(u);
			}
		}
		printf("\nStart enumeration\n");
		startTime = std::chrono::steady_clock::now();
		branchNew(0);
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime);
		result.time = duration.count();
	
	} else {

		// Core ordering
		printf("\nStart core ordering\n");
		coreOrdering(G);

		// Choose a proper side to enumerate
#ifdef ENUM_LARGER_DEGREE_SIDE
		uint32_t uSide = G.maxDeg[0] > G.maxDeg[1] ? 0 : 1;
#else
		uint32_t uSide = G.maxDeg[0] < G.maxDeg[1] ? 0 : 1;
#endif

		printf("\nStart enumeration\n");
		startTime = std::chrono::steady_clock::now();
		// Enumerate vertices in a decreasing core order
		for (uint32_t i = 0; i < G.n[uSide]; ++i) {
			uint32_t u = ordered[uSide][i];
			if (constructSets(u, uSide)) {
#ifdef DEBUG
				printf("-------------------- Enumeration starts from %d --------------------\n", G.vLabel[uSide][u]);
#endif

#ifdef SMALL_K
				if (k <= SMALL_K) {
					nonNbrC.clear();
					for (uint32_t i = C[uSide^1].frontPos(); i < C[uSide^1].backPos(); ++i) {
						uint32_t v = C[uSide^1][i];
						if (!G.connect(uSide, u, v)) {
							nonNbrC.push_back(v);
							C[uSide^1].popBack(v); --i;
							for (uint32_t w : G.nbr[uSide^1][v])
								--degC[uSide][w];
							X[uSide^1].pushBack(v);
						}
					}
					branchSmallK(0, 0, u, uSide);
				} 
				else
#endif
					branchNew(0);
			}

#ifdef DESTRUCT_SETS
			destructSets(u, uSide);
#endif
		}
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime);
		result.time = duration.count();

	}

	

	return result;
}

uint32_t biplex::biplexUpperBound(uint32_t u, uint32_t uSide) 
{
	std::vector<uint32_t>(k+1).swap(ptrB);
	uint32_t supSum = 0;

	for (uint32_t v : S[uSide]) {
		// if (v == u) continue;
		sup[v] = k - nonNbrS[uSide][v].size();
		supSum += sup[v];
	}

	if (G.nbr[uSide][u].size() < C[uSide^1].size()) {
		for (uint32_t w : G.nbr[uSide][u]) {
			if (!C[uSide^1].inside(w)) continue;
			uint32_t s = nonNbrS[uSide^1][w].size();
			B[s][ptrB[s]++] = w;
		}
	}
	else {
		for (uint32_t w : C[uSide^1]) {
			if (!G.connect(uSide, u, w)) continue;
			uint32_t s = nonNbrS[uSide^1][w].size();
			B[s][ptrB[s]++] = w;
		}
	}

	uint32_t ub = S[uSide^1].size() + k - nonNbrS[uSide][u].size();

	for (uint32_t i = 0; i <= k && i <= supSum; ++i) {
		for (uint32_t j = 0; j < ptrB[i] && i <= supSum; ++j) {

			uint32_t v = B[i][j], x = -1;

			for (uint32_t w : nonNbrS[uSide^1][v])
				if (x == -1u || sup[w] < sup[x])
					x = w;

			if (x == -1u || sup[x] > 0) {
				if (x != -1u) --sup[x];
				supSum -= i;
				++ub;
			}
		}
	}

	return ub;
}


void biplex::printSet(VertexSet* V, const std::string& name) {
	std::vector<uint32_t> sortedV[2];
	for (uint32_t i = 0; i <= 1; ++i) {
		sortedV[i].reserve(V[i].size());
		for (uint32_t v : V[i])
			sortedV[i].push_back(G.vLabel[i][v]);
		std::sort(sortedV[i].begin(), sortedV[i].end());
	}
	printf("%sL: ", name.c_str());
	for (uint32_t v : sortedV[0]) printf("%d, ", v);
	printf("\n%sR: ", name.c_str());
	for (uint32_t v : sortedV[1]) printf("%d, ", v);
	printf("\n");
}

void biplex::printSets() {
	printSet(S, "S");
	printSet(C, "C");
	printSet(X, "X");

	VertexSet SC[2] = {VertexSet(G.n[0]), VertexSet(G.n[1])};
	for (uint32_t i = 0; i <= 1; ++i) {
		for (uint32_t v : S[i])
			SC[i].pushBack(v);
		for (uint32_t v : C[i])
			SC[i].pushBack(v);
	}

	printSet(SC, "S&C");

}

void biplex::printResult(Result& result) {
	printf("Number of branches: %" PRIu64 ", bipartite: %" PRIu64 ", pivoting: %" PRIu64 "\n",
			result.numBranches,
			result.numBipartiteBranches,
			result.numPivotingBranches);
	printf("Number of biplexes: %" PRIu64 "\n", result.numBiPlexes);
	printf("Branch time: %" PRIu64 " ms\n", result.time);
/*
	std::cout << "Number of branches: " << result.numBranches << \
		", bipartite: " << result.numBipartiteBranches << \
		", pivoting: " << result.numPivotingBranches << std::endl;
	std::cout << "Number of biplexes: " << result.numBiPlexes << std::endl;
	std::cout << "Branch time: " << result.time << "ms" << std::endl;
*/
}

void biplex::branchSmallK(uint32_t dep, uint32_t begin, uint32_t u, uint32_t uSide)
{

	branchNew(0);

	if ((uint32_t)dep == k) return;

	std::pair<uint32_t, uint32_t> oldPosC, oldPosX;

	for (uint32_t i = begin; i < nonNbrC.size(); ++i) {
		uint32_t& v = nonNbrC[i];

		biplexUpdate(v, uSide^1, oldPosC, oldPosX);
		biplexUpdateDeg(v, uSide^1, oldPosC);

		branchSmallK(dep+1, i+1, u, uSide);

		biplexRestoreDeg(v, uSide^1, oldPosC);
		biplexRestore(v, uSide^1, oldPosC, oldPosX);

	}

}

uint32_t biplex::calcDegree(VertexSet* V, uint32_t u, uint32_t uSide)
{
	uint deg = 0;
	if (G.nbr[uSide][u].size() < V[uSide^1].size()) {
		for (uint32_t v : G.nbr[uSide][u])
			if (V[uSide^1].inside(v))
				++deg;
	}
	else {
		for (uint32_t v : V[uSide^1])
			if (G.connect(uSide, u, v))
				++deg;
	}
	return deg;	
}



void biplex::biplexPruneC(std::pair<uint32_t, uint32_t>& oldPosC) {

	oldPosC = std::make_pair(C[0].frontPos(), C[1].frontPos());

#ifdef PRUNE_ONCE	

	for (uint32_t side = 0; side <= 1; ++side) {
		for (uint32_t u : C[side]) {
#ifdef CALC_DEGREE
			degC[side][u] = calcDegree(C, u, side);
#endif
			int32_t degCu = degC[side][u];
			if (degCu + (S[side^1].size() - nonNbrS[side][u].size()) + k < lb) {
				C[side].popFront(u);
				for (uint32_t v : G.nbr[side][u])
					--degC[side^1][v];
			}
		}
	}

#else

#ifdef PRUNE_WITH_QUEUE

	for (uint32_t side = 0; side <= 1; ++side) {
		Q[side].clear();
		for (uint32_t u : C[side]) {
#ifdef CALC_DEGREE
			degC[side][u] = calcDegree(C, u, side);
#endif
			int32_t degCu = degC[side][u];
			if (degCu + (S[side^1].size() - nonNbrS[side][u].size()) + k < lb) {
				Q[side].pushBack(u);
			}
		}
	}

	while (Q[0].size() != 0 || Q[1].size() != 0) {
		for (uint32_t side = 0; side <= 1; ++side) {
			while (Q[side].size() != 0) {
				uint32_t u = Q[side][0];
				Q[side].popBack(u);
				C[side].popFront(u);
				for (uint32_t v : G.nbr[side][u]) {
					--degC[side^1][v];
					if (C[side^1].inside(v) && degC[side^1][v] + (S[side].size() - nonNbrS[side^1][v].size()) + k < lb)
				 		Q[side^1].pushBack(v);
				}
			}
		}
	}

#else
	while (true) {
		bool flag = false;
		for (uint32_t side = 0; side <= 1; ++side) {
			for (uint32_t u : C[side]) {
#ifdef CALC_DEGREE
				degC[side][u] = calcDegree(C, u, side);
#endif
				int32_t degCu = degC[side][u];
				if (degCu + (S[side^1].size() - nonNbrS[side][u].size()) + k < lb) {
					flag = true;
					C[side].popFront(u);
					for (uint32_t v : G.nbr[side][u])
						--degC[side^1][v];
				}
			}
		}
		if (!flag) break;
	}
#endif // PRUNE_WITH_QUEUE

#endif // PRUNE_ONCE

}

void biplex::biplexRestorePruneC(std::pair<uint32_t, uint32_t>& oldPosC) {
	biplexRestoreDegC(oldPosC);
	C[0].restore(oldPosC.first);
	C[1].restore(oldPosC.second);
}

/**
 * The main algorithm.
 * */
void biplex::branchNew(uint32_t dep) 
{

	++result.numBranches;

	if (result.numBiPlexes >= resultNumThres) return;

#ifdef DEBUG
	printf("\n---------- dep = %d ----------\n", dep);
	printSets();
#endif

	std::pair<uint32_t, uint32_t> oldPosC, oldPosX, pruneOldPosC;
	
	// ********************************** Prune C **********************************
#ifdef PRUNE_C
	biplexPruneC(pruneOldPosC);	
#endif


	// ************************ Prune branch & update result ***********************
	if ((S[0].size() + C[0].size() < lb) || (S[1].size() + C[1].size() < lb)) {
#ifdef PRUNE_C
		biplexRestorePruneC(pruneOldPosC);
#endif
		return;
	}

	if (C[0].size() == 0 && C[1].size() == 0) {
		if (X[0].size() == 0 && X[1].size() == 0) {
			if (S[0].size() >= lb && S[1].size() >= lb) {
				++result.numBiPlexes;
				if (result.numBiPlexes == resultNumThres) {
					printf("Done!\n");
					auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - startTime);
					result.time = duration.count();
					printResult(result);
					exit(0);
				}
#ifdef DEBUG
				printf("*** find: No.%d\n", result.numBiPlexes);
#endif
			}
		}
		// ++result.numBranches;
#ifdef PRUNE_C
		biplexRestorePruneC(pruneOldPosC);
#endif
		return;
	}


	// ************************* Calculate minimum degree *************************

	int32_t inf = G.n[0] + G.n[1];
	uint32_t Sp[2] = {-1u, -1u};
	int32_t nonDegPInS[2] = {-inf, -inf};
	uint32_t Cp[2] = {-1u, -1u};
	int32_t nonDegPInC[2] = {-inf, -inf};


	for (uint32_t side = 0; side <= 1; ++side) {
		for (uint32_t u : S[side]) {
#ifdef CALC_DEGREE
			degC[side][u] = calcDegree(C, u, side);
#endif
			int32_t nonDegU = C[side^1].size() - degC[side][u] + nonNbrS[side][u].size();
			if (nonDegU > nonDegPInS[side]) {
				nonDegPInS[side] = nonDegU;
				Sp[side] = u;
			}
		}

		// minimum degree of vertices in S < q - k
		if ((int32_t)(S[side^1].size() + C[side^1].size()) - nonDegPInS[side] < (int32_t)(lb - k))  {
#ifdef PRUNE_C
			biplexRestorePruneC(pruneOldPosC);
#endif
			return;
		}
	}

#ifndef WITH_PIVOTING

	for (uint32_t side = 0; side <= 1; ++side) {
		for (uint32_t u : C[side]) {
#ifdef CALC_DEGREE
			degC[side][u] = calcDegree(C, u, side);
#endif
			int32_t nonDegU = C[side^1].size() - degC[side][u] + nonNbrS[side][u].size();
			if (nonDegU > nonDegPInC[side]) {
				nonDegPInC[side] = nonDegU;
				Cp[side] = u;
			}
		}
	}

	// maximum non-degree of vertices in G[S & C] < k
	if (std::max(nonDegPInC[0], nonDegPInS[0]) <= (int32_t)k && 
		std::max(nonDegPInC[1], nonDegPInS[1]) <= (int32_t)k &&
		X[0].size() == 0 && X[1].size() == 0) {
		++result.numBiPlexes;
#ifdef DEBUG
		printf("*** find: No.%d\n", result.numBiPlexes);
#endif
#ifdef PRUNE_C
		biplexRestorePruneC(pruneOldPosC);
#endif
		return;
	}


	
	

#endif // WITH_PIVOTING


	// ********************************** Bipartite **********************************
#ifdef WITH_PIVOTING

	if (nonDegPInS[0] > (int32_t)k || nonDegPInS[1] > (int32_t)k) {
		

		uint32_t side = nonDegPInS[0] > nonDegPInS[1] ? 0 : 1;

		uint32_t u = -1;
		int32_t deg_u = G.n[side], nonDegSu = -1;


		for (uint32_t v : C[side^1]) {
#ifdef CALC_DEGREE
			// degC[side^1][v] = calcDegree(C, v, side^1);
#endif
			int32_t nonDegSv = nonNbrS[side^1][v].size();
			int32_t deg_v = degC[side^1][v] + S[side].size() - nonDegSv;
			if (((deg_v < deg_u) || (deg_v == deg_u && nonDegSv > nonDegSu)) && !G.connect(side, Sp[side], v))  {
				u = v;
				deg_u = deg_v;
				nonDegSu = nonDegSv;
#ifdef CHOOSE_FIRST_VERTEX
				break;
#endif
			}
		}



#else	// WITH_PIVOTING	
	if (true) {

		uint32_t side = S[1].size() + C[1].size() - nonDegPInC[0] < S[0].size() + C[0].size() - nonDegPInC[1] ? 1 : 0;
		uint32_t u = Cp[side^1];

#endif	// WITH_PIVOTING

#ifdef DEBUG
		printf("Bipartite: %d%c\n", G.vLabel[side^1][u], side^1 ? 'R' : 'L');
#endif


		assert(u != -1u);
		biplexUpdate(u, side^1, oldPosC, oldPosX);
		if (S[0].size() + C[0].size() >= lb && S[1].size() + C[1].size() >= lb) {
#ifdef UPPERBOUND_PRUNE
			if (biplexUpperBound(u, side^1) >= lb) 
#endif
			{
				biplexUpdateDeg(u, side^1, oldPosC);
				++result.numBipartiteBranches;				
				branchNew(dep + 1);
				biplexRestoreDeg(u, side^1, oldPosC);
			}
		}
		biplexRestore(u, side^1, oldPosC, oldPosX);
		if (S[0].size() + C[0].size() >= lb && S[1].size() + C[1].size() >= lb) {
			++result.numBipartiteBranches;				
			branchNew(dep + 1);
		}
		biplexBacktrack(u, side^1);

		
	}


		


	// ********************************** Pivoting **********************************
	else {

		// Choose a pivot
		uint32_t p[2] = {-1u, -1u};
		int32_t degCp[2] = {-1, -1};
		
		for (uint32_t side = 0; side <= 1; ++side) {
			for (uint32_t u : C[side]) {
#ifdef CALC_DEGREE
				degC[side][u] = calcDegree(C, u, side);
#endif
				int32_t degCu = degC[side][u];
				if (degCu > degCp[side]) {
					p[side] = u;
					degCp[side] = degCu;
#ifdef CHOOSE_FIRST_VERTEX
					break;
#endif
				}
			}
#ifndef CHOOSE_FIRST_VERTEX
			for (uint32_t u : X[side]) {
#ifdef CALC_DEGREE
				degC[side][u] = calcDegree(C, u, side);
#endif
				int32_t degCu = degC[side][u];
				if (degCu > degCp[side]) {

					bool flag = false;

					for (uint32_t v : nonNbrS[side][u]) {

						if (C[side].size() - degC[side^1][v] + nonNbrS[side^1][v].size() >= k) {
							flag = true;
							break;
						}
					}

					if (!flag) {
						p[side] = u;
						degCp[side] = degCu;
					}
				}
			}
#endif // CHOOSE_FIRST_VERTEX
		}

		uint32_t pSide = C[1].size() - degCp[0] < C[0].size() - degCp[1] ? 0 : 1;
		uint32_t pivot = p[pSide];
	
#ifdef DEBUG
		printf("Pivoting: %d%c\n", G.vLabel[pSide][pivot], pSide ? 'R' : 'L');
#endif

		// Generate branches
		if (cand.size() <= dep) cand.resize(dep * 2);
		std::vector<uint32_t>& Cand = cand[dep];
		
		Cand.clear();
		
		Cand.reserve(C[pSide^1].size() - degCp[pSide]);

		for (uint32_t u : C[pSide^1])
			if (!G.connect(pSide, pivot, u))
				Cand.push_back(u);


		// branch with non-neighbors of pivot
		for (uint32_t u : Cand) {

			biplexUpdate(u, pSide^1, oldPosC, oldPosX);
			if (S[0].size() + C[0].size() >= lb && S[1].size() + C[1].size() >= lb) {
#ifdef UPPERBOUND_PRUNE
				if (biplexUpperBound(u, pSide^1) >= lb)
#endif
				{
					++result.numPivotingBranches;
					biplexUpdateDeg(u, pSide^1, oldPosC);
					branchNew(dep + 1);
					biplexRestoreDeg(u, pSide^1, oldPosC);
				}
			}
			biplexRestore(u, pSide^1, oldPosC, oldPosX);
		}
		
		// Branch with pivot
		if (C[pSide].inside(pivot)) {
			biplexUpdate(pivot, pSide, oldPosC, oldPosX);
			if (S[0].size() + C[0].size() >= lb && S[1].size() + C[1].size() >= lb) {
#ifdef UPPERBOUND_PRUNE
				if (biplexUpperBound(pivot, pSide) >= lb)
#endif
				{
					++result.numPivotingBranches;
					biplexUpdateDeg(pivot, pSide, oldPosC);
					branchNew(dep + 1);
					biplexRestoreDeg(pivot, pSide, oldPosC);
				}
			}
			biplexRestore(pivot, pSide, oldPosC, oldPosX);
			biplexBacktrack(pivot, pSide);
		}

		// Backtrack from Cand
		for (uint32_t u : Cand) {
			biplexBacktrack(u, pSide^1);
		}	
	}


#ifdef PRUNE_C
	biplexRestorePruneC(pruneOldPosC);
#endif

#ifdef DEBUG
	printf("\n>>>>>>>>>> dep = %d <<<<<<<<<<\n", dep);
	printSets();	
#endif

}

void biplex::biplexUpdateDegC(std::pair<uint32_t, uint32_t>& oldPosC) {
#ifndef CALC_DEGREE
	// Update degC
	for (uint32_t side = 0; side <= 1; ++side) {
		uint32_t oldPos = side == 0 ? oldPosC.first : oldPosC.second;
		for (uint32_t i = oldPos; i < C[side].frontPos(); ++i) {
			uint32_t v = C[side][i];
			for (uint32_t w : G.nbr[side][v])
				--degC[side^1][w];
		}
	}
#endif
}



void biplex::biplexUpdate(uint32_t u, uint32_t uSide, 
	std::pair<uint32_t, uint32_t>& oldPosC, 
	std::pair<uint32_t, uint32_t>& oldPosX) 
{
	// Update C & X
	biplexUpdateCX(u, uSide, C, oldPosC);
	biplexUpdateCX(u, uSide, X, oldPosX);

	// Update S
	S[uSide].pushBack(u);
}

void biplex::biplexUpdateDeg(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC) 
{

	// Update nonNbrS of S, C and X
	for (uint32_t v : S[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].push_back(u);

	for (uint32_t v : C[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].push_back(u);

	for (uint32_t v : X[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].push_back(u);

	biplexUpdateDegC(oldPosC);

}


void biplex::biplexUpdateCX(uint32_t u, uint32_t uSide, VertexSet* V, std::pair<uint32_t, uint32_t>& oldPos) 
{
	oldPos = std::make_pair(V[0].frontPos(), V[1].frontPos());
	V[uSide].popFront(u);

	bool flagNonNbrSu = nonNbrS[uSide][u].size() >= k;
	for (uint32_t v : V[uSide^1]) {
		if (G.connect(uSide, u, v)) continue;
		if (flagNonNbrSu || nonNbrS[uSide^1][v].size() >= k) {
			V[uSide^1].popFront(v);
		}
	}

	for (uint32_t v : nonNbrS[uSide][u]) {
		if (nonNbrS[uSide^1][v].size() == k-1) {
			if (V[uSide].size() < G.nbr[uSide^1][v].size()) {
				for (uint32_t w : V[uSide]) {
					if (!G.connect(uSide^1, v, w)) {
						V[uSide].popFront(w);
					}
				}
			}
			else {
				uint32_t newFrontPos = V[uSide].backPos();
				for (uint32_t w : G.nbr[uSide^1][v])
					if (V[uSide].inside(w))
						V[uSide].swapByVal(w, V[uSide][--newFrontPos]);
				V[uSide].restore(newFrontPos);
			}
		}
	}
}

void biplex::biplexRestoreDegC(std::pair<uint32_t, uint32_t>& oldPosC) {
#ifndef CALC_DEGREE
	// Restore degC
	for (uint32_t side = 0; side <= 1; ++side) {
		uint32_t oldPos = side == 0 ? oldPosC.first : oldPosC.second;
		for (uint32_t i = oldPos; i < C[side].frontPos(); ++i) {
			uint32_t v = C[side][i];
			for (uint32_t w : G.nbr[side][v])
				++degC[side^1][w];

		}
	}

#endif
}

void biplex::biplexRestoreDeg(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC) 
{

	// Restore nonNbrS
	for (uint32_t v : S[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].pop_back();

	for (uint32_t v : C[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].pop_back();

	for (uint32_t v : X[uSide^1])
		if (!G.connect(uSide, u, v))
			nonNbrS[uSide^1][v].pop_back();


	biplexRestoreDegC(oldPosC);

}

void biplex::biplexRestore(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC, std::pair<uint32_t, uint32_t>& oldPosX) 
{	

	// Restore C
	C[0].restore(oldPosC.first);
	C[1].restore(oldPosC.second);

	// Restore X
	X[0].restore(oldPosX.first);
	X[1].restore(oldPosX.second);

#ifndef CALC_DEGREE
	if (C[uSide].inside(u)) {
		for (uint32_t v : G.nbr[uSide][u])
			--degC[uSide^1][v];
	}
#endif

	// Move u from C to X

	S[uSide].popBack(u);
	C[uSide].popFront(u);
	X[uSide].pushBack(u);

}


void biplex::biplexBacktrack(uint32_t u, uint32_t uSide)
{

#ifndef CALC_DEGREE
	if (!C[uSide].inside(u)) {
		for (uint32_t v : G.nbr[uSide][u])
			++degC[uSide^1][v];
	}
#endif

	C[uSide].pushFront(u);
	X[uSide].popBack(u);
}
