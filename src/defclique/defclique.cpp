#include "defclique.h"
#include "../utils/ordering.hpp"
#include "../utils/coloring.hpp"
#include "../utils/log.hpp"
#include <chrono>
#include <cstring>
#include <ratio>
#include <vector>
#include <cstdio>
#include <sstream>

//#define EDGE_REDUCTION
#define DEBUG_RESULT
// #define DEBUG_BRANCH
#define UPPERBOUND_FULL


namespace defclique {
	int k, nnbS;
    // long long nnbSub;
	Graph Sub;
	Coloring clr;
	VertexSet S, C, Ss, C1, D;
	std::vector<int> degC1, degC, degS, cnD;
	std::vector<int> q;
	int mode;
}

void defclique::logSet(VertexSet &V, const std::string &name) {
	std::vector<int> S(V.begin(), V.end());
	std::sort(S.begin(), S.end());
	std::stringstream ss;
	ss << "{";
	for (int i = 0; i < S.size(); ++i) {
		if (i > 0) ss << ",";
		ss << S[i];
	}
	ss << "}";
	log("%s: size=%d, content=%s", 
		name.c_str(), V.size(), ss.str().c_str());
}

Graph defclique::coreReduction(Graph& G, int k) {
	if (k <= 1) return G;

	log("Running core reduction with k=%d...", k);

	auto startTimePoint = std::chrono::steady_clock::now();

	int head = 0, tail = 0;
	std::vector<bool> vis(G.n);
	std::vector<int> deg(G.n);

	for (int u : G.V) {
		deg[u] = G.nbr[u].size();
		if (deg[u] < k) {
			vis[u] = true;
			q[tail++] = u;
		}
	}

	while (head < tail) {
		int u = q[head++];
		for (int v : G.nbr[u]) {
			if (vis[v]) continue;
			if (--deg[v] < k) {
				vis[v] = true;
				q[tail++] = v;
			}
		}
	}

	Graph C;

	for (int u : G.V) {
		if (vis[u]) continue;
		for (int v : G.nbr[u]) {
			if (vis[v]) continue;
			C.addEdge(u, v);
		}
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Core reduction done! Time spent: %ld ms", duration.count());
	log("Before: n=%d, m=%d; After: n=%d, m=%d", G.V.size(), G.m, C.V.size(), C.m);

	return C;
}

Graph defclique::edgeReduction(Graph &G, int k) {

	log("Running edge reduction with k=%d...", k);

	auto startTimePoint = std::chrono::steady_clock::now();

	std::vector<int> cn(G.m), q(G.m);
	std::vector<std::pair<int, int>> edges;
	std::vector<std::unordered_map<int, int>> eid(G.n);	
	std::vector<bool> vis(G.m);

	edges.reserve(G.m);

	int head = 0, tail = 0;

	auto countCommonNeighbor = [&](int u, int v) {
		int cnt = 0;
		if (G.nbr[u] > G.nbr[v]) std::swap(u, v);
		for (int w : G.nbr[u])
			if (G.connect(v, w))
				++cnt;
		return cnt;
	};

	auto removeEdge = [&](int id) {
		vis[id] = true;
		// cn[id] = 0;
		int u = edges[id].first, v = edges[id].second;
		if (G.nbr[u] > G.nbr[v]) std::swap(u, v);
		for (int w : G.nbr[u])
			if (G.connect(v, w)) {
				int i = eid[u][w], j = eid[v][w];
				if (vis[i] || vis[j]) continue;
				--cn[i], --cn[j];
			}		
	};

	for (int u : G.V)
		for (int v : G.nbr[u])
			if (u < v) {
				int i = edges.size();
				eid[u][v] = eid[v][u] = i;
				edges.push_back(std::make_pair(u, v));
				cn[i] = countCommonNeighbor(u, v);
				if (cn[i] < k) {
					q[tail++] = i;
					removeEdge(i);
				}
			}

	while (head < tail) {
		int i = q[head++];
		int u = edges[i].first, v = edges[i].second;
		if (G.nbr[u] > G.nbr[v]) std::swap(u, v);
		for (int w : G.nbr[u]) {
			if (G.connect(v, w)) {
				int i = eid[u][w], j = eid[v][w];
				if (!vis[i] && cn[i] < k) {
					q[tail++] = i;
					removeEdge(i);
				}
				if (!vis[j] && cn[j] < k) {
					q[tail++] = j;
					removeEdge(j);
				}
			}
		}
	}

	Graph E;

	for (int i = 0; i < edges.size(); ++i) {
		if (vis[i]) continue;
		E.addEdge(edges[i].first, edges[i].second);
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Edge reduction done! Time spent: %ld ms", duration.count());
	log("Before: n=%d, m=%d; After: n=%d, m=%d", G.V.size(), G.m, E.V.size(), E.m);

	return E;

}

void defclique::preprocessing(Graph &G, Ordering &o, int u, int mode) {
	S.clear();
	C.clear();
	C1.clear();
	nnbS = 0; //nnbSub = 0;
	S.push(u);

	int i = o.order[u];

	// Construct C1
	if (G.nbr[u].size() < o.numOrdered-i-1) {
		for (int v : G.nbr[u])
			if (o.order[v] > i)
				C1.push(v);
	}
	else {
		for (int j = i+1; j < o.numOrdered; ++j) {
			int v = o.ordered[j];
			if (G.connect(u, v))
				C1.push(v);
		}
	}
	if (Ss.size() >= C1.size()+1+k) {
		S.clear();
		C.clear();
		C1.clear();
		return;
	}

	Sub.subGraph(G, C1);

	for (int v : C1)
		degC1[v] = Sub.nbr[v].size();

	int head = 0, tail = 0;

	for (int v : C1) {
		if (degC1[v] < Ss.size()-k-1) {
			q[tail++] = v;
			sub(Sub, C1, degC1, v);
		}
	}

	while (head < tail) {
		int v = q[head++];
		if (Sub.nbr[v].size() < C1.size()) {
			for (int w : Sub.nbr[v]) {
				if (C1.inside(w) && degC1[w] < Ss.size()-k-1) {
					q[tail++] = w;
					sub(Sub, C1, degC1, w);
				}
			}
		}
		else {
			for (int w : C1)
				if (Sub.connect(v, w) && degC1[w] < Ss.size()-k-1) {
					q[tail++] = w;
					sub(Sub, C1, degC1, w);
				}
		}
	}

	if (mode == TWO_HOP) {
		// Construct C2
		for (int v : C1) {
			if (G.nbr[v].size() < o.numOrdered-i-1) {
				for (int w : G.nbr[v]) {
					if (o.order[w] > i && !C1.inside(w)) {
						degC1[w] = 0;
						C.push(w);
					}
				}
			}
			else {
				for (int j = i+1; j < o.numOrdered; ++j) {
					int w = o.ordered[j];
					if (G.connect(v, w) && !C1.inside(w)) {
						degC1[w] = 0;
						C.push(w);
					}
				}
			}
		}	

		for (int v : C1) {
			if (C.size() < G.nbr[v].size()) {
				for (int w : C)
					if (G.connect(v, w))
						++degC1[w];
			}
			else {
				for (int w : G.nbr[v])
					if (C.inside(w))
						++degC1[w];
			}
		}

		for (int v : C)
			if (degC1[v] < Ss.size()-k)
				C.pop(v);
	}

	Sub.subGraph(G, S, C1, C);

	degC1[u] = C1.size();

	for (int v : Sub.V)
		degS[v] = degC[v] = 0;

	for (int v : Sub.nbr[u])
		degS[v] = 1;

	for (int v : C) {
		for (int w : Sub.nbr[v])
			++degC[w];
	}

	for (int v : C1)
		C.push(v);

	for (int v : Sub.V)
		degC[v] += degC1[v];

}

VertexSet defclique::heuristic(Graph &G, int k) {

	log("Running heuristic algorithm...");

	auto startTimePoint = std::chrono::steady_clock::now();

	Ordering o, oSub;

	o.degeneracyOrdering(G);

	for (int i = o.numOrdered-1; i >= 0; --i) {
		// fprintf(stderr, "%d/%d\r", o.numOrdered-i, o.numOrdered);

		int u = o.ordered[i];
		if (o.value[u] < Ss.size()-k) break;

		preprocessing(G, o, u, ONE_HOP);

		oSub.degeneracyOrdering(Sub);

		int j = oSub.numOrdered;	

		while (C.size() > 0) {

			// Select a vertex from C with maximum degeneracy
			int u = oSub.ordered[--j];
			for (; !C.inside(u); u = oSub.ordered[--j]);


			// Add u from C to S
			nnbS += S.size() - degS[u];
			add(Sub, S, degS, u);
			sub(Sub, C, degC, u);

			// Prune C 
			int head = 0, tail = 0;
			for (int v : C) if (nnbS + S.size()-degS[v] > k || degC[v]+degS[v] < Ss.size()-k + nnbS) {
				sub(Sub, C, degC, v);
				q[tail++] = v;
			}

			while (head < tail) {
				int v = q[head++];
				if (Sub.nbr[v].size() < C.size()) {
					for (int w : Sub.nbr[v]) if (degC[w]+degS[w] < Ss.size()-k + nnbS && C.inside(w)) {
						sub(Sub, C, degC, w);
						q[tail++] = w;
					}
				}
				else {
					for (int w : C) if (degC[w]+degS[w] < Ss.size()-k + nnbS && Sub.connect(v, w)) {
						sub(Sub, C, degC, w);
						q[tail++] = w;
					}
				}
			}

			bool flagBreak = false;

			for (int v : S) if (degC[v] <= Ss.size()-k - S.size() + nnbS) {
				C.clear();
				flagBreak = true;
				break;
			}

			if (flagBreak) break;
			
		}

		if (Ss.size() < S.size()) {
			Ss.clear();
			for (int v : S) Ss.push(v);
		}

		S.clear(); 
		S.push(u);
		for (int v : Sub.V) 
			degS[v] = (int)Sub.connect(u, v);

		while (C1.size() > 0) {

			// Select a vertex from C with maximum degree
			int u = C1[C1.frontPos()];
			for (int v : C1) {
				if (degS[v] + degC1[v] > degS[u] + degC1[u]) {
					u = v;
				}
			}


			// Add u from C to S
			nnbS += S.size() - degS[u];
			add(Sub, S, degS, u);
			sub(Sub, C1, degC1, u);

			// Prune C 
			int head = 0, tail = 0;
			for (int v : C1) if (nnbS + S.size()-degS[v] > k || degC1[v]+degS[v] < Ss.size()-k + nnbS) {
				sub(Sub, C1, degC1, v);
				q[tail++] = v;
			}

			while (head < tail) {
				int v = q[head++];
				if (Sub.nbr[v].size() < C1.size()) {
					for (int w : Sub.nbr[v]) if (degC1[w]+degS[w] < Ss.size()-k + nnbS && C1.inside(w)) {
						sub(Sub, C1, degC1, w);
						q[tail++] = w;
					}
				}
				else {
					for (int w : C1) if (degC1[w]+degS[w] < Ss.size()-k + nnbS && Sub.connect(v, w)) {
						sub(Sub, C1, degC1, w);
						q[tail++] = w;
					}
				}
			}

			bool flagBreak = false;

			for (int v : S) if (degC1[v] <= Ss.size()-k - S.size() + nnbS) {
				C1.clear();
				flagBreak = true;
				break;
			}

			if (flagBreak) break;
			
		}


		if (Ss.size() < S.size()) {
			Ss.clear();
			for (int v : S) Ss.push(v);
		}
	}


	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Heuristic algorithm done! Time spent: %ld ms", duration.count());
	logSet(Ss, "S*");

	return Ss;
	
}


void defclique::run(const std::string &filename, int k, int mode) {

	log("Reading graph: %s ...", strrchr(filename.c_str(), '/')+1);

	auto startTimePoint = std::chrono::steady_clock::now();

	Graph G(filename);

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);

	log("Reading graph done! Time spent: %ld ms", duration.count());
	log("Graph info: n=%d, m=%d, maxdeg=%d", G.V.size(), G.m, G.maxDeg);


	S.reserve(G.n);
	C.reserve(G.n);
	C1.reserve(G.n);
	D.reserve(G.n);
	Ss.reserve(G.n);

	degS.resize(G.n);
	degC.resize(G.n);
	degC1.resize(G.n);
	cnD.resize(G.n);
	q.resize(G.n);

	Sub.nbrMap = G.nbrMap;

	defclique::k = k;
	defclique::mode = mode;

	Ss = heuristic(G, k);

	if (Ss.size() < k+1) {
		for (int v : G.V) {
			Ss.push(v);
			if (Ss.size() == k+1) break;
		}
	}
	

	Graph Core = coreReduction(G, Ss.size() - k);
#ifdef EDGE_REDUCTION
	Core = edgeReduction(Core, Ss.size() - k - 1);
#endif
	Ordering o = Ordering::DegeneracyOrdering(Core);


	std::string modeString = mode == REDUCTION_SEARCH ? "Reduction" : "Russian Doll";
	log("Running %s search ...", modeString.c_str());

	startTimePoint = std::chrono::steady_clock::now();

	int branchTimeCount = 0;

	for (int i = mode == REDUCTION_SEARCH ? 0 : o.numOrdered - 1; 
		i >= 0 && i < o.numOrdered; mode == REDUCTION_SEARCH ? ++i : --i) {

		int u = o.ordered[i];

#ifdef DEBUG_BRANCH
		log("********** New branch: u=%d **********\n", u);
#endif

		if (mode == RUSSIANDOLL_SEARCH && o.value[u] < Ss.size()-k) break;
		if (mode == REDUCTION_SEARCH && o.numOrdered-i <= Ss.size()) break;

		// S.clear();
		// C.clear();

		// S.push(u);

		// if (Core.nbr[u].size() < o.numOrdered-i-1) {
		// 	for (int v : Core.nbr[u])
		// 		if (o.order[v] > i)
		// 			C.push(v);
		// }
		// else {
		// 	for (int j = i+1; j < o.numOrdered; ++j) {
		// 		int v = o.ordered[j];
		// 		if (Core.connect(u, v))
		// 			C.push(v);
		// 	}
		// }

		// if (upperbound() <= Ss.size()-k)
		// 	continue;

		// if (Ss.size() < k+1) {
		// 	S.clear();
		// 	C.clear();
		// 	S.push(u);
		// 	for (int j = i+1; j < o.numOrdered; ++j) C.push(o.ordered[j]);
		// 	Sub.subGraph(Core, S, C);
		// 	nnbS = 0;
		// 	// nnbSub = (((long long)Sub.V.size() * (Sub.V.size()-1)) >> 1) - Sub.m;
		// 	for (int v : Sub.V) degS[v] = degC[v] = 0;
		// 	for (int v : Sub.nbr[u]) degS[v] = 1;
		// 	for (int v : C) {
		// 		for (int w : Sub.nbr[v])
		// 			++degC[w];
		// 	}
			
		// }
		// else
		preprocessing(Core, o, u, TWO_HOP);
		clr.graphColoring(Sub, Ss.size()-k+1);
		auto branchStartTimePoint = std::chrono::steady_clock::now();
		branch(0);
		branchTimeCount += std::chrono::duration_cast<std::chrono::milliseconds>(
			std::chrono::steady_clock::now() - branchStartTimePoint).count();
	}
	

	auto totalTimeCount = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint).count();


	log("%s search done! Branch time: %ld ms, total time: %ld ms", 
		modeString.c_str(), branchTimeCount, totalTimeCount);

	if (Ss.size() < k+2) {
		Ss.clear();
		log("Warning: unable to find a defective clique with size larger than k+2.")		
	}

	logSet(Ss, "S*");

#ifdef DEBUG_RESULT

	Sub.subGraph(G, Ss);

	int cnt = 0;
	for (int v : Sub.V) {
		for (int w : Sub.V)
			if (v < w && !Sub.connect(v, w)) {
				log("Missing edge in S*: (%d, %d)", v, w);
				++cnt;
			}
	}
	log("Number of missing edges in S*: %d", cnt);
#endif

}

#ifdef UPPERBOUND_FULL
int defclique::upperbound() {
	D.clear(); C1.clear();
	int s = nnbS;
	int posC = C.frontPos();
	for (int v : C) {
		cnD[clr.color[v]] = 0;
		if (S.size() == degS[v]) C1.push(v);
	}
	for (int v : C1) {
		if (cnD[clr.color[v]] == 0) {
			D.push(v);
			C.pop(v);
			cnD[clr.color[v]] = 1;
		}
	}

	while (C.size() > 0 && nnbS < k) {
		int u = -1, minv = 0x7fffffff;
		for (int v : C) {
			if (cnD[clr.color[v]] + S.size() - degS[v] < minv) {
				minv = cnD[clr.color[v]] + S.size() - degS[v];
				u = v;
			}
		}

		if ((s += minv) > k) break;

		D.push(u);
		C.pop(u);
		++cnD[clr.color[u]];
		// for (int v : C) {
		// 	if (clr.color[v] == clr.color[u]) 
		// 		++cnD[v];
		// }
	}
	C.restore(posC);
	return S.size() + D.size();
}

#elif defined(UPPERBOUND_COLOR) // Color upperbound
int defclique::upperbound() {
	int cntClr = 0;
	for (int v : S) cnD[clr.color[v]] = 0;
	for (int v : C) cnD[clr.color[v]] = 0;
	for (int v : S) if (!cnD[clr.color[v]]) { cnD[clr.color[v]] = 1; ++cntClr; }
	for (int v : C) if (!cnD[clr.color[v]]) { cnD[clr.color[v]] = 1; ++cntClr; }
	return cntClr + k;
}

#elif defined(UPPERBOUND_CORE) // Core upperbound
int defclique::upperbound() {
	C1.clear();
	for (int v : S) { C1.push(v); degC1[v] = degS[v] + degC[v]; }
	for (int v : C) { C1.push(v); degC1[v] = degS[v] + degC[v]; }
	int head = 0, tail = 0;
	for (int v : C1)
		if (degC1[v] < Ss.size()-k) {
			C1.pop(v);
			q[tail++] = v;
		}

	while (head < tail) {
		int u = q[head++];
		if (Sub.nbr[u].size() < C1.size()) {
			for (int v : Sub.nbr[u])
				if (C1.inside(v) && --degC1[v] < Ss.size()-k) {
					C1.pop(v);
					q[tail++] = v;
				}
		}
		else {
			for (int v : C1)
				if (Sub.connect(u, v) && --degC1[v] < Ss.size()-k) {
					C1.pop(v);
					q[tail++] = v;
				}
		}
	}

	if (C1.size() == 0) return Ss.size();
	return Ss.size()+1;

}
#else // No upperbound
int defclique::upperbound() {
	return Ss.size()+1;
}
#endif


void defclique::moveCToS(int v) {
	nnbS += S.size() - degS[v];
	add(Sub, S, degS, v);
	sub(Sub, C, degC, v);
}

void defclique::moveSToC(int v) {
	add(Sub, C, degC, v);
	sub(Sub, S, degS, v);
	nnbS -= S.size() - degS[v];
}

int defclique::updateC(int v) {
	int posC = C.frontPos();
	int sizeS = S.size() - (int)S.inside(v);
	for (int u : C) {
		if (u != v && nnbS + 2 * sizeS-degS[u]-degS[v] + (int)!Sub.connect(u, v) > k) {
			sub(Sub, C, degC, u);
			// nnbSub -= S.size()-degS[u] + C.size()-degC[u];
		}
	}
	return posC;
}

void defclique::restoreC(int pos) {
	for (int i = C.frontPos()-1; i >= pos; --i) {
		int u = C[i];
		// nnbSub += S.size()-degS[u] + C.size()-degC[u];
		add(Sub, C, degC, u);
	}
}

int defclique::update(int v) {
	int posC = updateC(v);
	moveCToS(v);
	return posC;
}

void defclique::restore(int v, int posC) {
	moveSToC(v);
	restoreC(posC);
}

bool defclique::branch(int dep) {

#ifdef DEBUG_BRANCH

	int cntNnbSub = 0, cntNnbS = 0;
	VertexSet V = S + C;
	for (int u : V) {
		for (int v : V) {
			if (u < v && !Sub.connect(u, v)) {
				++cntNnbSub;
				if (S.inside(u) && S.inside(v))
					++cntNnbS;
			}
		}
	}

	log("\n*** dep=%d, |S|=%d, |C|=%d, nnbS=%d(real=%d)", 
		dep, S.size(), C.size(), nnbS, cntNnbS);
	logSet(S, "S");
	logSet(C, "C");

#endif

	if (C.size() == 0) {
		if (S.size() > Ss.size()) {
#ifdef DEBUG_BRANCH
			log("*** New S*: size=%d", S.size());
#endif
			Ss.clear();
			for (int v : S) Ss.push(v);
				return mode == RUSSIANDOLL_SEARCH;
		}
		return false;
	}
	// if (nnbSub <= k) {
	// 	if (S.size() + C.size() > Ss.size()) {
	// 		Ss.clear();
	// 		for (int v : S) Ss.push(v);
	// 		for (int v : C) Ss.push(v);
	// 		return mode == RUSSIANDOLL_SEARCH;
	// 	}
	// 	return false;
	// }
	if (S.size() + C.size() <= Ss.size() || upperbound() <= Ss.size()) 
		return false;

	C1.clear();
	for (int v : C) {
		if (S.size() - degS[v] <= 1)
			C1.push(v);
	}

	// All connected
	int initPosC = C.frontPos();
	for (int v : C1) {
		if ((S.size() - degS[v]) + (C.size() - degC[v]) == 1) {
			moveCToS(v);
			C1.pop(v);
		}
	}

	if (C.size() == 0) branch(dep+1);
	else do {
		bool flagReturn = false;
		// 1 non-neighbor
		for (int v : C1) {
			if ((S.size() - degS[v]) + (C.size() - degC[v]) == 2) {
				int posC = update(v);
				if (branch(dep+1)) return true;
				restore(v, posC);
				flagReturn = true;
				break;
			}
		}

		if (flagReturn) break;

		// 2 non-neighbors
		for (int v : C1) {
			if ((S.size() - degS[v]) + (C.size() - degC[v]) == 3) {
				int posC = update(v);
				if (branch(dep+1)) return true;
				restore(v, posC);
				D.clear();
				for (int w : C) {
					if (v != w && S.size()-degS[w] <= 1 && !Sub.connect(v, w))
						D.push(w);
				}
				if (D.size() == 2) {
					int u = D[D.frontPos()], w = D[D.frontPos() + 1];
					if (2 * S.size() - degS[u] - degS[w] == 0 && Sub.connect(u, w)) {
						subC(v);
						for (int x : C)
							if (!Sub.connect(u, x) || !Sub.connect(w, x)) {
								subC(x);
							}
						moveCToS(u);
						moveCToS(w);
						if (branch(dep+1)) return true;
						moveSToC(w);
						moveSToC(u);
					}

				}
				else if (S.size() - degS[v] == 1 && D.size() == 1) {
					int u = D[D.frontPos()];
					if (S.size() == degS[u]) {
						subC(v);
						for (int w : C)
							if (!Sub.connect(u, w)) {
								subC(w);
							}
						moveCToS(u);
						if (branch(dep+1)) return true;
						moveSToC(u);
					}
				}					
				restoreC(posC);
				flagReturn = true;
				break;
			}
		}

		if (flagReturn) break;

		if (C.size() > C1.size()) { // Bipartite
			int u = C[C.frontPos()];
			for (int v : C) {
				if (degS[v] < degS[u])
					u = v;
			}
			int posC = update(u);
			if (branch(dep+1)) return true;
			restore(u, posC);

			subC(u);
			if (branch(dep+1)) return true;
			addC(u);
		}


		else { // Pivoting 
			int posC = C.frontPos();
			int u = C1[C1.frontPos()];
			for (int v : C1) {
				if (degC[v] > degC[u]) {
					u = v;
				}
			}

			std::vector<int> P1 = {u}, P2;
			bool flagNnbSu = S.size()-degS[u] == 1;
			for (int v : C) {
				if (v != u && !Sub.connect(u, v)) {
					if (flagNnbSu && S.size() == degS[v])
						P1.push_back(v);
					else
						P2.push_back(v);
				}
			}

			for (int v : P1) {
				int posC2 = update(v);
				if (branch(dep+1)) return true;
				restore(v, posC2);
				subC(v);
			}

			for (int v : P2) {
				int posC2 = update(v);
				for (int w : P2)
					if (C.inside(w) && Sub.connect(v, w)) {
						int posC3 = update(w);
						if (branch(dep+1)) return true;
						restore(w, posC3);
						subC(w);
					}
				restore(v, posC2);
				subC(v);
			}

			restoreC(posC);
		}


	} while (0);


	// backtrack all connected
	for (int i = C.frontPos()-1; i >= initPosC; --i)
		moveSToC(C[i]);

	return false;
}




void defclique::add(Graph &G, VertexSet &V, std::vector<int> &degV, int v) {
	if (V.inside(v)) return;
	V.push(v);
	for (int w : G.nbr[v])
		++degV[w];
}

void defclique::sub(Graph &G, VertexSet &V, std::vector<int> &degV, int v) {
	if (!V.inside(v)) return;
	V.pop(v);
	for (int w : G.nbr[v])
		--degV[w];
}

void defclique::addC(int v) {
	// nnbSub += S.size()-degS[v] + C.size()-degC[v];
	add(Sub, C, degC, v);
}

void defclique::subC(int v) {
	sub(Sub, C, degC, v);
	// nnbSub -= S.size()-degS[v] + C.size()-degC[v];
}
