#include "defclique.h"
#include "../utils/ordering.hpp"
#include "../utils/coloring.hpp"
#include "../utils/log.hpp"
#include <chrono>
#include <cstdint>
#include <ratio>
#include <vector>
#include <queue>
#include <cstdio>
#include <sstream>

// #define EDGE_REDUCTION


namespace defclique {
	int k, nnbS, nnbSub;
	Graph Sub;
	Coloring clr;
	VertexSet S, C, Ss, C1, D;
	std::vector<int> degC1, degC, degS, cnD;
	std::vector<int> q;
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

	std::queue<int> q;
	std::vector<bool> vis(G.n);
	std::vector<int> deg(G.n);

	for (int u : G.V) {
		deg[u] = G.nbr[u].size();
		if (deg[u] < k) {
			vis[u] = true;
			q.push(u);
		}
	}

	while (!q.empty()) {
		int u = q.front(); q.pop();
		for (int v : G.nbr[u]) {
			if (vis[v]) continue;
			if (--deg[v] < k) {
				vis[v] = true;
				q.push(v);
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

void defclique::preprocessing(Graph &G, Ordering &o, int u) {
	S.clear();
	C.clear();
	C1.clear();
	nnbS = nnbSub = 0;
	S.push(u);

	int i = o.order[u];

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

	nnbSub = ((Sub.V.size() * (Sub.V.size()-1)) >> 1) - Sub.m;

}

VertexSet defclique::heuristic(Graph &G, int k) {

	log("Running heuristic algorithm...");

	auto startTimePoint = std::chrono::steady_clock::now();

	S.clear();
	C.clear();
	std::vector<int>(G.n).swap(degS);

	// std::queue<int> q;
	Ordering o, oSub;

	o.degeneracyOrdering(G);
	int u = o.ordered[o.numOrdered - 1];
	add(G, S, degS, u);

	// Expanding Ss with vertices in N2(v) in a degree ordering
	for (int v : G.nbr[u]) {
		C.push(v);
		// for (int w : G.nbr[v]) 
		// 	if (!C.inside(w) && w != u)
		// 		C1.push(w);
	}

	Sub.subGraph(G, S, C);
	o.degeneracyOrdering(Sub);

	for (int i = o.numOrdered - 1; i >= 0; --i) {
		int v = o.ordered[i];
		if (!C.inside(v)) continue;
		nnbS += S.size() - degS[v];
		if (nnbS > k) break;
		add(Sub, S, degS, v);
	}

	Ss = S;

	logSet(Ss, "Initial S*");

	Graph Core = coreReduction(G, Ss.size() - k);
	o.degeneracyOrdering(Core);

	for (int i = o.numOrdered-1; i >= 0; --i) {
		// fprintf(stderr, "%d/%d\r", o.numOrdered-i, o.numOrdered);

		int u = o.ordered[i];
		if (o.value[u] < Ss.size()-k) break;

		preprocessing(Core, o, u);

		oSub.degeneracyOrdering(Sub);

		int j = oSub.numOrdered;	

		while (nnbSub > k) {

			// Select a vertex from C1 cup C2 with maximum core number
			int u = oSub.ordered[--j];
			for (; !C1.inside(u) && !C.inside(u); u = oSub.ordered[--j]);


			// Add u from C1|C2 to S
			nnbS += S.size() - degS[u];
			add(Sub, S, degS, u);

			C1.inside(u) ? sub(Sub, C1, degC1, u) : sub(Sub, C, degC, u);

			// Prune C1 & C2 
			for (int v : C1) {
				if (nnbS + S.size()-degS[v] > k || degS[v]+degC1[v] < Ss.size()-k) {
					sub(Sub, C1, degC1, v);
					nnbSub -= (S.size()-degS[v] + C1.size()-degC1[v] + C.size()-degC[v]);
				}
			}

			for (int v : C) {
				if (nnbS + S.size()-degS[v] > k || degS[v]+degC1[v] < Ss.size()-k) {
					sub(Sub, C, degC, v);
					nnbSub -= (S.size()-degS[v] + C1.size()-degC1[v] + C.size()-degC[v]);
				}
			}

			bool flagBreak = false;

			for (int v : S) {
				if (degS[v]+degC1[v] < Ss.size()-k) {
					C1.clear();
					C.clear();
					flagBreak = true;
					break;
				}
			}

			if (flagBreak) break;
			
		}


		// Update Ss
		if (Ss.size() < S.size() + C.size() + C1.size()) {
			Ss.clear();
			for (int v : S) Ss.push(v);
			for (int v : C) Ss.push(v);
			for (int v : C1) Ss.push(v);
		}
	}


	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Heuristic algorithm done! Time spent: %ld ms", duration.count());
	logSet(Ss, "S*");

	return Ss;
	
}

void defclique::russianDollSearch(Graph &G, int k) {
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

	defclique::k = k;

	Sub.nbrMap = G.nbrMap;

	Ss = heuristic(G, k);
	// return;
	// for (int i = 0; i < 26; ++i) Ss.push(i);
	Graph Core = coreReduction(G, Ss.size() - k);
#ifdef EDGE_REDUCTION
	Core = edgeReduction(Core, Ss.size() - k - 1);
#endif
	Ordering o = Ordering::DegeneracyOrdering(Core);
	// std::queue<int> q;

	log("Running russian doll search...");

	auto startTimePoint = std::chrono::steady_clock::now();

	for (int i = o.numOrdered-1; i >= 0; --i) {

		int u = o.ordered[i];
		if (o.value[u] < Ss.size()-k) break;

		if (Ss.size() < k+1) {
		// if (true) {
			S.clear();
			C.clear();
			S.push(u);
			nnbS = nnbSub = 0;
			for (int j = i+1; j < o.numOrdered; ++j)
				C.push(o.ordered[j]);
				// add(G, C, degC, o.ordered[j]);
			Sub.subGraph(Core, S, C);
			nnbSub = ((Sub.V.size() * (Sub.V.size()-1)) >> 1) - Sub.m;
			for (int v : Sub.V)
				degS[v] = degC[v] = 0;
			for (int v : Sub.nbr[u])
				degS[v] = 1;
			for (int v : C) {
				for (int w : Sub.nbr[v])
					++degC[w];
			}
		}
		else {
			preprocessing(Core, o, u);
			for (int v : C1)
				C.push(v);
			for (int v : Sub.V)
				degC[v] += degC1[v];
		}
		// clr = Coloring::graphColoring(Sub, Ss.size()-k+1);
		clr.graphColoring(Sub, Ss.size()-k+1);
		branch(0);
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Russian doll search done! Time spent: %ld ms", duration.count());
	logSet(Ss, "S*");

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

}

// void defclique::reductionSearch(Graph G, int k) {

// 	S.reserve(G.n);
// 	C.reserve(G.n);
// 	C1.reserve(G.n);
// 	D.reserve(G.n);

// 	defclique::k = k;
// 	Ss = heuristic(G, k);
// 	G = coreReduction(G, Ss.size() - k);
// 	Ordering o = Ordering::degreeOrdering(G);
// 	Coloring c = Coloring::graphColoring(G, Ss.size()-k+1);
// 	std::vector<bool> vis(c.numColors);

// 	log("Starting branch and reduction search...");

// 	degS.resize(G.n);
// 	degC.resize(G.n);
// 	degC1.resize(G.n);
// 	cnD.resize(G.n);

// 	auto startTimePoint = std::chrono::steady_clock::now();
	
// 	for (int i = 0; i < G.V.size() - Ss.size(); ++i) {
// 		S.clear();
// 		C.clear();
// 		C1.clear();
// 		for (int v : G.V) degS[v] = degC[v] = degC1[v] = 0;
// 		nnbS = nnbSub = 0;

// 		int u = o.ordered[i];

// 		if (Ss.size() < k + 1) {
// 			addG(G, S, degS, u);
// 			for (int v = 0; v < G.n; ++v)
// 				addG(G, C, degC, v);
// 		}
// 		else {
// 			int deg = 0, cn = 0;
// 			for (int v : G.nbr[u]) {
// 				if (o.order[v] > i) {
// 					++deg;
// 					if (!vis[c.color[v]]) {
// 						vis[c.color[v]] = true;
// 						++cn;
// 					}
// 				}
// 			}
// 			for (int v : G.nbr[u])
// 				vis[c.color[v]] = false;

// 			if (deg > Ss.size()-k && cn >= Ss.size()-k) {
// 				addG(G, S, degS, u);
// 				for (int v : G.nbr[u]) {
// 					if (o.order[v] > i)
// 						addG(G, C1, degC1, v);
// 				}
// 				// for (int v : C1) {
// 				// 	if (degC1[v] < Ss.size() - k) {
// 				// 		sub(G, C1, degC1, v);
// 				// 		q.push(v);
// 				// 	}
// 				// }
// 				// while (!q.empty()) {
// 				// 	int v = q.front(); q.pop();
// 				// 	for (int w : G.nbr[v]) {
// 				// 		if (degC1[w] < Ss.size()-k && C1.inside(w)) {
// 				// 			sub(G, C1, degC1, w);
// 				// 			q.push(w);
// 				// 		}
// 				// 	}
// 				// }
// 				while (true) {
// 					bool flagBreak = true;
// 					for (int v : C1) {
// 						if (degC1[v] < Ss.size()-k) {
// 							flagBreak = false;
// 							subG(G, C1, degC1, v);
// 						}
// 					}
// 					if (flagBreak) break;
// 				}
// 				for (int v : C1) {
// 					addG(G, C, degC, v);
// 					for (int w : G.nbr[v]) {
// 						if (o.order[w] > i && degC1[w] > Ss.size()-k) 
// 							addG(G, C, degC, w);
// 					}
// 				}
// 			}
// 		}

// 		Sub.subGraph(G, S + C);
// 		clr = Coloring::graphColoring(Sub, Ss.size()-k+1);
// 		for (int v : Sub.V) 
// 			// nnbSub += S.size()-degS[v] + C.size()-degC[v]-1;
// 			nnbSub += Sub.V.size() - Sub.nbr[v].size() - 1;
// 		nnbSub >>= 1;
// 		branch(0);

// 	}

// 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
// 		std::chrono::steady_clock::now() - startTimePoint);

// 	log("Branch and reduction search done! Time spent: %ld ms", duration.count());
// 	logSet(Ss, "S*");
// }


int defclique::upperbound() {
	D.clear();
	int s = nnbS;
	int posC = C.frontPos();
	for (int v : C) cnD[v] = 0;
	while (C.size() > 0) {
		int u = -1, minv = 0x7fffffff;
		for (int v : C) {
			if (cnD[v] + S.size() - degS[v] < minv) {
				minv = cnD[v] + S.size() - degS[v];
				u = v;
			}
		}

		if ((s += minv) > k) break;

		D.push(u);
		C.pop(u);
		for (int v : C) {
			if (clr.color[v] == clr.color[u]) 
				++cnD[v];
		}
	}
	C.restore(posC);
	return S.size() + D.size();
}


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
			nnbSub -= S.size()-degS[u] + C.size()-degC[u];
		}
	}
	return posC;
}

void defclique::restoreC(int pos) {
	for (int i = C.frontPos()-1; i >= pos; --i) {
		int u = C[i];
		nnbSub += S.size()-degS[u] + C.size()-degC[u];
		add(Sub, C, degC, u);
	}
}


bool defclique::branch(int dep) {
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
	// if (flagDebug)  {
	// 	log("\n****** Depth=%d, |S|=%d, |C|=%d, |C1|=%d, nnbS=%d(real=%d), nnbSub=%d(real=%d)\n", dep, S.size(), C.size(), C1.size(), nnbS, cntNnbS, nnbSub, cntNnbSub);
	// 	logSet(S, "S");
	// 	logSet(C, "C");	
	// }

	if (nnbSub <= k) {
		if (S.size() + C.size() > Ss.size()) {
			Ss.clear();
			for (int v : S) Ss.push(v);
			for (int v : C) Ss.push(v);
			return true;
		}
		return false;
	}
	if (upperbound() <= Ss.size()) 
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

	bool flagReturn = false;

	// 1 non-neighbor
	for (int v : C1) {
		if ((S.size() - degS[v]) + (C.size() - degC[v]) == 2) {
			int posC = updateC(v);
			moveCToS(v);
			if (branch(dep+1)) return true;
			moveSToC(v);
			restoreC(posC);
			flagReturn = true;
			break;
		}
	}

	if (!flagReturn) {
		// 2 non-neighbors
		for (int v : C1) {
			if ((S.size() - degS[v]) + (C.size() - degC[v]) == 3) {
				int posC = updateC(v);
				moveCToS(v);
				if (branch(dep+1)) return true;
				moveSToC(v);
				restoreC(posC);
				D.clear();
				for (int w : C) {
					if (v != w && S.size()-degS[w] <= 1 && !Sub.connect(v, w))
						D.push(w);
				}
				if (D.size() == 2) {
					int u = D[D.frontPos()], w = D[D.frontPos() + 1];
					if (2 * S.size() - degS[u] - degS[w] == 0 && Sub.connect(u, w)) {
						sub(Sub, C, degC, v);
						nnbSub -= (S.size()-degS[v] + C.size()-degC[v]);
						for (int x : C)
							if (!Sub.connect(u, x) || !Sub.connect(w, x)) {
								sub(Sub, C, degC, x);
								nnbSub -= (S.size()-degS[x] + C.size()-degC[x]);
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
						sub(Sub, C, degC, v);
						nnbSub -= (S.size()-degS[v] + C.size()-degC[v]);
						for (int w : C)
							if (!Sub.connect(u, w)) {
								sub(Sub, C, degC, w);
								nnbSub -= (S.size()-degS[w] + C.size()-degC[w]);
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
	}

	if (!flagReturn) {
		if (C.size() > C1.size()) {
			int u = C[C.frontPos()];
			for (int v : C) {
				if (degS[v] < degS[u])
					u = v;
			}
			int posC = updateC(u);
			moveCToS(u);
			if (branch(dep+1)) return true;
			moveSToC(u);
			restoreC(posC);

			sub(Sub, C, degC, u);
			nnbSub -= S.size()-degS[u] + C.size()-degC[u];
			if (branch(dep+1)) return true;
			nnbSub += S.size()-degS[u] + C.size()-degC[u];
			add(Sub, C, degC, u);
		}


		else {
			int posC = C.frontPos();
			int u = C1[C1.frontPos()];
			for (int v : C1) {
				if (degS[v] + degC[v] > degS[u] + degC[u]) {
					u = v;
				}
			}

			// VertexSet Cand = C;
			std::vector<int> Cand(C.begin(), C.end());
			for (int v : Cand) {
				if (S.size()-degS[u] + C.size()-degC[u] == 3) break;
				if (v != u && !Sub.connect(u, v)) {
					C1.pop(v);
					int posC2 = updateC(v);
					moveCToS(v);
					if (branch(dep+1)) return true;
					moveSToC(v);
					restoreC(posC2);
					sub(Sub, C, degC, v);
					nnbSub -= S.size()-degS[v] + C.size()-degC[v];
				}
			}
			if (S.size()-degS[u] + C.size()-degC[u] == 3) {
				int posC2 = updateC(u);
				moveCToS(u);
				if (branch(dep+1)) return true;
				moveSToC(u);
				restoreC(posC2);
				D.clear();
				for (int v : C) {
					if (u != v && S.size()-degS[v] <= 1 && !Sub.connect(u, v))
						D.push(v);
				}
				if (D.size() == 2) {
					int v = D[D.frontPos()], w = D[D.frontPos() + 1];
					if (2 * S.size() - degS[v] - degS[w] == 0 && Sub.connect(v, w)) {
						sub(Sub, C, degC, u);
						nnbSub -= (S.size()-degS[u] + C.size()-degC[u]);
						for (int x : C)
							if (!Sub.connect(v, x) || !Sub.connect(w, x)) {
								sub(Sub, C, degC, x);
								nnbSub -= S.size()-degS[x] + C.size()-degC[x];
							}
						moveCToS(v);
						moveCToS(w);
						if (branch(dep+1)) return true;
						moveSToC(w);
						moveSToC(v);
					}

				}
				else if (S.size() - degS[u] == 1 && D.size() == 1) {
					int v = D[D.frontPos()];
					if (S.size() == degS[v]) {
						sub(Sub, C, degC, u);
						nnbSub -= (S.size()-degS[u] + C.size()-degC[u]);
						for (int w : C)
							if (!Sub.connect(v, w)) {
								sub(Sub, C, degC, w);
								nnbSub -= S.size()-degS[w] + C.size()-degC[w];
							}
						moveCToS(v);
						if (branch(dep+1)) return true;
						moveSToC(v);
					}
				}
			}

			restoreC(posC);
		}
	}

	

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

void defclique::addSub(Graph &G, VertexSet &Sub, VertexSet &V, std::vector<int> &degV, int v) {
	if (V.inside(v)) return;
	V.push(v);
	if (Sub.size() < G.nbr[v].size()) {
		for (int w : Sub)
			if (G.connect(v, w))
				++degV[w];
	}
	else {
		for (int w : G.nbr[v])
			if (Sub.inside(w))
				++degV[w];
	}
}

void defclique::subSub(Graph &G, VertexSet &Sub, VertexSet &V, std::vector<int> &degV, int v) {
	if (!V.inside(v)) return;
	V.pop(v);
	if (Sub.size() < G.nbr[v].size()) {
		for (int w : Sub)
			if (G.connect(v, w))
				--degV[w];
	}
	else {
		for (int w : G.nbr[v])
			if (Sub.inside(w))
				--degV[w];
	}
}
