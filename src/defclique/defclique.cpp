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


namespace defclique {
	int k, nnbS, nnbSub;
	Graph Sub;
	Coloring clr;
	VertexSet S, C, Ss, C1, D;
	std::vector<int> degC1, degC, degS, cnD;
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
	std::queue<int> q;
	std::vector<bool> vis(G.n);
	std::vector<int> deg(G.degree);

	if (k <= 1) return G;

	log("Running core reduction with k=%d...", k);

	for (int u : G.V) {
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

	log("Core reduction done!");
	log("Before: n=%d, m=%d; After: n=%d, m=%d", G.V.size(), G.m, C.V.size(), C.m);

	return C;
}

Graph defclique::subGraph(Graph &G, VertexSet S) {
	Graph Sub;
	for (int u : S) {
		if (S.size() < G.degree[u]) {
			for (int v : S)
				if (u < v && G.connect(u, v)) {
					Sub.addEdge(u, v);
				}
		}
		else {
			for (int v : G.nbr[u])
				if (u < v && S.inside(v)) {
					Sub.addEdge(u, v);
				}
		}
	}
	return Sub;
}

VertexSet defclique::heuristic(Graph &G, int k) {

	log("Running heuristic algorithm...");

	auto startTimePoint = std::chrono::steady_clock::now();

	S.clear();
	std::vector<int>(G.n).swap(degS);

	std::queue<int> q;

	Ordering o = Ordering::degeneracyOrdering(G);
	int u = o.ordered[o.numOrdered - 1];
	add(G, S, degS, u);

	// Expanding Ss with vertices in N2(v) in a degree ordering

	for (int v : G.nbr[u]) {
		C.push(v);
		for (int w : G.nbr[v]) 
			if (w != u)
				C.push(w);
	}

	std::vector<int> N2(C.begin(), C.end());

	auto compByDegree = [&] (const int& u, const int& v) {
		return G.degree[u] > G.degree[v];
	};

	auto compByCore = [&] (const int& u, const int& v) {
		return o.order[u] > o.order[v];
	};

	std::sort(N2.begin(), N2.end(), compByCore);

	for (int v : N2) {
		nnbS += S.size() - degS[v];
		if (nnbS > k) break;
		add(G, S, degS, v);
	}

	Ss = S;

	logSet(Ss, "Initial S*");

	Graph Core = coreReduction(G, Ss.size() - k);

	o = Ordering::degeneracyOrdering(Core);

	for (int i = o.numOrdered-1; i >= 0; --i) {
		int u = o.ordered[i];
		S.clear();
		C.clear();
		C1.clear();
		std::vector<int>(Core.n).swap(degS);
		std::vector<int>(Core.n).swap(degC);
		std::vector<int>(Core.n).swap(degC1);
		nnbS = nnbSub = 0;
		add(Core, S, degS, u);


		// Construct C1
		for (int v : Core.nbr[u]) {
			if (o.order[v] > i) {
				add(Core, C1, degC1, v);
			}
		}
		// Prune C1
		while (true) {
			bool flagBreak = true;
			for (int v : C1) {
				if (degC1[v] < Ss.size() - k) {
					flagBreak = false;
					sub(Core, C1, degC1, v);
				}
			}
			if (flagBreak) break;
		}

		// Construct C2
		for (int v : C1) {
			for (int w : Core.nbr[v]) {
				if (o.order[w] > i && !C1.inside(w) && degC1[w] >= Ss.size()-k) {
					add(Core, C, degC, w);
				}
			}
		}	

		Sub = subGraph(Core, S + C1 + C);
		Ordering oSub = Ordering::degeneracyOrdering(Sub);
		for (int v : Sub.V)
			// nnbSub += S.size()-degS[v] + C1.size()-degC1[v] + C.size()-degC[v] - 1;
			nnbSub += Sub.V.size() - Sub.degree[v] - 1;

		nnbSub >>= 1;

		int j = oSub.numOrdered;	

		while (nnbSub > k) {

			// Select a vertex from C1 cup C2 with maximum core number
			for (u = oSub.ordered[--j]; !C1.inside(u) && !C.inside(u); u = oSub.ordered[--j]);


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
		if (Ss.size() < S.size() + C1.size() + C.size()) {
			Ss = S + C1 + C;
		}
	}


	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Heuristic algorithm done! Time spent: %ld ms", duration.count());
	logSet(Ss, "S*");

	return Ss;
	
}

void defclique::russianDollSearch(Graph G, int k) {

	S.reserve(G.n);
	C.reserve(G.n);
	C1.reserve(G.n);
	D.reserve(G.n);

	degS.resize(G.n);
	degC.resize(G.n);
	degC1.resize(G.n);
	cnD.resize(G.n);

	defclique::k = k;
	Ss = heuristic(G, k);
	Graph Raw = G;
	G = coreReduction(G, Ss.size() - k);
	Ordering o = Ordering::degeneracyOrdering(G);
	std::queue<int> q;

	log("Running russian doll search...");

	auto startTimePoint = std::chrono::steady_clock::now();

	for (int i = o.numOrdered-1; i >= 0; --i) {
		S.clear();
		C.clear();
		C1.clear();
		std::vector<int>(G.n).swap(degS);
		std::vector<int>(G.n).swap(degC);
		std::vector<int>(G.n).swap(degC1);
		
		nnbS = nnbSub = 0;

		int u = o.ordered[i];

		add(G, S, degS, u);

		if (Ss.size() < k+1) {
			for (int j = i+1; j < o.numOrdered; ++j)
				add(G, C, degC, o.ordered[j]);
		}
		else {
			for (int v : G.nbr[u]) {
				if (o.order[v] > i) {
					add(G, C1, degC1, v);
				}
			}
			while (true) {
				bool flagBreak = true;
				for (int v : C1) {
					if (degC1[v] < Ss.size()-k) {
						flagBreak = false;
						sub(G, C1, degC1, v);
					}
				}
				if (flagBreak) break;
			}
			for (int v : C1) {
				add(G, C, degC, v);
				for (int w : G.nbr[v]) {
					if (o.order[w] > i && degC1[w] > Ss.size()-k) {
						add(G, C, degC, w);
					}
				}
			}
		}

		Sub = subGraph(G, S + C);
		clr = Coloring::graphColoring(Sub, Ss.size()-k+1);
		for (int v : Sub.V) 
			// nnbSub += S.size()-degS[v] + C.size()-degC[v]-1;
			nnbSub += Sub.V.size() - Sub.degree[v] - 1;

		nnbSub >>= 1;

		branch(0);
	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Russian doll search done! Time spent: %ld ms", duration.count());
	logSet(Ss, "S*");

	Graph subSs = subGraph(Raw, Ss);

	int cnt = 0;
	for (int v : subSs.V) {
		for (int w : subSs.V)
			if (v < w && !subSs.connect(v, w)) {
				log("Missing edge in S*: (%d, %d)", v, w);
				++cnt;
			}
	}
	log("Number of missing edges in S*: %d", cnt);

}

void defclique::reductionSearch(Graph G, int k) {

	S.reserve(G.n);
	C.reserve(G.n);
	C1.reserve(G.n);
	D.reserve(G.n);

	defclique::k = k;
	Ss = heuristic(G, k);
	G = coreReduction(G, Ss.size() - k);
	Ordering o = Ordering::degreeOrdering(G);
	Coloring c = Coloring::graphColoring(G, Ss.size()-k+1);
	std::vector<bool> vis(c.numColors);

	log("Starting branch and reduction search...");

	degS.resize(G.n);
	degC.resize(G.n);
	degC1.resize(G.n);
	cnD.resize(G.n);

	auto startTimePoint = std::chrono::steady_clock::now();
	
	for (int i = 0; i < G.V.size() - Ss.size(); ++i) {
		S.clear();
		C.clear();
		C1.clear();
		for (int v : G.V) degS[v] = degC[v] = degC1[v] = 0;
		nnbS = nnbSub = 0;

		int u = o.ordered[i];

		if (Ss.size() < k + 1) {
			add(G, S, degS, u);
			for (int v = 0; v < G.n; ++v)
				add(G, C, degC, v);
		}
		else {
			int deg = 0, cn = 0;
			for (int v : G.nbr[u]) {
				if (o.order[v] > i) {
					++deg;
					if (!vis[c.color[v]]) {
						vis[c.color[v]] = true;
						++cn;
					}
				}
			}
			for (int v : G.nbr[u])
				vis[c.color[v]] = false;

			if (deg > Ss.size()-k && cn >= Ss.size()-k) {
				add(G, S, degS, u);
				for (int v : G.nbr[u]) {
					if (o.order[v] > i)
						add(G, C1, degC1, v);
				}
				// for (int v : C1) {
				// 	if (degC1[v] < Ss.size() - k) {
				// 		sub(G, C1, degC1, v);
				// 		q.push(v);
				// 	}
				// }
				// while (!q.empty()) {
				// 	int v = q.front(); q.pop();
				// 	for (int w : G.nbr[v]) {
				// 		if (degC1[w] < Ss.size()-k && C1.inside(w)) {
				// 			sub(G, C1, degC1, w);
				// 			q.push(w);
				// 		}
				// 	}
				// }
				while (true) {
					bool flagBreak = true;
					for (int v : C1) {
						if (degC1[v] < Ss.size()-k) {
							flagBreak = false;
							sub(G, C1, degC1, v);
						}
					}
					if (flagBreak) break;
				}
				for (int v : C1) {
					add(G, C, degC, v);
					for (int w : G.nbr[v]) {
						if (o.order[w] > i && degC1[w] > Ss.size()-k) 
							add(G, C, degC, w);
					}
				}
			}
		}

		Sub = subGraph(G, S + C);
		clr = Coloring::graphColoring(Sub, Ss.size()-k+1);
		for (int v : Sub.V) 
			// nnbSub += S.size()-degS[v] + C.size()-degC[v]-1;
			nnbSub += Sub.V.size() - Sub.degree[v] - 1;
		nnbSub >>= 1;
		branch(0);

	}

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

	log("Branch and reduction search done! Time spent: %ld ms", duration.count());
	logSet(Ss, "S*");
}


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
	// int cntNnbSub = 0, cntNnbS = 0;
	// VertexSet V = S + C;
	// for (int u : V) {
	// 	for (int v : V) {
	// 		if (u < v && !Sub.connect(u, v)) {
	// 			++cntNnbSub;
	// 			if (S.inside(u) && S.inside(v))
	// 				++cntNnbS;
	// 		}
	// 	}
	// }
	// log("Depth=%d, |S|=%d, |C|=%d, |C1|=%d, nnbS=%d(real=%d), nnbSub=%d(real=%d)", dep, S.size(), C.size(), C1.size(), nnbS, cntNnbS, nnbSub, cntNnbSub);
	if (nnbSub <= k) {
		if (S.size() + C.size() > Ss.size()) {
			// log("New S*!");
			Ss = S + C;
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
	for (int v : C1) {
		if ((S.size() - degS[v]) + (C.size() - degC[v]) == 1) {
			moveCToS(v);
			C1.pop(v);
		}
	}

	// 1 non-neighbor
	for (int v : C1) {
		if ((S.size() - degS[v]) + (C.size() - degC[v]) == 2) {
			int posC = updateC(v);
			moveCToS(v);
			if (branch(dep+1)) return true;
			moveSToC(v);
			restoreC(posC);
			return false;
		}
	}

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
				if (S.size()-degS[w] <= 1 && !Sub.connect(v, w))
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
			return false;
		}
	}


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

		VertexSet Cand = C;
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
				if (S.size()-degS[v] <= 1 && !Sub.connect(u, v))
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
		return false;
	}

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
