#include "kdbb.h"
#include "pmc/pmc.h"
#include <chrono>
#include <chrono>
#include <omp.h>
#include <queue>
#include <algorithm>

namespace kdbb {
	VertexSet S, C;
	std::vector<int> degS, degC;
	int k, nnbS, lb, numBranches, numBound;
	Graph G;
	std::vector<int> bin;
}


Graph kdbb::preprocessing(Graph &G, int k, int lb) {
	fprintf(stderr, "Preprocessing......");
	Graph C = coreReduction(G, lb-k);
	C = edgeReduction(C, lb-k-1);
	fprintf(stderr, "Done, before: n=%d, m=%d; after: n=%d, m=%d\n", G.V.size(), G.m, C.V.size(), C.m);
	return C;
}

Graph kdbb::coreReduction(Graph& G, int k) {
	if (k <= 1) return G;

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

	return C;
}



Graph kdbb::edgeReduction(Graph &G, int k) {
	std::vector<int> cn(G.m), q(G.m);
	std::vector<std::pair<int, int>> edges;
	std::vector<std::unordered_map<int, int>> eid(G.n);	
	std::vector<bool> vis(G.m);

	edges.reserve(G.m);

	int head = 0, tail = 0;

	auto countCommonNeighbor = [&](int u, int v) {
		int cnt = 0;
		if (G.nbr[u].size() > G.nbr[v].size()) std::swap(u, v);
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

	return E;

}


int kdbb::fastLB(std::string filename) {
	auto startTimePoint = std::chrono::steady_clock::now();
    pmc::pmc_graph G(filename);
    input in;
    in.graph = filename;
    in.threads = std::min(8, omp_get_max_threads());
    G.compute_cores();
    in.ub = G.get_max_core() + 1;
    std::vector<int> C;
    pmc::pmc_heu maxclique(G, in);
    in.lb = maxclique.search(G, C);
    if (in.lb == in.ub) {
	    auto duration = std::chrono::duration_cast<chrono::milliseconds>(
			std::chrono::steady_clock::now() - startTimePoint);
    	fprintf(stderr, "PMC result: size=%d, time=%ld ms\n", in.lb, duration.count());
    	return in.lb;
    }
    if (G.num_vertices() < in.adj_limit) {
        G.create_adj();
        pmc::pmcx_maxclique finder(G,in);
        finder.search_dense(G,C);
    }
    else {
        pmc::pmcx_maxclique finder(G,in);
        finder.search(G,C);
    }
    auto duration = std::chrono::duration_cast<chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);
    fprintf(stderr, "PMC result: size=%d, time=%ld ms\n", C.size(), duration.count());
    return C.size();
}


int kdbb::run(std::string filename, int k) {
	// TODO: FastLB
	Graph inputG(filename);
	lb = fastLB(filename);
	lb = std::max(lb, k+1);
	kdbb::k = k;
	G = preprocessing(inputG, k, lb);
	S.reserve(G.n);
	C.reserve(G.n);
	degS.resize(G.n);
	degC.resize(G.n);
	bin.resize(G.maxDeg+1);
	auto startTimePoint = std::chrono::steady_clock::now();
	for (int v : G.V) C.push(v);
	nnbS = numBranches = numBound = 0;
	branch(0, -1);
	auto duration = std::chrono::duration_cast<chrono::milliseconds>(
		std::chrono::steady_clock::now() - startTimePoint);
	fprintf(stderr, "KDBB result: size=%d, time=%ld ms, numBranches=%d, numBound=%d\n", lb, duration.count(), numBranches, numBound);
	return lb;
}

void printSet(VertexSet &V, const std::string &name) {
	std::vector<int> S(V.begin(), V.end());
	std::sort(S.begin(), S.end());
	std::stringstream ss;
	ss << "{";
	for (int i = 0; i < S.size(); ++i) {
		if (i > 0) ss << ",";
		ss << S[i];
	}
	ss << "}";
	fprintf(stderr, "%s: size=%d, content=%s\n", 
		name.c_str(), V.size(), ss.str().c_str());
}

void kdbb::branch(int dep, int u) {
	++numBranches;

	auto numNbrS = [&](int v) {
		int cnt = 0;
		if (S.size() < G.nbr[v].size()) {
			for (int w : S)
				if (G.connect(w, v))
					++cnt;
		}
		else {
			for (int w : G.nbr[v])
				if (S.inside(w) && G.connect(v, w))
					++cnt;
		}
		return cnt;
	};

	auto numNbrC = [&](int v) {
		int cnt = 0;
		if (C.size() < G.nbr[v].size()) {
			for (int w : C)
				if (G.connect(w, v))
					++cnt;
		}
		else {
			for (int w : G.nbr[v])
				if (C.inside(w) && G.connect(v, w))
					++cnt;
		}
		return cnt;
	};

	auto numCommNbrC = [&](int u, int v) {
		int cnt = 0;
		if (G.nbr[u].size() > G.nbr[v].size()) std::swap(u, v);
		if (G.nbr[u].size() < C.size()) {
			for (int w : G.nbr[u])
				if (C.inside(w) && G.connect(u, w) && G.connect(v, w))
					++cnt;
		}
		else {
			for (int w : C)
				if (G.connect(u, w) && G.connect(v, w))
					++cnt;
		}
		return cnt;
	};

	auto candibound = [&]() {
		int cb = S.size(), maxNonDeg = 0, nnbCnt = nnbS;
		for (int v : C) {
			maxNonDeg = std::max(maxNonDeg, S.size()-degS[v]+1);
			++bin[S.size()-degS[v]];
		}

		for (int i = 0; i < maxNonDeg; ++i) {
			if (bin[i]*i + nnbCnt <= k) {
				cb += bin[i];
				nnbCnt += bin[i] * i;
			}
			else {
				cb += (k-nnbCnt) / i;
				break;
			}
		}

		for (int i = 0; i < maxNonDeg; ++i)
			bin[i] = 0;
		return cb;
	};

	// fprintf(stderr, "*** dep=%d, |S|=%d, |C|=%d, nnbS=%d, u=%d, lb=%d\n", 
	// 	dep, S.size(), C.size(), nnbS, u, lb);

	// printSet(S, "S");
	// printSet(C, "C");


	if (nnbS > k) return;

	int posC = C.frontPos();

	std::vector<std::pair<int, int>> removedEdges;

	for (int v : C) {
		degS[v] = numNbrS(v);
		degC[v] = numNbrC(v);
	}

	// prune C
	if (u != -1) {
		if (S.inside(u)) {
			for (int v : C) {
				int comm = numCommNbrC(u, v);
				if (S.size()+1 + comm + std::min(k-nnbS-(S.size()-degS[v]), C.size()-comm-1) <= lb)
					C.pop(v);
			}
		}

		for (int w : C) {
			bool flag = false;
			if (G.connect(u, w)) flag = true;
			else {
				for (int x : G.nbr[w]) if (!S.inside(x) && !C.inside(x)) {
					flag = true;
					break;
				}
			}
			if (flag) {
			// if (true) {
				if (S.size()+1 + degC[w] + std::min(k-nnbS-(S.size()-degS[w]), C.size()-degC[w]-1) <= lb) {
					C.pop(w);
					continue;
				}

				if (G.nbr[w].size() < C.size()) {	
					for (int u : G.nbr[w]) if (C.inside(u) && G.connect(u, w)) {
						int comm = numCommNbrC(u, w);
						if (S.size()+2 + comm + std::min(k-nnbS-(2*S.size()-degS[u]-degS[w]), C.size()-comm-2) <= lb) {
							removedEdges.push_back(make_pair(u, w));
							G.nbrMap[u].erase(w);
							G.nbrMap[w].erase(u);
						} 
					}
				}
				else {
					for (int u : C) if (G.connect(u, w)) {
						int comm = numCommNbrC(u, w);
						if (S.size()+2 + comm + std::min(k-nnbS-(2*S.size()-degS[u]-degS[w]), C.size()-comm-2) <= lb) {
							removedEdges.push_back(make_pair(u, w));
							G.nbrMap[u].erase(w);
							G.nbrMap[w].erase(u);
						}
					}
				}
			}
		}
	}

	do {

		if (C.size() == 0) {
			lb = std::max(lb, S.size());
			break;
		}

		if (candibound() <= lb) {
			++numBound;
			break;
		}

		int v = C[C.frontPos()];
		int degSv = degS[v];
		
		nnbS += S.size() - degSv;
		C.pop(v);
		S.push(v);
		branch(dep+1, v);
		S.pop(v);
		nnbS -= S.size() - degSv;
		branch(dep+1, v);
		C.push(v);

	} while (false);

	C.restore(posC);

	for (auto e : removedEdges) {
		G.nbrMap[e.first].insert(e.second);
		G.nbrMap[e.second].insert(e.first);
	}
}
