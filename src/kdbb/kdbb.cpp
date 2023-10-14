#include "kdbb.h"
#include <vector>
#include <queue>
#include <algorithm>

namespace kdbb {
	VertexSet S, C;
	std::vector<int> degS, degC;
	int k, nnbS, nnbSub;
	std::vector<std::vector<int>> bin;
}

int kdbb::upperbound() {
	for (int i = 0; i <= k; ++i) bin[i].clear();
	for (int v : C) bin[S.size() - degS[v]].push_back(v);
	int nnbCnt = nnbS, vCnt = 0;
	for (int i = 0; i <= k; ++i) {
		for (int v : bin[i]) {
			if (nnbCnt + i > k) {
				return S.size() + vCnt;
			}
			else {
				nnbCnt += i;
				++vCnt;
			}
		}
	}
	return S.size() + vCnt;
}

Graph kdbb::preprocessing(Graph &G, int k, int lb) {
	Graph C = coreReduction(G, lb-k);
	C = edgeReduction(C, lb-k-1);
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

	auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
		std::chrono::steady_clock::now() - startTimePoint);

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

	return E;

}


int kdbb::KDBB(Graph &G, int k) {
	// TODO: FastLB
	int lb = FastLB();
	G = preprocessing(G, k, lb);

}

void kdbb::branch(int v) {
	
}