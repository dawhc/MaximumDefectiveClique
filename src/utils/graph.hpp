#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cstdint>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "hash.hpp"
#include "fastio.hpp"
#include "vertexset.hpp"

struct Graph {
	int n, m, maxDeg;
	std::vector<CuckooHash> nbrMap;
	std::vector<std::vector<int>> nbr;
	VertexSet V;

	Graph() {
		n = m = maxDeg = 0;
	}

	Graph(int n): n(n) {
		m = maxDeg = 0;
		resize(n);
	}

	Graph(const std::string& dataset) {
		loadFromFile(dataset);
	}

	void clear() {
		std::vector<std::vector<int>>().swap(nbr);
		std::vector<CuckooHash>().swap(nbrMap);
		V.clear();
		n = m = maxDeg = 0;
	}

	void resize(int size) {
		nbr.resize(size);
		nbrMap.resize(size);
		V.reserve(size);
	}

	void loadFromFile(const std::string& filename) {
		clear();
		FastIO fio(filename, "r");
		n = fio.getUInt();
		int numEdges = fio.getUInt();
		resize(n);
		for (int i = 0; i < numEdges; ++i) {
			int u = fio.getUInt();
			int v = fio.getUInt();
			addEdge(u, v);
		}
	}

	void addEdge(int u, int v) {
		if (std::max(u, v) >= n) {
			n = std::max(u, v) + 1;
			resize(n);
		} 
		V.push(u);
		V.push(v);
		if (nbrMap[u].find(v)) return;
		++m;
		nbr[u].push_back(v);
		nbr[v].push_back(u);
		nbrMap[u].insert(v);
		nbrMap[v].insert(u);
		maxDeg = std::max(maxDeg, (int)nbr[u].size());
		maxDeg = std::max(maxDeg, (int)nbr[v].size());
	}

	bool connect(int u, int v) const {
		return nbrMap[u].find(v);
	}
};

struct LinkedGraph {

	struct Edge {
		int u, v, next, prev;
	};
	
	int n, m;

	std::vector<Edge> edges;
	std::vector<int> first;
	std::vector<CuckooHash> nbrMap;

	LinkedGraph(int n, int m = 0): n(n), m(m), first(std::vector<int>(n, -1)) {
		edges.reserve(m);
		m = 0;
	}

	void addEdge(int u, int v) {
		int w = std::max(u, v);
		if (w > first.size())
			first.resize(w + 1);
		edges.push_back((Edge){u, v, first[u], -1});
		first[u] = edges[first[u]].prev = edges.size() - 1;
		edges.push_back((Edge){v, u, first[v], -1});
		first[v] = edges[first[v]].prev = edges.size() - 1;
		nbrMap[u].insert(v);
		nbrMap[v].insert(u);
		++m;
	}

	void removeEdge(int eid) {
		if (eid & 1) --eid;
		Edge &e1 = edges[eid], &e2 = edges[eid+1];
		nbrMap[e1.u].remove(e1.v);
		nbrMap[e2.u].remove(e2.v);
		if (e1.next != -1) edges[e1.next].prev = e1.prev;
		if (e1.prev != -1) edges[e1.prev].next = e1.next;
		if (e2.next != -1) edges[e2.next].prev = e2.prev;
		if (e2.prev != -1) edges[e2.prev].next = e2.next;
	}

	bool connect(int u, int v) const {
		return nbrMap[u].find(v);
	}

};

#endif