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
	std::vector<int> degree;
	std::vector<CuckooHash> nbrMap;
	std::vector<std::vector<int>> nbr;
	VertexSet V;

	Graph() {
		n = m = maxDeg = 0;
	}

	Graph(int size) {
		n = m = maxDeg = 0;
		resize(size);
	}

	Graph(const std::string& dataset) {
		loadFromFile(dataset);
	}

	void clear() {
		std::vector<int>().swap(degree);
		std::vector<std::vector<int>>().swap(nbr);
		std::vector<CuckooHash>().swap(nbrMap);
		V.clear();
		n = m = maxDeg = 0;
	}

	void resize(int size) {
		degree.resize(size);
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
		++degree[u];
		++degree[v];
		maxDeg = std::max(maxDeg, degree[u]);
		maxDeg = std::max(maxDeg, degree[v]);
	}

	bool connect(int u, int v) const {
		return nbrMap[u].find(v);
	}
};

#endif