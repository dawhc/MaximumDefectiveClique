#ifndef BIGRAPH_HPP
#define BIGRAPH_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

#include "hash.hpp"
#include "fastio.hpp"


struct BiGraph {
	uint32_t nLabels[2], n[2], m;
	uint32_t maxDeg[2];

	std::vector<CuckooHash> nbrMap[2];
	std::vector<uint32_t> degree[2];
	std::vector<std::vector<uint32_t>> nbr[2];	
	std::unordered_map<uint32_t, uint32_t> labelId[2];
	std::vector<uint32_t> vLabel[2];

	BiGraph(const std::string& dataset) {
		clear();
		loadFromFile(dataset);
	}

	BiGraph() {
		clear();
	}

	void clear() {
		m = 0;
		for (uint32_t i = 0; i <= 1; ++i) {
			maxDeg[i] = nLabels[i] = n[i] = 0;
			degree[i].clear();
			nbr[i].clear();
			nbrMap[i].clear();
			vLabel[i].clear();
			labelId[i].clear();
		}
	}

	uint32_t getIndex(uint32_t label, uint32_t side) {
		if (labelId[side].count(label)) return labelId[side][label];
		vLabel[side].push_back(label);
		degree[side].push_back(0);
		nbr[side].push_back(std::vector<uint32_t>());
		nbrMap[side].push_back(CuckooHash());
		return labelId[side][label] = n[side]++;
	}

	void loadFromFile(const std::string& filename) {
		FastIO fio(filename, "r");
		uint32_t m = fio.getUInt();
		nLabels[0] = fio.getUInt();
		nLabels[1] = fio.getUInt();

		for (uint32_t i = 0; i < m; ++i) {
			uint32_t u = fio.getUInt();
			uint32_t v = fio.getUInt();
			addEdgeWithLabel(u, v);
		}

	}

	void addEdgeWithLabel(uint32_t ulabel, uint32_t vlabel) {
		nLabels[0] = std::max(nLabels[0], ulabel + 1);
		nLabels[1] = std::max(nLabels[1], vlabel + 1);
		addEdge(getIndex(ulabel, 0), getIndex(vlabel, 1));
	}

	void addEdge(uint32_t u, uint32_t v) {
		if (nbrMap[0][u].find(v)) return;
		++m;
		nbr[0][u].push_back(v);
		nbr[1][v].push_back(u);
		nbrMap[0][u].insert(v);
		nbrMap[1][v].insert(u);
		++degree[0][u];
		++degree[1][v];
		maxDeg[0] = std::max(maxDeg[0], degree[0][u]);
		maxDeg[1] = std::max(maxDeg[1], degree[1][v]);
	}

	bool connect(uint32_t uSide, uint32_t u, uint32_t v) const {
		return nbrMap[uSide][u].find(v);
	}

};

#endif
