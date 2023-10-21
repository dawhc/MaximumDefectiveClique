#ifndef ORDERING_HPP
#define ORDERING_HPP

#pragma once

#include "linearheap.hpp"
#include "graph.hpp"
//#include "coloring.hpp"
#include <vector>

class Ordering {
	LinearHeap vHeap;
public:
	std::vector<int> order, ordered, value;
	int numOrdered, capacity;
	Ordering(int size=0):capacity(0) {
		resize(size);
		numOrdered = 0;
	}
	void push(int u, int val) {
		ordered[numOrdered] = u;
		order[u] = numOrdered++;
		value[u] = val;
	}
	void resize(int size) {
		capacity = size;
		order.resize(size);
		ordered.resize(size);
		value.resize(size);
	}
	void degeneracyOrdering(Graph &G) {
		if (capacity < G.n) resize(G.n);
		vHeap.build(G.V, [&](int v) { return (int)G.nbr[v].size(); });
		numOrdered = 0;
		while (!vHeap.empty()) {
			int u = vHeap.top(); 
			vHeap.pop();
			for (int v : G.nbr[u]) {
				if (vHeap[v] <= vHeap[u]) continue;
				if (vHeap.inside(v)) vHeap.dec(v);
			}
			push(u, vHeap[u]);
		}
	}
	void degreeOrdering(Graph &G) {
		if (capacity < G.n) resize(G.n);
		vHeap.build(G.V, [&](int v) { return (int)G.nbr[v].size(); });
		numOrdered = 0;
		while (!vHeap.empty()) {
			int u = vHeap.top();
			vHeap.pop();
			push(u, vHeap[u]);
		}
	}
	static Ordering DegeneracyOrdering(Graph &G) {
		Ordering o;
		o.degeneracyOrdering(G);
		return o;
	}
	static Ordering DegreeOrdering(Graph &G) {
		Ordering o;
		o.degreeOrdering(G);
		return o;
	}
	// static Ordering colorOrdering(Graph &G, int tau) {
	// 	Ordering o(G.n);
	// 	Coloring c = Coloring::graphColoring(G, tau);
	// 	std::vector<int> cnt(G.n);
	// 	std::vector<bool> vis(c.numColors);
	// 	std::vector<std::vector<int>> bin(G.maxDeg + 1);
	// 	int maxValue = 0;
	// 	for (int u : G.V) {
	// 		for (int v : G.nbr[u])
	// 			if (!vis[c.color[v]]) {
	// 				vis[c.color[v]] = true;
	// 				maxValue = std::max(maxValue, ++cnt[u]);
	// 			}
	// 		bin[cnt[u]].push_back(u);
	// 		for (int v : G.nbr[u])
	// 			vis[c.color[v]] = false;
	// 	}

	// 	for (int cn = 0; cn <= maxValue; ++cn) {
	// 		for (int u : bin[cn])
	// 			o.push(u, cnt[u]);
	// 	}
	// 	return o;
	// }
};

#endif // ORDERING_HPP