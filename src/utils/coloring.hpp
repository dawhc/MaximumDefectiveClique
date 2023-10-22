#ifndef COLORING_HPP
#define COLORING_HPP

#pragma once

#include "graph.hpp"
#include "ordering.hpp"
#include <cstdint>
#include <vector>

constexpr int uncolored = -1;

class Coloring {
	std::vector<bool> vis;
	std::vector<std::vector<int>> bin;
public:
	std::vector<int> color;
	int numColors, capacity;
	Ordering o;

	Coloring(int size): numColors(0) { resize(size); }

	Coloring(): capacity(0), numColors(0) {}

	void colorVertex(int v, int c) {
		color[v] = c;
		numColors = std::max(numColors, c+1);
	}

	void resize(int newCapacity) {
		capacity = newCapacity;
		color.resize(newCapacity, uncolored);
		vis.resize(newCapacity, false);
	}

	void graphColoring(Graph &G, int tau) {
		if (tau < 0) tau = 0;

		if (tau > bin.size()) bin.resize(tau);

		int n = 0;

		for (int v : G.V) {
			if (v < capacity) color[v] = uncolored;
			n = std::max(n, v+1);
		}

		if (n > capacity) resize(n);

		o.degeneracyOrdering(G);

		for (int i = o.numOrdered - 1; i >= 0; --i) {
			int u = o.ordered[i];

			// Coloring
			for (int v : G.nbr[u]) {
				if (color[v] != uncolored)
					vis[color[v]] = true;
			}
			for (int i = 0; i <= G.nbr[u].size(); ++i) {
				if (!vis[i]) {
					colorVertex(u, i);
					break;
				}
			}
			for (int v : G.nbr[u]) {
				if (color[v] != uncolored)
					vis[color[v]] = false;
			}

			// Recoloring
			if (color[u] >= tau) {
				for (int i = 0; i < tau; ++i) bin[i].clear();
				for (int v : G.nbr[u]) {
					if (color[v] != uncolored && color[v] < tau)
						bin[color[v]].push_back(v);
				}
				for (int i = 0; i < tau; ++i) {
					if (bin[i].size() == 1) {
						int v = bin[i][0];
						bool flag = false;
						for (int w : G.nbr[v]) {
							if (color[w] != uncolored && color[w] < tau)
								vis[color[w]] = true;
						}
						for (int j = 0; j < tau; ++j) {
							if (!vis[j] && j != i) {
								flag = true;
								colorVertex(v, j);
								colorVertex(u, i);
								break;
							}
						}
						for (int w : G.nbr[v]) {
							if (color[w] != uncolored && color[w] < tau)
								vis[color[w]] = false;
						}
						if (flag) break;
					}
				}
			}

		}
	}
	static Coloring GraphColoring(Graph &G, int tau) {
		Coloring c(G.n);
		c.graphColoring(G, tau);
		return c;	
	}
};

#endif // COLORING_HPP