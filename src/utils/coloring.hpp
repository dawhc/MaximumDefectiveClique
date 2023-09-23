#ifndef COLORING_HPP
#define COLORING_HPP

#pragma once

#include "graph.hpp"
#include "ordering.hpp"
#include <cstdint>
#include <vector>

constexpr int uncolored = -1;

class Coloring {
public:
	std::vector<int> color;
	int numColors;

	Coloring(Graph& G): color(G.n, uncolored), numColors(0) {}

	Coloring(): numColors(0) {}

	void colorVertex(int v, int c) {
		color[v] = c;
		numColors = std::max(numColors, c+1);
	}

	static Coloring graphColoring(Graph& G, int tau) {
		if (tau < 0) tau = 0;
		Coloring c(G);
		
		std::vector<bool> vis(G.n, false);
		std::vector<std::vector<int>> bin;

		Ordering o = Ordering::degeneracyOrdering(G);
		for (int i = o.numOrdered - 1; i >= 0; --i) {
			int u = o.ordered[i];

			// Coloring
			for (int v : G.nbr[u]) {
				if (c.color[v] != uncolored)
					vis[c.color[v]] = true;
			}
			for (int i = 0; i <= G.degree[u]; ++i) {
				if (!vis[i]) {
					c.colorVertex(u, i);
					break;
				}
			}
			for (int v : G.nbr[u]) {
				if (c.color[v] != uncolored)
					vis[c.color[v]] = false;
			}

			// Recoloring
			if (c.color[u] >= tau) {
				std::vector<std::vector<int>>(tau).swap(bin);
				for (int v : G.nbr[u]) {
					if (c.color[v] != uncolored && c.color[v] < tau)
						bin[c.color[v]].push_back(v);
				}
				for (int i = 0; i < tau; ++i) {
					if (bin[i].size() == 1) {
						int v = bin[i][0];
						bool flag = false;
						for (int w : G.nbr[v]) {
							if (c.color[w] != uncolored && c.color[w] < tau)
								vis[c.color[w]] = true;
						}
						for (int j = 0; j < tau; ++j) {
							if (!vis[j] && j != i) {
								flag = true;
								c.colorVertex(v, j);
								c.colorVertex(u, i);
								break;
							}
						}
						for (int w : G.nbr[v]) {
							if (c.color[w] != uncolored && c.color[w] < tau)
								vis[c.color[w]] = false;
						}
						if (flag) break;
					}
				}
			}

		}

		return c;
	}
};

#endif // COLORING_HPP