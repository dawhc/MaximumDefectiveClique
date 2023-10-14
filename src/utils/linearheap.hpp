#ifndef LINEAR_HEAP
#define LINEAR_HEAP

#include "vertexset.hpp"
#include <vector>
#include <algorithm>
#include <cassert>
#include <functional>

class LinearHeap {
	std::vector<int> rank, pos, bin;
	std::vector<int> keys, vals;
	int ptr;
public:

	LinearHeap(const std::vector<int>& values)
		: rank(values.size() + 1), pos(values.size() + 1), keys(values.size() + 1), vals(values), ptr(0) 
	{
		int maxValue = vals.empty() ? 0 : *std::max_element(vals.begin(), vals.end());
		bin.resize(maxValue + 2);
		for (int v : vals) ++bin[v+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int i = 0; i < vals.size(); ++i) {
			rank[i] = bin[vals[i]]++;
			pos[rank[i]] = i;
			keys[i] = i;
		}
	}

	LinearHeap(const std::vector<int>& keys, const std::vector<int>& values)
		: pos(keys.size() + 1), keys(keys), vals(values), ptr(0)
	{
		int maxValue = 0, maxKey = 0;
		for (int k : keys) {
			maxValue = std::max(maxValue, vals[k]);
			maxKey = std::max(maxKey, k);
		}
		bin.resize(maxValue + 2);
		rank.resize(maxKey + 1);
		for (int k : keys) ++bin[vals[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : keys) {
			rank[k] = bin[vals[k]]++;
			pos[rank[k]] = k;
		}
	}

	LinearHeap(const VertexSet& V, const std::vector<int>& values)
		: pos(V.size() + 1), keys(V.begin(), V.end()), vals(values), ptr(0)
	{
		int maxValue = 0, maxKey = keys.empty() ? 0 : *std::max_element(keys.begin(), keys.end());
		for (int k : keys) {
			maxValue = std::max(maxValue, vals[k]);
			maxKey = std::max(maxKey, k);
		}
		bin.resize(maxValue + 2);
		rank.resize(maxKey + 1);
		for (int k : V) ++bin[vals[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : V) {
			rank[k] = bin[vals[k]]++;
			pos[rank[k]] = k;
		}
	}

	LinearHeap(const VertexSet &V, const std::function<int(int)> &valueFunc)
	: pos(V.size() + 1), keys(V.begin(), V.end()), ptr(0)
	{
		int maxValue = 0, maxKey = keys.empty() ? 0 : *std::max_element(keys.begin(), keys.end());
		vals.resize(maxKey + 1);
		for (int k : keys) {
			vals[k] = valueFunc(k);
			maxValue = std::max(maxValue, vals[k]);
		}
		bin.resize(maxValue + 2);
		rank.resize(maxKey + 1);
		for (int k : V) ++bin[vals[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : V) {
			rank[k] = bin[vals[k]]++;
			pos[rank[k]] = k;
		}
	}

	void inc(int k) {
		//if (!inside(k)) return;
		assert(vals[k] < bin.size());
		int k2 = pos[--bin[vals[k]]];
		rank[k2] = rank[k];
		rank[k] = bin[vals[k]++];
		pos[rank[k]] = k;
		pos[rank[k2]] = k2;
	}
	
	void dec(int k) {
		//if (!inside(k)) return;
		assert(vals[k] > 0);
		if (vals[pos[ptr]] == vals[k])
			bin[vals[k]-1] = ptr;
		int k2 = pos[bin[--vals[k]]];
		rank[k2] = rank[k];
		rank[k] = bin[vals[k]]++;
		pos[rank[k]] = k;
		pos[rank[k2]] = k2;

	}

	int top() {
		return pos[ptr];
	}

	void pop() {
		++ptr;
	}

	bool empty() {
		return ptr >= keys.size();
	}

	bool inside(int k) {
		return rank[k] >= ptr;
	}

	int size() {
		return keys.size() - ptr;
	}

	int operator[] (int k) {
		return vals[k];
	}

};

#endif