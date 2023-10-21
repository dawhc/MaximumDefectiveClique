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
	int ptr, capKey, capVal;

	void resizeKey(int capacity) {
		capKey = capacity;
		pos.resize(capacity + 1);
		rank.resize(capacity);
		// keys.resize(capacity);
		vals.resize(capacity);
	}
	void resizeVal(int capacity) {
		capVal = capacity;
		bin.resize(capacity + 2);
	}
public:
	LinearHeap(): ptr(0), capKey(0), capVal(0) {}
	LinearHeap(const std::vector<int> &values): capKey(0), capVal(0) {
		build(values);
	}

	void build(const std::vector<int> &values) {
		if (values.size() > capKey) resizeKey(values.size());
		ptr = 0;
		int maxValue = 0;
		for (int i = 0; i < values.size(); ++i) { 
			keys[i] = i;
			vals[i] = values[i];
			maxValue = std::max(maxValue, vals[i]); 
		}
		// if (capVal < maxValue) resizeVal(maxValue);
		bin.assign(maxValue + 2, 0);
		for (int v : values) ++bin[v+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int i = 0; i < vals.size(); ++i) {
			rank[i] = bin[values[i]]++;
			pos[rank[i]] = i;
		}
	}

	LinearHeap(const std::vector<int> &keys, const std::vector<int> &values) {
		build(keys, values);
	}
	void build(const std::vector<int> &keys, const std::vector<int> &values) {
		this->keys.assign(keys.begin(), keys.end());
		ptr = 0;
		int maxValue = 0, maxKey = 0;
		for (int k : keys) {
			maxValue = std::max(maxValue, vals[k]);
			maxKey = std::max(maxKey, k);
		}
		if (maxKey >= capKey) resizeKey(maxKey+1);
		bin.assign(maxValue + 2, 0);
		for (int k : keys) ++bin[values[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : keys) {
			vals[k] = values[k];
			rank[k] = bin[values[k]]++;
			pos[rank[k]] = k;
		}
	}

	LinearHeap(const VertexSet &V, const std::vector<int> &values) {
		build(V, values);
	}

	void build(const VertexSet &V, const std::vector<int> &values) {
		keys.assign(V.begin(), V.end());
		ptr = 0;
		int maxValue = 0, maxKey = 0;
		for (int k : keys) {
			maxValue = std::max(maxValue, values[k]);
			maxKey = std::max(maxKey, k);
		}
		bin.assign(maxValue + 2, 0);
		if (maxKey >= capKey) resizeKey(maxKey+1);
		for (int k : V) ++bin[values[k]+1];
		for (int v = 1; v <= maxValue; ++v) bin[v] += bin[v-1];
		for (int k : V) {
			vals[k] = values[k];
			rank[k] = bin[values[k]]++;
			pos[rank[k]] = k;
		}
	}

	LinearHeap(const VertexSet &V, const std::function<int(int)> &valueFunc) {
		build(V, valueFunc);
	}

	void build(const VertexSet &V, const std::function<int(int)> &valueFunc) {
		keys.assign(V.begin(), V.end());
		ptr = 0;
		int maxValue = 0, maxKey = 0;
		for (int k : V) { maxKey = std::max(maxKey, k); }
		if (maxKey >= capKey) resizeKey(maxKey+1);
		for (int k : keys) {
			vals[k] = valueFunc(k);
			maxValue = std::max(maxValue, vals[k]);
		}
		// if (maxValue > capVal) resizeVal(maxValue+1);
		bin.assign(maxValue + 2, 0);
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