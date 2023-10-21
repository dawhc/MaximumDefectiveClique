#ifndef VERTEX_SET
#define VERTEX_SET

#include <cstdint>
#include <vector>
#include <algorithm>
#include <cassert>
#include <string>
#include <sstream>

class VertexSet {
	int lp, rp, capacity;
	std::vector<int> s, pos;

  void swapByPos(int i, int j) {
		std::swap(s[i], s[j]);
		pos[s[i]] = i;
		pos[s[j]] = j;
	}

public:
	VertexSet(int n) {
		capacity = 0;
		reserve(n);
		clear();
	}

	VertexSet() {
		capacity = 0;
		clear();
	}
 
  void swapByVal(int u, int v) {
    swapByPos(pos[u], pos[v]);
  }

	void reserve(int new_capacity) {
		if (capacity >= new_capacity) return;
		s.resize(new_capacity);
		pos.resize(new_capacity);
		for (int i = capacity; i < new_capacity; ++i)
			s[i] = pos[i] = i;
		capacity = new_capacity;
	}

	void push(int v) {
		if (v >= capacity) reserve(v << 1);
		if (pos[v] < lp) pushFront(v);
		else if (pos[v] >= rp) pushBack(v);
	}

	void pushFront(int v) {
		if (pos[v] < lp) swapByPos(pos[v], --lp);
	}

	void pushBack(int v) {
		if (pos[v] >= rp) swapByPos(pos[v], rp++);
	}

	void pop(int v) {
		popFront(v);
	}

	void popFront(int v) {
		if (inside(v)) swapByPos(pos[v], lp++);
	}

	void popBack(int v) {
		if (inside(v)) swapByPos(pos[v], --rp);
	}


	bool inside(int v) const {
		return v < capacity && pos[v] >= lp && pos[v] < rp;
	}

  int size() const {
		return rp - lp;
	}

	int frontPos() const {
		return lp;
	}

	int backPos() const {
		return rp;
	}

	void restore(int pos) {
		lp = pos;
	}

	void clear() {
		lp = rp = 0;
	}

	int operator [] (int index) const {
		return s[index];
	}

	const int* begin() const {
		return s.data() + lp;
	}

	const int* end() const {
		return s.data() + rp;
	}

	friend VertexSet operator + (VertexSet a, VertexSet b) {
		VertexSet c = a.capacity > b.capacity ? std::move(a) : std::move(b);
		if (a.capacity > b.capacity) {
			for (int v : b)
				c.push(v);
		}
		else {
			for (int v : a)
				c.push(v);
		}
		return c;
	}

	std::string toString(const std::string &sep = ",") {
		std::stringstream ss;
		ss << "{";
		for (int i = lp; i < rp; ++i) {
			if (i > lp) ss << sep;
			ss << s[i];
		}
		ss << "}";
		return ss.str();
	}

};

#endif
