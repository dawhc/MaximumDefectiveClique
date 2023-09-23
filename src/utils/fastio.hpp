#ifndef FASTIO_HPP
#define FASTIO_HPP

#include <cstdio>
#include <string>

class FastIO {
	const int PAGESIZE = 1<<12;
	char *buf, *sp, *tp;
	FILE *_f;

	size_t nextPage() {
		size_t n = fread(buf, 1, PAGESIZE, _f);
		tp = (sp = buf) + n;
		return n;
	}
public:
	FastIO(FILE *f) {
		_f = f;
		buf = new char[PAGESIZE];
	}
	FastIO(const std::string& filename, const std::string& mode) {
		_f = fopen(filename.c_str(), mode.c_str());
		sp = tp = buf = new char[PAGESIZE];
	}
	~FastIO() {
		fclose(_f);
		delete[] buf;
	}
	bool empty() {
		// while (true) {
			// if (sp == tp && nextPage() == 0) return true;
			// for (; sp != tp && isspace(*sp); sp++);
			// if (!isspace(*sp)) break;
		// }

		if (sp == tp && nextPage() == 0) return true;
		return false;
	}
	char getChar() {
		
		if (empty()) return EOF;
		return *sp++;
	}
	unsigned int getUInt() {
		char c;
		unsigned int val = 0;
		for (c = getChar(); isspace(c); c = getChar());
		for (; !isspace(c); c = getChar())
			val = (val<<1) + (val<<3) + (c-'0');
		return val;
	}
	int getInt() {
		char c;
		int val = 0, sign = 1;
		for (c = getChar(); isspace(c); c = getChar());
		if (c == '-') sign = -1;
		for (; !isspace(c); c = getChar())
			val = (val<<1) + (val<<3) + (c-'0');
		return val * sign;
	}
};

#endif