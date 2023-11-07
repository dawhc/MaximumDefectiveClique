#ifndef HASH_HPP
#define HASH_HPP

#include <cassert>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <emmintrin.h>

constexpr int unfilled = -1;
constexpr int buff_size = sizeof(int);

class CuckooHash {
private:
	/* data */
	int capacity;
	int mask;
	int size;
	std::vector<int> hashtable;
	// int *hashtable;

	// void rehash(int **_table) {
	// 	int oldcapacity = capacity;
	// 	mask = mask == 0 ? 1 : ((mask << 1) | 1);
	// 	capacity = (mask + 1) * buff_size;
	// 	int *newhash = new int[capacity];
	// 	memset((newhash), unfilled, sizeof(int) * capacity);
	// 	for (int i = 0; i < oldcapacity; ++i){
	// 		if ((*_table)[i] != unfilled) insert((*_table)[i], &newhash);
	// 	}
	// 	std::swap((*_table), newhash);
	// 	delete[] newhash;
	// }
    
	// void insert(const int &_u, int **_table) {
		
	// 	int hs = hash1(_u);
	// 	for (int i = 0; i < buff_size; ++i) {
	// 		if ((*_table)[hs * buff_size + i] == unfilled){
	// 			(*_table)[hs * buff_size + i] = _u;
	// 			return;
	// 		}
	// 	}
	// 	hs = hash2(_u);
	// 	for (int i = 0; i < buff_size; ++i) {
	// 		if ((*_table)[hs * buff_size + i] == unfilled){
	// 			(*_table)[hs * buff_size + i] = _u;
	// 			return;
	// 		}
	// 	}

	// 	bool use_hash1 = true;
	// 	int u = _u;
	// 	for (int i = 0; i < mask; ++i) {
	// 		int replaced;
	// 		if (use_hash1) hs = hash1(u);
	// 		else hs = hash2(u);
	// 		int j = 0;
	// 		for (; j < buff_size; ++j) {
	// 			if ((*_table)[hs * buff_size + j] == unfilled) break;
	// 		}
	// 		if (buff_size == j) {
	// 			replaced = (*_table)[hs * buff_size];
	// 			j = 1;
	// 			for (; j < buff_size; j++) {
	// 				(*_table)[hs * buff_size + j - 1] = (*_table)[hs * buff_size + j];
	// 			}
	// 			(*_table)[hs * buff_size + j - 1] = u;
	// 		}
	// 		else {
	// 			replaced = (*_table)[hs * buff_size + j];
	// 			(*_table)[hs * buff_size + j] = u;
	// 		}
	// 		use_hash1 = hs == hash2(replaced);
	// 		u = replaced;
	// 		if (u == unfilled) return;
	// 	}
	// 	rehash(_table);
	// 	insert(u, _table);
	// }

	void rehash(std::vector<int> &table) {
		int oldcapacity = capacity;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		capacity = (mask + 1) * buff_size;
		//int *newhash = new int[capacity];
		//memset((newhash), unfilled, sizeof(int) * capacity);
		std::vector<int> newhash(capacity, unfilled);
		for (int i = 0; i < oldcapacity; ++i){
			if (table[i] != unfilled) insert(table[i], newhash);
		}
		//table = std::move(newhash);
		table.swap(newhash);
		//table = newhash;
		//std::swap((*_table), newhash);
		//delete[] newhash;
	}
    
	void insert(int u, std::vector<int>& table) {
		
		int hs = hash1(u);
		for (int i = 0; i < buff_size; ++i) {
			if (table[hs * buff_size + i] == unfilled){
				table[hs * buff_size + i] = u;
				return;
			}
		}
		hs = hash2(u);
		for (int i = 0; i < buff_size; ++i) {
			if (table[hs * buff_size + i] == unfilled){
				table[hs * buff_size + i] = u;
				return;
			}
		}

		bool use_hash1 = true;
		for (int i = 0; i < mask; ++i) {
			int replaced;
			if (use_hash1) hs = hash1(u);
			else hs = hash2(u);
			int j = 0;
			for (; j < buff_size; ++j) {
				if (table[hs * buff_size + j] == unfilled) break;
			}
			if (buff_size == j) {
				replaced = table[hs * buff_size];
				j = 1;
				for (; j < buff_size; j++) {
					table[hs * buff_size + j - 1] = table[hs * buff_size + j];
				}
				table[hs * buff_size + j - 1] = u;
			}
			else {
				replaced = table[hs * buff_size + j];
				table[hs * buff_size + j] = u;
			}
			use_hash1 = hs == hash2(replaced);
			u = replaced;
			if (u == unfilled) return;
		}
		rehash(table);
		insert(u, table);
	}

	int hash1(const int &x) const { return x & mask;}
	int hash2(const int &x) const { return ~x & mask;}

public:
	CuckooHash(/* args */) {
		clear();
	}
	~CuckooHash() {
		// if (hashtable) delete[] hashtable;
	}
	void clear() {
		capacity = mask = size = 0;
	}

	void reserve(int size) {
		if (capacity >= size) return;
		mask = mask == 0 ? 1 : ((mask << 1) | 1);
		while (size >= mask * buff_size) mask = (mask << 1) | 1;
		capacity = (mask + 1) * buff_size;
		std::vector<int>(capacity, unfilled).swap(hashtable);
		// if (hashtable) delete[] hashtable;
		// hashtable = new int[capacity];
		// memset(hashtable, unfilled, sizeof(int) * capacity);
	}

	void insert(int u) {
		if (size == capacity) rehash(hashtable);
		if (find(u)) return;
		insert(u, hashtable);
		size++;
	}

	bool find(int u) const {
		if (size == 0) return false;
		int hs1 = hash1(u);
		int hs2 = hash2(u);

		const int* hashtable_ptr = hashtable.data();
	
		assert(buff_size == 4 && sizeof(int) == 4);
		__m128i cmp = _mm_set1_epi32(u);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable_ptr[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable_ptr[buff_size * hs2]);
        __m128i flag = _mm_or_si128(_mm_cmpeq_epi32(cmp, b1), _mm_cmpeq_epi32(cmp, b2));

		return _mm_movemask_epi8(flag) != 0;
	}
	void erase(int u) {
		if (size == 0) return;
		int hs1 = hash1(u);
		int hs2 = hash2(u);
		int* hashtable_ptr = hashtable.data();
	
		assert(buff_size == 4 && sizeof(int) == 4);
		__m128i cmp = _mm_set1_epi32(u);
		__m128i b1 = _mm_load_si128((__m128i*)&hashtable_ptr[buff_size * hs1]);
		__m128i b2 = _mm_load_si128((__m128i*)&hashtable_ptr[buff_size * hs2]);
		__m128i flag1 = _mm_cmpeq_epi32(cmp, b1);
		__m128i flag2 = _mm_cmpeq_epi32(cmp, b2);
		if (_mm_movemask_epi8(flag1) != 0) {
			__m128i data = _mm_or_si128(b1, flag1);
			_mm_store_si128((__m128i*)&hashtable_ptr[buff_size * hs1], data);
		}
		else if (_mm_movemask_epi8(flag2) != 0) {
			__m128i data = _mm_or_si128(b2, flag2);
			_mm_store_si128((__m128i*)&hashtable_ptr[buff_size * hs2], data);
		}
	}
	int getcapacity() {return capacity;}
	int getsize() {return size;}
	int getmask() {return mask;}

	//int *gethashtable() {return hashtable;}

	bool operator[](const int &u) const {
		return find(u);
	}
};

#endif	