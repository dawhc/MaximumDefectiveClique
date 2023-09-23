#ifndef BIPLEX_H
#define BIPLEX_H

#include "bigraph.hpp"
#include "vertexset.hpp"
#include "hash.hpp"
#include <cstdint>
#include <vector>

namespace biplex {

	struct Result {
		uint64_t numBiPlexes, numBranches, numPivotingBranches, numBipartiteBranches, time;
	};

	void coreOrdering(const BiGraph& G);

	BiGraph coreReduction(const BiGraph& G, uint32_t alpha, uint32_t beta);

	BiGraph butterflyReduction(const BiGraph& G, uint32_t threshold);

	bool pruneCX(uint32_t u, uint32_t uSide);

	bool constructSets(uint32_t u, uint32_t uSide);

	void destructSets(uint32_t u, uint32_t uSide);

	Result enumeration(const BiGraph& G, uint32_t q, uint32_t k, uint64_t resultNum);

	//uint32_t biplexLowerBound(const BiGraph& G);

	uint32_t biplexUpperBound(uint32_t u, uint32_t uSide);

	uint32_t calcDegree(VertexSet* V, uint32_t u, uint32_t uSide);

	void biplexPruneC(std::pair<uint32_t, uint32_t>& oldPosC);

	void biplexRestorePruneC(std::pair<uint32_t, uint32_t>& oldPosC);

	void biplexUpdate(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC, std::pair<uint32_t, uint32_t>& oldPosX);

	void biplexUpdateDeg(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC);

	void biplexUpdateDegC(std::pair<uint32_t, uint32_t>& oldPosC);

	void biplexUpdateCX(uint32_t u, uint32_t uSide, VertexSet* V, std::pair<uint32_t, uint32_t>& oldPos);

	void biplexRestore(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC, std::pair<uint32_t, uint32_t>& oldPosX);

	void biplexRestoreDeg(uint32_t u, uint32_t uSide, std::pair<uint32_t, uint32_t>& oldPosC);

	void biplexRestoreDegC(std::pair<uint32_t, uint32_t>& oldPosC);

	void biplexBacktrack(uint32_t u, uint32_t uSide);

	void branchSmallK(uint32_t dep, uint32_t begin, uint32_t u, uint32_t uSide);

	void branchNew(uint32_t dep);

	void printSet(VertexSet* V, const std::string& name);

	void printSets();

	void printResult(Result& result);

}

#endif
