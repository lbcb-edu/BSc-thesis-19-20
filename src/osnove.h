#pragma once

#include <map>
#include <vector>
#include <memory>

namespace hera {

	extern uint32_t MinMCPaths, HardNodeLimit, SoftNodeLimit;
	extern uint32_t NumDFSNodes, MaxMCIterations;
	extern uint32_t MinPathsinGroup;
	extern float SImin, OHmax;

	extern bool test_short_length, test_contained_reads, test_low_quality;

	extern bool print_output;

	extern double PathGroupHalfSize;

	extern std::string logFile;

    class Sequence;
    class Overlap;
    class Node;
    class Edge;

    using MapIdToOvlp = std::map<std::pair<std::string, int>, std::vector<std::shared_ptr<Overlap>>>;
    using MapIdToSeq = std::map<std::string, std::shared_ptr<Sequence>>;
    using MapIdToNode = std::map<std::string, std::shared_ptr<Node>>;
    using vOvlp = std::vector<std::unique_ptr<Overlap>>;


}