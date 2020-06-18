#include "osnove.h"

namespace hera {

	extern void parseProcessFastq(const std::string& strFastq, MapIdToSeq& mIdToSeq);
	extern void parseProcessFasta(const std::string& strFasta, MapIdToSeq& mIdToSeq);
	extern void parseProcessPaf(const std::string& strPaf, MapIdToOvlp& mIdToOvl);
	extern void parseProcessPaf(const std::string& strPaf, vOvlp& vOvlp);

}