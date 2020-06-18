#include <unordered_set>

#include <bioparser/bioparser.hpp>

#include "parsiranje.h"
#include "../sequences/sequences.h"
#include "preklapanje.h"



namespace hera {


    void parseProcessFastq(const std::string& strFastq, MapIdToSeq& mIdToSeq) {
        std::vector<std::unique_ptr<Sequence>> aReads;
        auto fastqParser = bioparser::createParser<bioparser::FastqParser, hera::Sequence>(strFastq);
        fastqParser->parse(aReads, -1);
        for (uint32_t i = 0; i < aReads.size(); i++) {
            std::unique_ptr<Sequence> newSeq_ptr = std::make_unique<Sequence>(aReads[i]->seq_strName.c_str(), (uint32_t)(aReads[i]->seq_strName.length()),
                                                                aReads[i]->seq_strData.c_str(), (uint32_t)(aReads[i]->seq_strData.length()));

            mIdToSeq.emplace(newSeq_ptr->seq_strName, std::move(newSeq_ptr));
        }
    }

    void parseProcessFasta(const std::string& strFasta, MapIdToSeq& mIdToSeq) {
        std::vector<std::unique_ptr<Sequence>> aReads;
        auto fastaParser = bioparser::createParser<bioparser::FastaParser, hera::Sequence>(strFasta);
        fastaParser->parse(aReads, -1);
        for (uint32_t i = 0; i < aReads.size(); i++) {
            mIdToSeq.emplace(aReads[i]->seq_strName, std::move(aReads[i]));
        }
    }

    void parseProcessPaf(const std::string& strPaf, MapIdToOvlp& mIdToOvl) {
        std::vector<std::unique_ptr<hera::Overlap>> aExtensions;
        auto pafParser = bioparser::createParser<bioparser::PafParser, Overlap>(strPaf);
        pafParser->parse(aExtensions, -1);

        // preprocessing - mark contained reads
        std::unordered_set<std::string> usContained;
        for (const auto& ext : aExtensions) {
            if (ext->ext_oType == ET_CONTAINED) {
                usContained.emplace(ext->ext_strName);
            }
        }

        for (uint32_t i = 0; i < aExtensions.size(); i++) {
            if (aExtensions[i]->ext_oType != ET_INVALID && usContained.find(aExtensions[i]->ext_strName) == usContained.end()) {
                mIdToOvl[{aExtensions[i]->ext_strTarget, aExtensions[i]->ext_oType}].emplace_back(std::move(aExtensions[i]));
            }
        }
    }

    extern void parseProcessPaf(const std::string& strPaf, vOvlp& vOvlp) {
        std::vector<std::unique_ptr<hera::Overlap>> aExtensions;
        auto pafParser = bioparser::createParser<bioparser::PafParser, Overlap>(strPaf);
        pafParser->parse(aExtensions, -1);

        for (uint32_t i = 0; i < aExtensions.size(); i++) {
            vOvlp.emplace_back(std::move(aExtensions[i]));
        }
    }
}