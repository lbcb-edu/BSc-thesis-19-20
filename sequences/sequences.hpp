#include <cstdint>
#include <iostream>

class Fast {
    public:
        std::string name_;
        std::string sequence_;
        std::string quality_;

        Fast(const char* name, std::uint32_t name_length,
            const char* sequence, std::uint32_t sequence_length);
        
        Fast(const char* name, std::uint32_t name_length,
            const char* sequence, std::uint32_t sequence_length,
            const char* quality, std::uint32_t quality_length);
};

class Paf {
    public:
        std::string q_name_;
        std::uint32_t q_name_length_;
        std::uint32_t q_length_;
        std::uint32_t q_begin_;
        std::uint32_t q_end_;
        char orientation_;
        std::string t_name_;
        std::uint32_t t_name_length_;
        std::uint32_t t_length_;
        std::uint32_t t_begin_;
        std::uint32_t t_end_;
        std::uint32_t matching_bases_;
        std::uint32_t overlap_length_;
        std::uint32_t mapping_quality_;

        Paf(const char* q_name, std::uint32_t q_name_length,
            std::uint32_t q_length,
            std::uint32_t q_begin,
            std::uint32_t q_end,
            char orientation,
            const char* t_name, std::uint32_t t_name_length,
            std::uint32_t t_length,
            std::uint32_t t_begin,
            std::uint32_t t_end,
            std::uint32_t matching_bases,
            std::uint32_t overlap_length,
            std::uint32_t mapping_quality);
};