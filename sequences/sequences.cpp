#include "sequences.hpp"
#include <iostream>

Fast::Fast(const char* name, std::uint32_t name_length,
            const char* sequence, std::uint32_t sequence_length) :
            name_(name, name_length),
            sequence_(sequence, sequence_length) {};

Fast::Fast(const char* name, std::uint32_t name_length,
            const char* sequence, std::uint32_t sequence_length,
            const char* quality, std::uint32_t quality_length) :
            name_(name, name_length),
            sequence_(sequence, sequence_length),
            quality_(quality, quality_length) {};

Paf::Paf(const char* q_name, std::uint32_t q_name_length,
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
            std::uint32_t mapping_quality) :
            q_name_(q_name, q_name_length), 
            q_length_(q_length), q_begin_(q_begin_), q_end_(q_end),
            orientation_(orientation), t_name_(t_name, t_name_length),
            t_length_(t_length), t_begin_(t_begin), t_end_(t_end),
            matching_bases_(matching_bases), overlap_length_(overlap_length), 
            mapping_quality_(mapping_quality) {};


