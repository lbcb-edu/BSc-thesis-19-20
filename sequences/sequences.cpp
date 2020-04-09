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
