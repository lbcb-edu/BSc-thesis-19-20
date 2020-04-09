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