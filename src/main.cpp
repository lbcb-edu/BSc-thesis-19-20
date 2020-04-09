#include <iostream>
#include <vector>

#include "../sequences/sequences.hpp"

#include <bioparser/bioparser.hpp>

int main(int argc, char* argv[]) {
    // std::cout << "In principio erat Verbum" << std::endl;

    std::vector<std::unique_ptr<Fast>> fasta_objects;
    auto fasta_parser = bioparser::createParser<bioparser::FastaParser, Fast>(argv[1]);
    fasta_parser->parse(fasta_objects, -1);

    std::cout << "Size: " << fasta_objects.size() << std::endl;
    for (int i = 0; i < fasta_objects.size(); i++) {
        std::unique_ptr<Fast>& tmp = fasta_objects[i];
        std::string s = tmp->name_;
        std::cout << s << std::endl;
        break;
    }

    return 0;
}