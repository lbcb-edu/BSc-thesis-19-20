#include <iostream>
#include <vector>
#include <fstream>

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

    std::vector<std::unique_ptr<Paf>> paf_objects;
    auto paf_parser = bioparser::createParser<bioparser::PafParser, Paf>(argv[2]);
    paf_parser->parse(paf_objects, -1);

    std::cout << "paf: " << paf_objects.size() << std::endl;

    for (int i = 0; i < paf_objects.size(); i++) {
        std::unique_ptr<Paf>& tmp = paf_objects[i];
        std::cout << tmp->q_name_ << " -> " << tmp->t_name_ << " : " << tmp->overlap_length_ << std::endl;
    }


    return 0;
}
