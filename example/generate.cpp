#include <iostream>
#include <fstream>
#include <vector>

/**
 * ./generate eschericiha_coli_reference.fasta
 * */

int main(int argc, char* argv[]) {
    std::ifstream f;
    f.open(argv[1]);
    std::ofstream kontinzi("kontizi.fasta");
    std::ofstream ocitanja("ocitanja.fastq");
    std::ofstream ukupno("ukupno.fasta");

    std::vector<std::string> v;
    std::string line;
    std::string baci, ime;
    int cnt = 0;
    while (std::getline(f, line)) {
        if (cnt++ == 0) {
            ime = line;
            continue;
        }
        baci += line;
        v.push_back(line);
    }
    f.close();

    kontinzi << ">kontig1";
    kontinzi << "\n";
    kontinzi << baci.substr(0, 20000);
    kontinzi << "\n";
    kontinzi << ">kontig2";
    kontinzi << "\n";
    kontinzi << baci.substr(25000, 20000);
    kontinzi << "\n";
    kontinzi << ">kontig3";
    kontinzi << "\n";
    kontinzi << baci.substr(50000, 20000);
    kontinzi.close();

    std::string tmp = baci.substr(0, 5000);
    ocitanja << "@ocitanje1";
    ocitanja << "\n";
    ocitanja << baci.substr(18000, 5000);
    ocitanja << "\n";
    ocitanja << "+\n";
    ocitanja << tmp;
    ocitanja << "\n";
    ocitanja << "@ocitanje2";
    ocitanja << "\n";
    ocitanja << baci.substr(21000, 5000);
    ocitanja << "\n";
    ocitanja << "+\n";
    ocitanja << tmp;
    ocitanja << "\n";
    ocitanja << "@ocitanje3";
    ocitanja << "\n";
    ocitanja << baci.substr(43000, 5000);
    ocitanja << "\n";
    ocitanja << "+\n";
    ocitanja << tmp;
    ocitanja << "\n";
    ocitanja << "@ocitanje4";
    ocitanja << "\n";
    ocitanja << baci.substr(46000, 5000);
    ocitanja << tmp;
    ocitanja << "\n";
    ocitanja.close();

    ukupno << baci.substr(0, 70000);
    ukupno.close();
    return 0;
}