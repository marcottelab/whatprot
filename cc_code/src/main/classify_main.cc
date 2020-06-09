// Author: Matthew Beauregard Smith (UT Austin)

#include <cstring>
#include <iostream>

#include "main/ann_main.h"
#include "main/hmm_main.h"

namespace {
using fluoroseq::ann_main;
using fluoroseq::hmm_main;
using std::cout;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "bad inputs\n";
        return 1;
    }
    char* mode = argv[1];
    if (0 == strcmp(mode, "hmm")) {
        return hmm_main(argc, argv);
    } else if (0 == strcmp(mode, "ann")) {
        return ann_main(argc, argv);
    }
    return 0;
}