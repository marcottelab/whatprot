// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_DYE_SEQS_IO_H
#define FLUOROSEQ_IO_DYE_SEQS_IO_H

#include <string>

#include "common/dye_seq.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void read_dye_seqs(const std::string& filename,
                   int* num_channels,
                   int* num_dye_seqs,
                   SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_DYE_SEQS_IO_H
