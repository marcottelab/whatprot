// Author: Matthew Beauregard Smith (UT Austin)
#include "generate_radiometries.h"

#include <random>
#include <vector>

#include "common/dye_seq.h"
#include "common/radiometry.h"
#include "common/error_model.h"
#include "simulation/generate_radiometry.h"

namespace fluoroseq {

void generate_radiometries(
        const ErrorModel& error_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int radiometries_per_dye_seq,
        std::default_random_engine* generator,
        std::vector<SourcedData<Radiometry, SourceCount<int>>>* radiometries) {
    radiometries->reserve(dye_seqs.size() * radiometries_per_dye_seq);
    for (const SourcedData<DyeSeq, SourceCount<int>>& dye_seq : dye_seqs) {
        for (int i = 0; i < radiometries_per_dye_seq; i++) {
            radiometries->push_back(
                    SourcedData<Radiometry, SourceCount<int>>(
                            Radiometry(num_timesteps, num_channels),
                            dye_seq.source));
            generate_radiometry(error_model,
                                dye_seq.value,
                                num_timesteps,
                                num_channels,
                                generator,
                                &radiometries->back().value);
        }
    }
}

}  // namespace fluoroseq
