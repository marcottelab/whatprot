// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "dedup_dye_tracks.h"

#include <unordered_map>
#include <utility>  // for std::move
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

#include "common/dye_track.h"
#include "common/sourced_data.h"
#include "keyvalue.h"
#include "mapreduce.h"

namespace fluoroseq {

namespace {
using MAPREDUCE_NS::KeyValue;
using MAPREDUCE_NS::MapReduce;
using std::move;
using std::unordered_map;
using std::vector;
}  // namespace

void dedup_dye_tracks(
        int num_timesteps,
        int num_channels,
        vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks_in,
        vector<SourcedData<DyeTrack,
                           SourceCountHitsList<int>>>* dye_tracks_out) {
#ifdef USE_MPI
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#else  // USE_MPI
    int mpi_size = 1;
#endif  // USE_MPI
    MapReduce mr(MPI_COMM_WORLD);
    // map function puts key/value pairs into MapReduce system.
    mr.map(mpi_size,  // nmap - this ensures one callback per process.
           &map_dye_tracks,  // callback function
           (void*) dye_tracks_in);  // ptr - passed through to callback.
    // collate function brings together key/value pairs with the same key.
    mr.collate(NULL);  // Collating with default hash function.
    // reduce converts keys with multiple values into keys with one value.
    mr.reduce(&reduce_dye_tracks,  // callback function
              NULL);  // not providing a pointer to be passed to callback.
    // gather brings the key/value pairs back to one process.
    mr.gather(1);  // parameter indicates nprocs to gather to
    // scan goes over all the key/value pairs and puts them into the
    // dye_tracks_out vector.
    OutputInfo output_info(num_timesteps, num_channels, dye_tracks_out);
    mr.scan(&output_dye_tracks,  // callback function
            (void*) &output_info);  // ptr - passed through to callback.
}

void map_dye_tracks(int index, KeyValue* kv, void* ptr) {
    vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks =
            (vector<SourcedData<DyeTrack, SourceCount<int>>>*) ptr;
    for (SourcedData<DyeTrack, SourceCount<int>>& dye_track : *dye_tracks) {
        kv->add((char*) &dye_track.value.counts[0],  // key
                dye_track.value.counts.size() * sizeof(short),  // keybytes
                (char*) &dye_track.source,  // value
                sizeof(SourceCount<int>));  // valuebytes
    }
}

void reduce_dye_tracks(char* key,
                       int keybytes,
                       char* multivalue,
                       int nvalues,
                       int* valuebytes,
                       KeyValue* kv,
                       void* ptr) {
    unordered_map<int, SourceCountHits<int>> source_map;
    if (multivalue == NULL) {
        MapReduce* mr = (MapReduce*) valuebytes;
        int nblocks;
        uint64_t nvalues_total = mr->multivalue_blocks(nblocks);
        source_map.reserve(nvalues_total);
        for (int i = 0; i < nblocks; i++) {
            nvalues = mr->multivalue_block(i, &multivalue, &valuebytes);
            reduce_dye_tracks_helper(multivalue, nvalues, &source_map);
        }
    } else {
        source_map.reserve(nvalues);
        reduce_dye_tracks_helper(multivalue, nvalues, &source_map);
    }
    int num_sources = source_map.size();
    SourceCountHits<int>* sources = new SourceCountHits<int>[num_sources];
    int index = 0;
    for (const auto& source_map_entry : source_map) {
        sources[index] = source_map_entry.second;
        index++;
    }
    kv->add(key,
            keybytes,
            (char*) sources,  // value
            num_sources * sizeof(SourceCountHits<int>));  // valuebytes
}

void reduce_dye_tracks_helper(
        char* multivalue,
        int nvalues,
        unordered_map<int, SourceCountHits<int>>* source_map) {
    SourceCount<int>* source_counts = (SourceCount<int>*) multivalue;
    for (int i = 0; i < nvalues; i++) {
        int source = source_counts[i].source;
        int count = source_counts[i].count;
        if (source_map->find(source) == source_map->end()) {
            (*source_map)[source] = SourceCountHits<int>(source, count, 1);
        } else {
            (*source_map)[source].hits++;
        }
    }
}

OutputInfo::OutputInfo(
        int num_timesteps,
        int num_channels,
        std::vector<SourcedData<
                DyeTrack, SourceCountHitsList<int>>>* dye_tracks_out)
        : num_timesteps(num_timesteps),
          num_channels(num_channels),
          dye_tracks_out(dye_tracks_out) {}
        

void output_dye_tracks(char* key,
                       int keybytes,
                       char* value,
                       int valuebytes,
                       void* ptr) {
    short* counts = (short*) key;
    SourceCountHits<int>* val_sources = (SourceCountHits<int>*) value;
    int num_sources = valuebytes / sizeof(SourceCountHits<int>);
    OutputInfo* output_info = (OutputInfo*) ptr;
    int num_timesteps = output_info->num_timesteps;
    int num_channels = output_info->num_channels;
    vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks =
            output_info->dye_tracks_out;
    SourceCountHits<int>** sources = new SourceCountHits<int>*[num_sources];
    for (int i = 0; i < num_sources; i++) {
        sources[i] = new SourceCountHits<int>(val_sources[i]);
    }
    dye_tracks->push_back(
            move(SourcedData<DyeTrack, SourceCountHitsList<int>>(
                    move(DyeTrack(num_timesteps, num_channels, counts)),
                    move(SourceCountHitsList<int>(num_sources, sources)))));
}

}  // namespace fluoroseq
