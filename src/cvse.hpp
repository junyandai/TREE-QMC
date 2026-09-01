#pragma once

#include "tree.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#if ENABLE_TOB

namespace cvse {

enum class CVQuartetSamplingMode {
    ThreeFixOneAlter,
    TwoFixTwoAlter
};

enum class CVSEWeightingMode {
    HardTopK,
    Soft
};

struct CVSEConfig {
    std::size_t outer_folds = 5;
    std::size_t inner_folds = 4;
    std::size_t repeats = 1;
    std::size_t max_quartets_per_edge = 0; // 0 means no cap
    double pseudocount = 0.5;
    unsigned seed = 42;


    CVSEWeightingMode weighting_mode = CVSEWeightingMode::HardTopK;

    // Coverage filter, applied separately inside each training split.
    std::size_t min_observed_loci = 0;
    double min_coverage_fraction = 0.0;
};

struct CVSECandidate {
    std::array<index_t, 4> taxa;
    int expected_topology = -1;
};

struct CVSEResult {
    bool non_tree_like = false;
    double mean_outer_gain = 0.0;
    double se_outer_gain = 0.0;

    //   HardTopK -> median selected k
    //   Soft     -> median number of candidates with nonzero weight
    // across successful outer fits.
    std::size_t median_selected_k = 0;
    std::size_t outer_fit_count = 0;

    // Soft-mode tuning diagnostic. Remains 1.0 in hard mode.
    double mean_selected_alpha = 1.0;

    std::vector<double> selection_frequency;
    std::size_t representative_candidate =
        std::numeric_limits<std::size_t>::max();

    // Distinct selected candidates classified on the final full-data fit:
    // [wrong topology 1, wrong topology 2, unrestricted/asymmetric].
    std::array<std::size_t, 3> selected_alternative_types{{0, 0, 0}};

    // Full-data diagnostic counts only; CV eligibility remains split-specific.
    std::size_t eligible_candidate_count = 0;
    std::size_t filtered_candidate_count = 0;
};

CVQuartetSamplingMode parse_cv_quartet_sampling_mode(
    const std::string &mode_string
);

CVSEWeightingMode parse_cvse_weighting_mode(
    const std::string &mode_string
);

std::string cvse_weighting_mode_name(
    CVSEWeightingMode mode
);

std::vector<CVSECandidate> generate_3f1a_cvse_candidates(
    const std::array<std::vector<index_t>, 4>& partitions,
    std::size_t maximum,
    unsigned seed
);

std::vector<CVSECandidate> generate_2f2a_cvse_candidates(
    const std::vector<index_t>& left,
    const std::vector<index_t>& right,
    std::size_t maximum,
    unsigned seed
);

CVSEResult run_cvse_for_edge(
    std::vector<Tree *> &gene_trees,
    const std::vector<CVSECandidate> &candidates,
    const CVSEConfig &config
);

std::array<weight_t, 3> full_cvse_counts_for_candidate(
    std::vector<Tree *> &gene_trees,
    const CVSECandidate &candidate
);

} // namespace cvse

#endif // ENABLE_TOB