
#include "cvse.hpp"

#if ENABLE_TOB

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cvse {

namespace {

struct CVSEModels {
    std::array<double, 3> tree_probability{{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0
    }};
    std::array<double, 3> unrestricted_probability{{
        1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0
    }};
    double training_gain = 0.0;          // mean gain per observed training locus
    std::size_t observed_loci = 0;       // usable loci for this candidate in this fit
    bool eligible = false;               // determined from this fit's training loci only
};

} // namespace

static std::string cvse_lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

CVQuartetSamplingMode parse_cv_quartet_sampling_mode(const std::string &mode_string) {
    const std::string mode = cvse_lowercase(mode_string);

    if (mode == "3f1a" || mode == "3f1a-hard" || mode == "3f1a-soft" ||
        mode == "3-fix-1-alter" || mode == "three_fix_one_alter" ||
        mode == "three-fix-one-alter") {
        return CVQuartetSamplingMode::ThreeFixOneAlter;
    }

    if (mode == "2f2a" || mode == "2f2a-hard" || mode == "2f2a-soft" ||
        mode == "2-fix-2-alter" || mode == "two_fix_two_alter" ||
        mode == "two-fix-two-alter") {
        return CVQuartetSamplingMode::TwoFixTwoAlter;
    }

    throw std::invalid_argument(
        "Unknown CV quartet sampling mode '" + mode_string +
        "'. Expected 3f1a[-hard|-soft] or 2f2a[-hard|-soft]."
    );
}

CVSEWeightingMode parse_cvse_weighting_mode(const std::string &mode_string) {
    const std::string mode = cvse_lowercase(mode_string);


    if (mode == "hard" || mode == "topk" || mode == "hard-top-k" ||
        mode == "3f1a" || mode == "2f2a" ||
        mode == "3f1a-hard" || mode == "2f2a-hard" ||
        mode == "3-fix-1-alter" || mode == "three_fix_one_alter" ||
        mode == "three-fix-one-alter" || mode == "2-fix-2-alter" ||
        mode == "two_fix_two_alter" || mode == "two-fix-two-alter") {
        return CVSEWeightingMode::HardTopK;
    }

    if (mode == "soft" || mode == "3f1a-soft" || mode == "2f2a-soft") {
        return CVSEWeightingMode::Soft;
    }

    throw std::invalid_argument(
        "Unknown CVSE weighting mode '" + mode_string +
        "'. Expected hard/topk, 3f1a[-hard|-soft], or 2f2a[-hard|-soft]."
    );
}

std::string cvse_weighting_mode_name(CVSEWeightingMode mode) {
    switch (mode) {
        case CVSEWeightingMode::HardTopK:
            return "hard-top-k";
        case CVSEWeightingMode::Soft:
            return "soft";
    }
    return "unknown";
}

static int expected_topology_after_sorting(
    const std::array<index_t, 4> &group_order_taxa) {

    std::array<index_t, 4> sorted_taxa = group_order_taxa;
    std::sort(sorted_taxa.begin(), sorted_taxa.end());

    int first_position = -1;
    int second_position = -1;

    for (int i = 0; i < 4; ++i) {
        if (sorted_taxa[i] == group_order_taxa[0]) first_position = i;
        if (sorted_taxa[i] == group_order_taxa[1]) second_position = i;
    }

    if (first_position < 0 || second_position < 0 ||
        first_position == second_position) {
        return -1;
    }

    if (first_position > second_position) {
        std::swap(first_position, second_position);
    }

    if ((first_position == 0 && second_position == 1) ||
        (first_position == 2 && second_position == 3)) {
        return 0;
    }
    if ((first_position == 0 && second_position == 2) ||
        (first_position == 1 && second_position == 3)) {
        return 1;
    }
    return 2;
}

static bool make_cvse_candidate(const std::array<index_t, 4> &group_order_taxa,
                         CVSECandidate &candidate) {
    std::array<index_t, 4> sorted_taxa = group_order_taxa;
    std::sort(sorted_taxa.begin(), sorted_taxa.end());

    if (std::adjacent_find(sorted_taxa.begin(), sorted_taxa.end()) !=
        sorted_taxa.end()) {
        return false;
    }

    const int expected_topology =
        expected_topology_after_sorting(group_order_taxa);
    if (expected_topology < 0) return false;

    candidate.taxa = sorted_taxa;
    candidate.expected_topology = expected_topology;
    return true;
}

static void cap_cvse_candidates(std::vector<CVSECandidate> &candidates,
                         std::size_t maximum,
                         unsigned seed) {
    if (maximum == 0 || candidates.size() <= maximum) return;

    std::mt19937 rng(seed);
    std::shuffle(candidates.begin(), candidates.end(), rng);
    candidates.resize(maximum);

    std::sort(candidates.begin(), candidates.end(),
              [](const CVSECandidate &left, const CVSECandidate &right) {
                  if (left.taxa != right.taxa) return left.taxa < right.taxa;
                  return left.expected_topology < right.expected_topology;
              });
}

std::vector<CVSECandidate> generate_3f1a_cvse_candidates(
    const std::array<std::vector<index_t>, 4>& partitions,
    std::size_t maximum,
    unsigned seed) {

    for (const auto& partition : partitions) {
        if (partition.empty()) return {};
    }


    std::array<index_t, 4> anchors = {
        partitions[0][0],
        partitions[1][0],
        partitions[2][0],
        partitions[3][0]
    };

    std::set<std::array<index_t, 4>> seen;
    std::vector<CVSECandidate> candidates;
    candidates.reserve(
        partitions[0].size() + partitions[1].size() +
        partitions[2].size() + partitions[3].size()
    );

    for (std::size_t alter = 0; alter < 4; ++alter) {
        for (index_t taxon : partitions[alter]) {
            std::array<index_t, 4> group_order_taxa = anchors;
            group_order_taxa[alter] = taxon;

            CVSECandidate candidate;
            if (!make_cvse_candidate(group_order_taxa, candidate)) continue;

            if (seen.insert(candidate.taxa).second) {
                candidates.push_back(candidate);
            }
        }
    }

    cap_cvse_candidates(candidates, maximum, seed);
    return candidates;
}

std::vector<CVSECandidate> generate_2f2a_cvse_candidates(
    const std::vector<index_t>& left,
    const std::vector<index_t>& right,
    std::size_t maximum,
    unsigned seed) {

    if (left.size() < 2 || right.size() < 2) return {};

    const index_t left_anchor = left[0];
    const index_t right_anchor = right[0];

    std::vector<CVSECandidate> candidates;
    candidates.reserve((left.size() - 1) * (right.size() - 1));

    for (std::size_t i = 1; i < left.size(); ++i) {
        for (std::size_t j = 1; j < right.size(); ++j) {
            const std::array<index_t, 4> group_order_taxa = {
                left_anchor, left[i],
                right_anchor, right[j]
            };

            CVSECandidate candidate;
            if (make_cvse_candidate(group_order_taxa, candidate)) {
                candidates.push_back(candidate);
            }
        }
    }

    cap_cvse_candidates(candidates, maximum, seed);
    return candidates;
}

static std::vector<std::vector<std::size_t>> make_cvse_folds(
    const std::vector<std::size_t> &loci,
    std::size_t requested_folds,
    unsigned seed) {

    if (loci.empty()) return {};

    const std::size_t fold_count =
        std::max<std::size_t>(1, std::min(requested_folds, loci.size()));

    std::vector<std::size_t> shuffled = loci;
    std::mt19937 rng(seed);
    std::shuffle(shuffled.begin(), shuffled.end(), rng);

    std::vector<std::vector<std::size_t>> folds(fold_count);
    for (std::size_t i = 0; i < shuffled.size(); ++i) {
        folds[i % fold_count].push_back(shuffled[i]);
    }
    return folds;
}

static std::vector<std::size_t> concatenate_other_cvse_folds(
    const std::vector<std::vector<std::size_t>> &folds,
    std::size_t held_out) {

    std::vector<std::size_t> result;
    for (std::size_t i = 0; i < folds.size(); ++i) {
        if (i == held_out) continue;
        result.insert(result.end(), folds[i].begin(), folds[i].end());
    }
    return result;
}

static std::array<double, 3> count_cvse_outcomes(
    const std::vector<std::int8_t> &outcomes,
    std::size_t number_of_loci,
    std::size_t candidate_index,
    const std::vector<std::size_t> &loci) {

    std::array<double, 3> counts{{0.0, 0.0, 0.0}};
    const std::size_t offset = candidate_index * number_of_loci;

    for (std::size_t locus : loci) {
        const int topology = static_cast<int>(outcomes[offset + locus]);
        if (topology >= 0 && topology < 3) counts[topology] += 1.0;
    }
    return counts;
}

static double cvse_log_likelihood(const std::array<double, 3> &counts,
                           const std::array<double, 3> &probability) {
    double result = 0.0;
    for (int topology = 0; topology < 3; ++topology) {
        if (counts[topology] <= 0.0) continue;
        result += counts[topology] *
                  std::log(std::max(probability[topology], 1.0e-300));
    }
    return result;
}

static std::array<double, 3> fit_cvse_tree_probability(
    const std::array<double, 3> &counts,
    int expected_topology,
    double pseudocount) {

    const double total = counts[0] + counts[1] + counts[2];
    if (expected_topology < 0 || expected_topology > 2) {
        return {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
    }

    const double smoothed_total = total + 3.0 * pseudocount;
    if (smoothed_total <= 0.0) {
        return {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
    }

    const double discordant =
        total - counts[expected_topology] + 2.0 * pseudocount;
    double x = discordant / (2.0 * smoothed_total);
    x = std::max(1.0e-12, std::min(1.0 / 3.0, x));

    std::array<double, 3> probability{{x, x, x}};
    probability[expected_topology] = 1.0 - 2.0 * x;
    return probability;
}

static std::array<double, 3> fit_cvse_unrestricted_probability(
    const std::array<double, 3> &counts,
    double pseudocount) {

    const double total = counts[0] + counts[1] + counts[2];
    const double denominator = total + 3.0 * pseudocount;
    if (denominator <= 0.0) {
        return {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
    }

    return {{
        (counts[0] + pseudocount) / denominator,
        (counts[1] + pseudocount) / denominator,
        (counts[2] + pseudocount) / denominator
    }};
}

static std::size_t cvse_minimum_required_observations(
    std::size_t training_locus_count,
    const CVSEConfig &config) {

    const double coverage_fraction =
        std::max(0.0, std::min(1.0, config.min_coverage_fraction));

    const std::size_t minimum_by_fraction =
        static_cast<std::size_t>(
            std::ceil(
                coverage_fraction * static_cast<double>(training_locus_count) -
                1.0e-12
            )
        );

    return std::max(config.min_observed_loci, minimum_by_fraction);
}

static std::vector<CVSEModels> fit_cvse_models(
    const std::vector<std::int8_t> &outcomes,
    std::size_t number_of_loci,
    const std::vector<CVSECandidate> &candidates,
    const std::vector<std::size_t> &training_loci,
    const CVSEConfig &config) {

    std::vector<CVSEModels> models(candidates.size());

    const std::size_t minimum_required =
        cvse_minimum_required_observations(training_loci.size(), config);

    for (std::size_t candidate = 0; candidate < candidates.size(); ++candidate) {
        const auto counts = count_cvse_outcomes(
            outcomes, number_of_loci, candidate, training_loci
        );

        const double observed = counts[0] + counts[1] + counts[2];
        models[candidate].observed_loci =
            static_cast<std::size_t>(observed);

        // Even when both configured thresholds are disabled, a model still
        // needs at least one observed training quartet to be estimable.
        models[candidate].eligible =
            models[candidate].observed_loci > 0 &&
            models[candidate].observed_loci >= minimum_required;

        models[candidate].tree_probability = fit_cvse_tree_probability(
            counts, candidates[candidate].expected_topology, config.pseudocount
        );
        models[candidate].unrestricted_probability =
            fit_cvse_unrestricted_probability(counts, config.pseudocount);

        const double raw_gain =
            cvse_log_likelihood(
                counts,
                models[candidate].unrestricted_probability
            ) -
            cvse_log_likelihood(
                counts,
                models[candidate].tree_probability
            );

        /*
         * Among candidates that pass the training-split coverage filter, rank
         * by mean gain per usable training locus. This prevents a candidate
         * from ranking higher merely because it has greater taxon coverage.
         */
        models[candidate].training_gain =
            models[candidate].eligible && observed > 0.0
                ? raw_gain / observed
                : -std::numeric_limits<double>::infinity();
    }

    return models;
}


static std::pair<double, double> cvse_mean_and_se(
    const std::vector<double> &values) {

    if (values.empty()) return {0.0, 0.0};

    const double mean =
        std::accumulate(values.begin(), values.end(), 0.0) /
        static_cast<double>(values.size());

    if (values.size() < 2) return {mean, 0.0};

    double sum_squared_difference = 0.0;
    for (double value : values) {
        const double difference = value - mean;
        sum_squared_difference += difference * difference;
    }

    const double variance =
        sum_squared_difference / static_cast<double>(values.size() - 1);
    const double standard_error =
        std::sqrt(variance / static_cast<double>(values.size()));

    return {mean, standard_error};
}



static std::vector<std::size_t> rank_cvse_candidates(
    const std::vector<CVSEModels> &models) {

    std::vector<std::size_t> ranking;
    ranking.reserve(models.size());

    for (std::size_t candidate = 0; candidate < models.size(); ++candidate) {
        if (models[candidate].eligible) {
            ranking.push_back(candidate);
        }
    }

    std::stable_sort(
        ranking.begin(), ranking.end(),
        [&](std::size_t left, std::size_t right) {
            return models[left].training_gain > models[right].training_gain;
        }
    );
    return ranking;
}

static std::size_t choose_cvse_exception_count(
    const std::vector<std::int8_t> &outcomes,
    std::size_t number_of_loci,
    const std::vector<CVSECandidate> &candidates,
    const std::vector<std::size_t> &outer_training_loci,
    const CVSEConfig &config,
    unsigned seed) {

    if (candidates.empty() || outer_training_loci.size() < 4) return 0;

    const std::size_t inner_fold_count =
        std::min(config.inner_folds, outer_training_loci.size());
    if (inner_fold_count < 2) return 0;

    const auto inner_folds = make_cvse_folds(
        outer_training_loci, inner_fold_count, seed
    );

    const std::size_t candidate_count = candidates.size();

    struct InnerFit {
        std::vector<std::size_t> validation_loci;
        std::vector<CVSEModels> models;
        std::vector<std::size_t> ranking;
    };

    std::vector<InnerFit> inner_fits;
    inner_fits.reserve(inner_folds.size());

    std::size_t common_max_k = candidate_count;

    for (std::size_t fold = 0; fold < inner_folds.size(); ++fold) {
        const auto inner_training_loci =
            concatenate_other_cvse_folds(inner_folds, fold);

        InnerFit fit;
        fit.validation_loci = inner_folds[fold];
        fit.models = fit_cvse_models(
            outcomes,
            number_of_loci,
            candidates,
            inner_training_loci,
            config
        );
        fit.ranking = rank_cvse_candidates(fit.models);

        common_max_k = std::min(common_max_k, fit.ranking.size());
        inner_fits.push_back(std::move(fit));
    }

    if (common_max_k == 0) return 0;

    std::vector<std::vector<double>> validation_gains(common_max_k + 1);

    for (const InnerFit &fit : inner_fits) {
        for (std::size_t locus : fit.validation_loci) {
            std::vector<double> candidate_gains(candidate_count, 0.0);
            std::size_t eligible_observed = 0;

            for (std::size_t candidate = 0;
                 candidate < candidate_count; ++candidate) {
                if (!fit.models[candidate].eligible) continue;

                const int topology = static_cast<int>(
                    outcomes[candidate * number_of_loci + locus]
                );
                if (topology < 0 || topology > 2) continue;

                ++eligible_observed;

                const double tree_log_probability = std::log(std::max(
                    fit.models[candidate].tree_probability[topology],
                    1.0e-300
                ));
                const double unrestricted_log_probability = std::log(std::max(
                    fit.models[candidate].unrestricted_probability[topology],
                    1.0e-300
                ));

                candidate_gains[candidate] =
                    unrestricted_log_probability - tree_log_probability;
            }

            if (eligible_observed == 0) continue;

            validation_gains[0].push_back(0.0);
            double cumulative_gain = 0.0;

            for (std::size_t k = 1; k <= common_max_k; ++k) {
                const std::size_t candidate = fit.ranking[k - 1];
                cumulative_gain += candidate_gains[candidate];

                validation_gains[k].push_back(
                    cumulative_gain /
                    static_cast<double>(eligible_observed)
                );
            }
        }
    }

    std::vector<double> mean_gain(
        common_max_k + 1,
        -std::numeric_limits<double>::infinity()
    );
    std::vector<double> se_gain(common_max_k + 1, 0.0);

    for (std::size_t k = 0; k <= common_max_k; ++k) {
        if (validation_gains[k].empty()) continue;
        const auto mean_se = cvse_mean_and_se(validation_gains[k]);
        mean_gain[k] = mean_se.first;
        se_gain[k] = mean_se.second;
    }

    const auto best_iterator =
        std::max_element(mean_gain.begin(), mean_gain.end());

    if (best_iterator == mean_gain.end() ||
        !std::isfinite(*best_iterator)) {
        return 0;
    }

    const std::size_t best_k =
        static_cast<std::size_t>(best_iterator - mean_gain.begin());

    if (best_k == 0 || mean_gain[best_k] <= 0.0) return 0;

    const double one_se_threshold =
        mean_gain[best_k] - se_gain[best_k];

    for (std::size_t k = 0; k <= common_max_k; ++k) {
        if (mean_gain[k] >= one_se_threshold) return k;
    }

    return best_k;
}

static std::vector<double> make_cvse_hard_weights(
    const std::vector<CVSEModels> &models,
    std::size_t selected_k) {

    std::vector<double> weights(models.size(), 0.0);
    const auto ranking = rank_cvse_candidates(models);
    const std::size_t effective_k = std::min(selected_k, ranking.size());

    for (std::size_t rank = 0; rank < effective_k; ++rank) {
        weights[ranking[rank]] = 1.0;
    }
    return weights;
}

/*
 * Soft evidence weighting
 * -----------------------
 *
 * Hard CVSE ranks candidates and makes a top-k decision:
 *     weight_q in {0, 1}.
 *
 * Here we retain every sampled/eligible candidate but give it a continuous
 * training-derived weight in [0, 1].  Let g_q be the candidate's mean
 * training gain and g_max the largest positive gain in that training split.
 * For a relative threshold alpha in [0, 1],
 *
 *     tau = alpha * g_max
 *     w_q = max(0, (g_q - tau) / g_max).
 *
 * Thus:
 *   alpha = 1  -> all weights are zero (exact all-tree/null model);
 *   alpha = 0  -> positive-gain candidates are weighted in proportion to
 *                 their training gain;
 *   0 < alpha < 1 -> weak candidates are set to zero and stronger candidates
 *                    are smoothly downweighted instead of hard-selected.
 *
 * The same alpha is comparable across CV folds because it is relative to that
 * split's own strongest positive training gain.
 */
static const std::vector<double> &cvse_soft_weight_alpha_grid() {
    static const std::vector<double> grid{
        1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50,
        0.40, 0.30, 0.20, 0.10, 0.05, 0.00
    };
    return grid;
}

static std::vector<double> make_cvse_soft_weights_from_scores(
    const std::vector<CVSEModels> &models,
    const std::vector<double> &scores,
    double alpha) {

    std::vector<double> weights(models.size(), 0.0);
    if (scores.size() != models.size()) return weights;

    double maximum_positive_gain = 0.0;
    for (std::size_t candidate = 0; candidate < models.size(); ++candidate) {
        if (!models[candidate].eligible) continue;
        const double score = scores[candidate];
        if (!std::isfinite(score)) continue;
        if (score > maximum_positive_gain) {
            maximum_positive_gain = score;
        }
    }

    if (maximum_positive_gain <= 0.0) return weights;

    alpha = std::max(0.0, std::min(1.0, alpha));
    const double threshold = alpha * maximum_positive_gain;

    for (std::size_t candidate = 0; candidate < models.size(); ++candidate) {
        if (!models[candidate].eligible) continue;

        const double score = scores[candidate];
        if (!std::isfinite(score) || score <= threshold || score <= 0.0) {
            continue;
        }

        weights[candidate] =
            std::max(0.0, std::min(
                1.0,
                (score - threshold) / maximum_positive_gain
            ));
    }

    return weights;
}

static std::vector<double> make_cvse_soft_weights(
    const std::vector<CVSEModels> &models,
    double alpha) {

    std::vector<double> scores(models.size(), 0.0);
    for (std::size_t candidate = 0; candidate < models.size(); ++candidate) {
        scores[candidate] = models[candidate].training_gain;
    }

    return make_cvse_soft_weights_from_scores(models, scores, alpha);
}

static double choose_cvse_soft_weight_alpha(
    const std::vector<std::int8_t> &outcomes,
    std::size_t number_of_loci,
    const std::vector<CVSECandidate> &candidates,
    const std::vector<std::size_t> &outer_training_loci,
    const CVSEConfig &config,
    unsigned seed) {

    /*
     * alpha = 1 is the exact all-tree model, so it is the conservative
     * fallback whenever inner CV cannot establish a positive predictive gain.
     */
    if (candidates.empty() || outer_training_loci.size() < 4) return 1.0;

    const std::size_t inner_fold_count =
        std::min(config.inner_folds, outer_training_loci.size());
    if (inner_fold_count < 2) return 1.0;

    const auto inner_folds = make_cvse_folds(
        outer_training_loci, inner_fold_count, seed
    );

    struct InnerFit {
        std::vector<std::size_t> validation_loci;
        std::vector<CVSEModels> models;
        std::vector<std::vector<double>> weights_by_alpha;
    };

    const auto &alpha_grid = cvse_soft_weight_alpha_grid();
    std::vector<InnerFit> inner_fits;
    inner_fits.reserve(inner_folds.size());

    for (std::size_t fold = 0; fold < inner_folds.size(); ++fold) {
        const auto inner_training_loci =
            concatenate_other_cvse_folds(inner_folds, fold);

        InnerFit fit;
        fit.validation_loci = inner_folds[fold];
        fit.models = fit_cvse_models(
            outcomes,
            number_of_loci,
            candidates,
            inner_training_loci,
            config
        );

        fit.weights_by_alpha.reserve(alpha_grid.size());
        for (double alpha : alpha_grid) {
            fit.weights_by_alpha.push_back(
                make_cvse_soft_weights(fit.models, alpha)
            );
        }

        inner_fits.push_back(std::move(fit));
    }

    /*
     * For every alpha, collect paired held-out edge gains.
     *
     * The denominator remains the number of ALL training-eligible candidates
     * observable at this validation locus.  It does NOT depend on alpha or on
     * the number/sum of nonzero weights.  This preserves the conservative
     * edge-level scale and avoids the small-selected-set amplification problem
     * found in the earlier implementation.
     */
    std::vector<std::vector<double>> validation_gains(alpha_grid.size());

    for (const InnerFit &fit : inner_fits) {
        for (std::size_t locus : fit.validation_loci) {
            std::vector<double> candidate_gains(candidates.size(), 0.0);
            std::vector<bool> candidate_observed(candidates.size(), false);
            std::size_t eligible_observed = 0;

            for (std::size_t candidate = 0;
                 candidate < candidates.size(); ++candidate) {
                if (!fit.models[candidate].eligible) continue;

                const int topology = static_cast<int>(
                    outcomes[candidate * number_of_loci + locus]
                );
                if (topology < 0 || topology > 2) continue;

                candidate_observed[candidate] = true;
                ++eligible_observed;

                const double tree_log_probability = std::log(std::max(
                    fit.models[candidate].tree_probability[topology],
                    1.0e-300
                ));
                const double unrestricted_log_probability = std::log(std::max(
                    fit.models[candidate].unrestricted_probability[topology],
                    1.0e-300
                ));

                candidate_gains[candidate] =
                    unrestricted_log_probability - tree_log_probability;
            }

            if (eligible_observed == 0) continue;

            for (std::size_t alpha_index = 0;
                 alpha_index < alpha_grid.size(); ++alpha_index) {
                double weighted_gain = 0.0;
                const auto &weights = fit.weights_by_alpha[alpha_index];

                for (std::size_t candidate = 0;
                     candidate < candidates.size(); ++candidate) {
                    if (!candidate_observed[candidate]) continue;
                    if (weights[candidate] <= 0.0) continue;

                    weighted_gain +=
                        weights[candidate] * candidate_gains[candidate];
                }

                validation_gains[alpha_index].push_back(
                    weighted_gain /
                    static_cast<double>(eligible_observed)
                );
            }
        }
    }

    std::vector<double> mean_gain(
        alpha_grid.size(),
        -std::numeric_limits<double>::infinity()
    );
    std::vector<double> se_gain(alpha_grid.size(), 0.0);

    for (std::size_t alpha_index = 0;
         alpha_index < alpha_grid.size(); ++alpha_index) {
        if (validation_gains[alpha_index].empty()) continue;
        const auto mean_se =
            cvse_mean_and_se(validation_gains[alpha_index]);
        mean_gain[alpha_index] = mean_se.first;
        se_gain[alpha_index] = mean_se.second;
    }

    const auto best_iterator =
        std::max_element(mean_gain.begin(), mean_gain.end());

    if (best_iterator == mean_gain.end() ||
        !std::isfinite(*best_iterator)) {
        return 1.0;
    }

    const std::size_t best_index =
        static_cast<std::size_t>(best_iterator - mean_gain.begin());

    /*
     * alpha = 1 has exactly zero gain.  If no soft-weighted model improves
     * held-out prediction, keep the all-tree model.
     */
    if (mean_gain[best_index] <= 0.0) return 1.0;

    /*
     * One-standard-error rule:
     * choose the MOST strongly regularized model (largest alpha) whose mean
     * held-out gain is within one SE of the best alpha.
     *
     * The grid is ordered from alpha=1 down to alpha=0, so the first model
     * meeting the threshold is the simplest one.
     */
    const double one_se_threshold =
        mean_gain[best_index] - se_gain[best_index];

    for (std::size_t alpha_index = 0;
         alpha_index < alpha_grid.size(); ++alpha_index) {
        if (mean_gain[alpha_index] >= one_se_threshold) {
            return alpha_grid[alpha_index];
        }
    }

    return alpha_grid[best_index];
}


static int classify_cvse_selected_alternative(
    const std::array<double, 3> &counts,
    int expected_topology,
    double pseudocount) {

    const double total = counts[0] + counts[1] + counts[2];
    if (total <= 0.0) return 2;

    std::array<int, 2> wrong_topologies{{-1, -1}};
    int position = 0;
    for (int topology = 0; topology < 3; ++topology) {
        if (topology != expected_topology) {
            wrong_topologies[position++] = topology;
        }
    }

    const auto wrong_first_probability = fit_cvse_tree_probability(
        counts, wrong_topologies[0], pseudocount
    );
    const auto wrong_second_probability = fit_cvse_tree_probability(
        counts, wrong_topologies[1], pseudocount
    );
    const auto unrestricted_probability =
        fit_cvse_unrestricted_probability(counts, pseudocount);

    const double log_n = std::log(std::max(1.0, total));
    const double wrong_first_bic =
        -2.0 * cvse_log_likelihood(counts, wrong_first_probability) + log_n;
    const double wrong_second_bic =
        -2.0 * cvse_log_likelihood(counts, wrong_second_probability) + log_n;
    const double unrestricted_bic =
        -2.0 * cvse_log_likelihood(counts, unrestricted_probability) +
        2.0 * log_n;

    if (wrong_first_bic <= wrong_second_bic &&
        wrong_first_bic <= unrestricted_bic) {
        return 0;
    }
    if (wrong_second_bic <= unrestricted_bic) return 1;
    return 2;
}

CVSEResult run_cvse_for_edge(
    std::vector<Tree *> &gene_trees,
    const std::vector<CVSECandidate> &candidates,
    const CVSEConfig &config) {

    CVSEResult result;
    result.selection_frequency.assign(candidates.size(), 0.0);

    const std::size_t locus_count = gene_trees.size();
    if (candidates.empty() || locus_count < 4) return result;

    /*
     * Build outcomes for every candidate and every locus once. Missing
     * quartets remain -1. Coverage eligibility is always determined inside
     * each training split; validation/test presence never controls eligibility.
     */
    std::vector<std::int8_t> outcomes(
        candidates.size() * locus_count,
        static_cast<std::int8_t>(-1)
    );

    for (std::size_t candidate = 0;
         candidate < candidates.size(); ++candidate) {
        index_t quartet[4] = {
            candidates[candidate].taxa[0],
            candidates[candidate].taxa[1],
            candidates[candidate].taxa[2],
            candidates[candidate].taxa[3]
        };

        for (std::size_t locus = 0; locus < locus_count; ++locus) {
            const index_t topology = gene_trees[locus]->get_quartet(quartet);
            if (topology >= 0 && topology < 3) {
                outcomes[candidate * locus_count + locus] =
                    static_cast<std::int8_t>(topology);
            }
        }
    }

    std::vector<std::size_t> all_loci(locus_count);
    std::iota(all_loci.begin(), all_loci.end(), 0);

    std::vector<double> outer_fold_gains;
    std::vector<double> selected_alpha_values;

    /*
     * Legacy result name:
     *   HardTopK -> selected top-k count
     *   Soft     -> number of candidates with nonzero weight
     */
    std::vector<std::size_t> active_candidate_counts;
    std::vector<std::size_t> selection_counts(candidates.size(), 0);
    std::vector<double> accumulated_selection_weight(candidates.size(), 0.0);

    for (std::size_t repeat = 0; repeat < config.repeats; ++repeat) {
        const auto outer_folds = make_cvse_folds(
            all_loci,
            std::min(config.outer_folds, locus_count),
            config.seed + static_cast<unsigned>(104729 * repeat)
        );

        if (outer_folds.size() < 2) break;

        for (std::size_t fold = 0; fold < outer_folds.size(); ++fold) {
            const auto outer_training_loci =
                concatenate_other_cvse_folds(outer_folds, fold);
            const auto &outer_test_loci = outer_folds[fold];

            const unsigned tuning_seed =
                config.seed + static_cast<unsigned>(
                    13007 * repeat + 8191 * fold + 17
                );

            std::size_t selected_k = 0;
            double selected_alpha = 1.0;

            if (config.weighting_mode == CVSEWeightingMode::HardTopK) {
                selected_k = choose_cvse_exception_count(
                    outcomes,
                    locus_count,
                    candidates,
                    outer_training_loci,
                    config,
                    tuning_seed
                );
            } else {
                selected_alpha = choose_cvse_soft_weight_alpha(
                    outcomes,
                    locus_count,
                    candidates,
                    outer_training_loci,
                    config,
                    tuning_seed
                );
            }

            /* Refit using only the complete outer-training set. */
            const auto outer_models = fit_cvse_models(
                outcomes,
                locus_count,
                candidates,
                outer_training_loci,
                config
            );

            std::vector<double> weights;
            if (config.weighting_mode == CVSEWeightingMode::HardTopK) {
                weights = make_cvse_hard_weights(outer_models, selected_k);
            } else {
                weights = make_cvse_soft_weights(outer_models, selected_alpha);
            }

            const std::size_t active_candidate_count =
                static_cast<std::size_t>(std::count_if(
                    weights.begin(), weights.end(),
                    [](double weight) { return weight > 0.0; }
                ));

            std::vector<double> locus_gains;

            for (std::size_t locus : outer_test_loci) {
                double gain = 0.0;
                std::size_t eligible_observed = 0;

                for (std::size_t candidate = 0;
                     candidate < candidates.size(); ++candidate) {
                    if (!outer_models[candidate].eligible) continue;

                    const int topology = static_cast<int>(
                        outcomes[candidate * locus_count + locus]
                    );
                    if (topology < 0 || topology > 2) continue;

                    ++eligible_observed;
                    if (weights[candidate] <= 0.0) continue;

                    const double candidate_gain =
                        std::log(std::max(
                            outer_models[candidate]
                                .unrestricted_probability[topology],
                            1.0e-300
                        )) -
                        std::log(std::max(
                            outer_models[candidate]
                                .tree_probability[topology],
                            1.0e-300
                        ));

                    gain += weights[candidate] * candidate_gain;
                }

                /*
                 * Fixed denominator: all OUTER-TRAINING-eligible candidates
                 * observable at this test locus, independent of k/alpha and
                 * independent of the sum of selected weights.
                 */
                if (eligible_observed > 0) {
                    locus_gains.push_back(
                        gain / static_cast<double>(eligible_observed)
                    );
                }
            }

            if (!locus_gains.empty()) {
                const auto fold_mean_se = cvse_mean_and_se(locus_gains);
                outer_fold_gains.push_back(fold_mean_se.first);
                active_candidate_counts.push_back(active_candidate_count);

                if (config.weighting_mode == CVSEWeightingMode::Soft) {
                    selected_alpha_values.push_back(selected_alpha);
                }

                for (std::size_t candidate = 0;
                     candidate < weights.size(); ++candidate) {
                    if (weights[candidate] > 0.0) {
                        ++selection_counts[candidate];
                        accumulated_selection_weight[candidate] +=
                            weights[candidate];
                    }
                }
            }
        }
    }

    result.outer_fit_count = outer_fold_gains.size();
    const auto outer_mean_se = cvse_mean_and_se(outer_fold_gains);
    result.mean_outer_gain = outer_mean_se.first;
    result.se_outer_gain = outer_mean_se.second;

    if (!selected_alpha_values.empty()) {
        result.mean_selected_alpha =
            std::accumulate(
                selected_alpha_values.begin(),
                selected_alpha_values.end(),
                0.0
            ) /
            static_cast<double>(selected_alpha_values.size());
    }

    if (!active_candidate_counts.empty()) {
        std::sort(
            active_candidate_counts.begin(),
            active_candidate_counts.end()
        );
        result.median_selected_k =
            active_candidate_counts[active_candidate_counts.size() / 2];
    }

    if (result.outer_fit_count > 0) {
        for (std::size_t candidate = 0;
             candidate < candidates.size(); ++candidate) {
            result.selection_frequency[candidate] =
                static_cast<double>(selection_counts[candidate]) /
                static_cast<double>(result.outer_fit_count);
        }
    }

    const bool any_selected =
        std::any_of(
            selection_counts.begin(),
            selection_counts.end(),
            [](std::size_t count) { return count > 0; }
        );

    /*
     * Contract only when the selected exception model improves held-out
     * composite prediction by more than one standard error.
     */
    result.non_tree_like =
        any_selected &&
        result.outer_fit_count > 0 &&
        result.mean_outer_gain > result.se_outer_gain;

    /* Full-data fits below are diagnostics/annotation only. */
    const auto full_models = fit_cvse_models(
        outcomes,
        locus_count,
        candidates,
        all_loci,
        config
    );

    result.eligible_candidate_count =
        static_cast<std::size_t>(std::count_if(
            full_models.begin(), full_models.end(),
            [](const CVSEModels &model) { return model.eligible; }
        ));
    result.filtered_candidate_count =
        candidates.size() - result.eligible_candidate_count;

    if (any_selected) {
        result.representative_candidate = 0;
        bool found_representative = false;

        for (std::size_t candidate = 0;
             candidate < candidates.size(); ++candidate) {
            if (selection_counts[candidate] == 0) continue;

            if (!found_representative) {
                result.representative_candidate = candidate;
                found_representative = true;
                continue;
            }

            const std::size_t current = result.representative_candidate;
            if (accumulated_selection_weight[candidate] >
                    accumulated_selection_weight[current] ||
                (accumulated_selection_weight[candidate] ==
                     accumulated_selection_weight[current] &&
                 selection_counts[candidate] > selection_counts[current]) ||
                (accumulated_selection_weight[candidate] ==
                     accumulated_selection_weight[current] &&
                 selection_counts[candidate] == selection_counts[current] &&
                 full_models[candidate].training_gain >
                     full_models[current].training_gain)) {
                result.representative_candidate = candidate;
            }
        }

        if (!found_representative) {
            result.representative_candidate =
                std::numeric_limits<std::size_t>::max();
        }
    }

    for (std::size_t candidate = 0;
         candidate < candidates.size(); ++candidate) {
        if (selection_counts[candidate] == 0) continue;

        const auto counts = count_cvse_outcomes(
            outcomes,
            locus_count,
            candidate,
            all_loci
        );

        const int type = classify_cvse_selected_alternative(
            counts,
            candidates[candidate].expected_topology,
            config.pseudocount
        );

        ++result.selected_alternative_types[type];
    }

    return result;
}

std::array<weight_t, 3> full_cvse_counts_for_candidate(
    std::vector<Tree *> &gene_trees,
    const CVSECandidate &candidate) {

    std::array<weight_t, 3> counts{{0, 0, 0}};
    index_t quartet[4] = {
        candidate.taxa[0], candidate.taxa[1],
        candidate.taxa[2], candidate.taxa[3]
    };

    for (Tree *tree : gene_trees) {
        const index_t topology = tree->get_quartet(quartet);
        if (topology >= 0 && topology < 3) counts[topology] += 1;
    }
    return counts;
}


} // namespace cvse

#endif // ENABLE_TOB