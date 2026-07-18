#ifndef QCF_WRITER_HPP
#define QCF_WRITER_HPP

#include "csvparser.hpp"
#include "tree.hpp"

#include <array>
#include <string>
#include <unordered_set>

class QCFWriter
{
private:
    CSVWriter writer;
    Dict *dict;

    /*
     * Deduplicate by branch ID and sorted quartet.
     */
    std::unordered_set<std::string> seen;

    static std::string make_key(
        index_t branch_id,
        const std::array<index_t, 4> &taxa
    )
    {
        return std::to_string(branch_id) + ":" +
               std::to_string(taxa[0]) + ":" +
               std::to_string(taxa[1]) + ":" +
               std::to_string(taxa[2]) + ":" +
               std::to_string(taxa[3]);
    }

public:
    QCFWriter(
        const std::string &filename,
        Dict *dictionary
    )
        : writer(filename),
          dict(dictionary)
    {
        if (dict == nullptr)
        {
            throw std::invalid_argument(
                "QCFWriter requires a non-null Dict pointer."
            );
        }

        writer.write_row(
            "branch_id",
            "taxon_A",
            "taxon_B",
            "taxon_C",
            "taxon_D",
            "count_AB_CD",
            "count_AC_BD",
            "count_AD_BC",
            "resolved_gene_tree_count"
        );
    }

    void write(
        index_t branch_id,
        const index_t *sorted_indices,
        const std::array<weight_t, 3> &qcfs
    )
    {
        const std::array<index_t, 4> taxa = {
            sorted_indices[0],
            sorted_indices[1],
            sorted_indices[2],
            sorted_indices[3]
        };

        const std::string key = make_key(branch_id, taxa);

        /*
         * The same quartet may be requested repeatedly because of:
         *
         * - p-value caching;
         * - repeated random starts;
         * - neighbor search;
         * - std::sort comparator calls.
         */
        if (!seen.insert(key).second)
        {
            return;
        }

        writer.write_row(
            branch_id,
            dict->index2label(taxa[0]),
            dict->index2label(taxa[1]),
            dict->index2label(taxa[2]),
            dict->index2label(taxa[3]),
            qcfs[0],
            qcfs[1],
            qcfs[2],
            qcfs[0] + qcfs[1] + qcfs[2]
        );
    }

    void flush()
    {
        writer.flush();
    }
};

#endif