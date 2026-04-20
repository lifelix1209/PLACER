#include <cassert>
#include <vector>

#include "decision_policy.h"

int main() {
    using namespace placer;

    {
        EventGenotypeInput input;
        input.alt_struct_reads = 3;
        input.ref_span_reads = 3;
        input.min_depth = 3;
        input.min_gq = 20;
        input.error_rate = 0.10;
        input.event_length = 120;
        input.alt_observed_lengths = {118, 121, 123};

        const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);
        assert(decision.best_gt == "0/1");
        assert(decision.depth == 6);
        assert(decision.gq >= 20);
        assert(decision.pass);

        const EventExistenceEvidence evidence = build_event_existence_evidence(input);
        assert(evidence.best_gt == "0/1");
        assert(evidence.gq >= 20);
        assert(evidence.score > 0.0);
    }

    {
        EventGenotypeInput input;
        input.alt_struct_reads = 2;
        input.ref_span_reads = 8;
        input.min_depth = 3;
        input.min_gq = 20;
        input.error_rate = 0.10;
        input.event_length = 120;
        input.alt_observed_lengths = {220, 260};

        const EventGenotypeDecision decision = genotype_event_from_alt_vs_ref(input);
        assert(decision.best_gt == "0/0");
        assert(decision.depth == 10);
        assert(decision.gq == 0);
        assert(!decision.pass);

        const EventExistenceEvidence evidence = build_event_existence_evidence(input);
        assert(evidence.best_gt == "0/0");
        assert(evidence.gq == 0);
        assert(evidence.score <= 0.0);
    }

    return 0;
}
