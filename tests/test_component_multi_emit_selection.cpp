#ifdef NDEBUG
#undef NDEBUG
#endif
#include "pipeline.h"

#include <cassert>
#include <vector>

int main() {
    using namespace placer;

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate truth_near;
        truth_near.pos = 17680845;
        truth_near.score = 4.19922;
        truth_near.emit_te = true;
        candidates.push_back(truth_near);

        ComponentFinalCallCandidate distant_resolved;
        distant_resolved.pos = 17682680;
        distant_resolved.score = 5.95314;
        distant_resolved.emit_te = true;
        candidates.push_back(distant_resolved);

        const auto selected = select_component_final_call_indices(candidates);
        assert(selected.size() == 2);
        assert(selected[0] == 0);
        assert(selected[1] == 1);
    }

    {
        std::vector<ComponentFinalCallCandidate> candidates;

        ComponentFinalCallCandidate weaker;
        weaker.pos = 1000;
        weaker.score = 2.0;
        weaker.emit_te = true;
        candidates.push_back(weaker);

        ComponentFinalCallCandidate stronger_same_locus;
        stronger_same_locus.pos = 1040;
        stronger_same_locus.score = 3.0;
        stronger_same_locus.emit_te = true;
        candidates.push_back(stronger_same_locus);

        ComponentFinalCallCandidate rejected;
        rejected.pos = 1400;
        rejected.score = 9.0;
        rejected.emit_te = false;
        candidates.push_back(rejected);

        const auto selected = select_component_final_call_indices(candidates);
        assert(selected.size() == 1);
        assert(selected[0] == 1);
    }

    return 0;
}
