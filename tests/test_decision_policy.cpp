#include "decision_policy.h"

#include <cassert>

int main() {
    using namespace placer;

    {
        InsertionEvidence e;
        e.tier = 1;
        e.support_reads = 1;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = false;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(d.pass);
        assert(!d.low_confidence);
    }

    {
        InsertionEvidence e;
        e.tier = 3;
        e.support_reads = 10;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = false;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(!d.pass);
    }

    {
        InsertionEvidence e;
        e.tier = 2;
        e.support_reads = 2;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = true;
        p.low_conf_min_support_reads = 3;
        p.low_conf_max_tier = 3;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(d.pass);
        assert(!d.low_confidence);
    }

    {
        InsertionEvidence e;
        e.tier = 3;
        e.support_reads = 5;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = true;
        p.low_conf_min_support_reads = 4;
        p.low_conf_max_tier = 3;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(d.pass);
        assert(d.low_confidence);
    }

    {
        InsertionEvidence e;
        e.tier = 3;
        e.support_reads = 2;
        InsertionAcceptanceParams p;
        p.emit_low_confidence_calls = true;
        p.low_conf_min_support_reads = 4;
        p.low_conf_max_tier = 3;
        const auto d = evaluate_insertion_acceptance(e, p);
        assert(!d.pass);
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = false;
        in.has_proxy_signal = false;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.95;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_CERTAIN");
        assert(!d.promote_top1_name);
        assert(d.qc_proxy_tag.empty());
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = true;
        in.has_proxy_signal = true;
        in.proxy_certain = true;
        in.te_confidence_prob = 0.10;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_CERTAIN");
        assert(d.promote_top1_name);
        assert(d.qc_proxy_tag == "TE_PROXY_CERTAIN");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = true;
        in.has_proxy_signal = true;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.50;
        TeOpenSetParams p;
        p.conf_certain_min = 0.85;
        p.conf_uncertain_min = 0.35;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "TE_UNCERTAIN");
        assert(!d.promote_top1_name);
        assert(d.qc_proxy_tag == "TE_PROXY_WEAK");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = true;
        in.te_gate_uncertain_path = true;
        in.has_proxy_signal = false;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.90;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.status == "NON_TE");
        assert(d.qc_proxy_tag == "TE_PROXY_NONE");
    }

    {
        TeOpenSetInput in;
        in.te_gate_pass = false;
        in.te_gate_uncertain_path = false;
        in.has_proxy_signal = true;
        in.proxy_certain = false;
        in.te_confidence_prob = 0.20;
        TeOpenSetParams p;
        const auto d = classify_te_open_set(in, p);
        assert(d.add_te_gate_fail_tag);
    }

    return 0;
}
