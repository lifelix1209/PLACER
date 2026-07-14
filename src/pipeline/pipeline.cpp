#include "pipeline.h"
#include "conformal_selector.h"
#include "decision_policy.h"
#include "local_interval_cache.h"
#include "parallel_executor_internal.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <htslib/sam.h>
#include <abpoa.h>


namespace placer {

constexpr int32_t kFinalCallDedupDistanceBp = 50;

namespace {

#include "pipeline_local_alignment_helpers.inc"

#include "pipeline_event_helpers.inc"

#include "pipeline_breakpoint_helpers.inc"

#include "pipeline_window_helpers.inc"

}  // namespace

#include "pipeline_call_selection.inc"

#include "pipeline_breakpoint_stage.inc"

#include "pipeline_event_evidence_stage.inc"

#include "pipeline_consensus_stage.inc"

#include "pipeline_segmentation_stage.inc"

#include "pipeline_finalization_stage.inc"

#include "pipeline_entrypoints.inc"

#include "pipeline_hypothesis_emission_stage.inc"

#include "pipeline_bin_processing_stage.inc"

#include "pipeline_factory.inc"

}  // namespace placer
