TE 插入片段共识方案 v3.4（TE-guided + Anchor-locked + Deterministic P_i + Deterministic theta）

目标不变：

P_i 对实现细节不敏感

anchor-locked 把 TE 内部重复导致的多解压成唯一且可证伪选择

不要求 TE 分类稳定作为门票，仅要求模板可识别

新增一个你必须承认的目标：

插入方向 theta ∈ {FWD, REV} 必须显式判别，否则 ref_side 和 phi(p) 没有可比性

1. 适用与回退

v3.4 适用于 TE reference 质量高、覆盖全的场景

模板不可识别则回退 v2.1 reference-free

无 QS/QE/QL 或无法一致恢复 query 绝对坐标时：split-SA 仅诊断，不进核心链路

备注：SA tag 的定义里包含 CIGAR，所以理论上可以从 SA 的 CIGAR 推回 query span；但是否可靠取决于比对器对 clipping 的约定一致性，你不能假设总一致，需要一致性检查（见 3.3.2）
samtools.github.io
+1

2. SemiGlobalAlign（写死语义）
   2.1 打分模型

仿射 gap，cost 越小越好（你原样保留）

match = 0

mismatch = 1

gap_open = 2

gap_extend = 1

注：仿射 gap 的 DP 形式是标准套路，不是你发明的，所以别把它当“实现细节”。
Rhodes College Computer Science

2.2 半全局边界（代码级约束）

query 必须全长对齐

target(TE) 只允许 free-end，不允许 free-start

写死：

DP[0][0] = 0

DP[0][j>0] = +inf

终点在末行 i=n 上从 j 中选

2.3 带宽（确定性）

你现在的 band = clamp(ceil(b0 + b1\*sqrt(n)), ...) 没问题，但我要你明确：这不是“算法必然”，只是你用回归版本号固化的工程带宽。真正的理论支撑是“带状 DP 可以把复杂度从 O(mn) 降到 O(k·min(m,n)) 级别”，经典脉络来自 Ukkonen 等近似匹配思想。
ScienceDirect

2.4 成本归一化

保留你定义：

raw_cost = mismatch_cnt1 + gap_open_cnt2 + gap_extend_cnt\*1

aligned_cols = alignment_column_count

norm_cost = raw_cost / max(1, aligned_cols)

3. QS/QE/QL 与 split-SA 坐标
   3.1 QS 生成规则（固定版）

你原样保留即可（略）。

依据：SAM CIGAR 的 H/S 定义、read 长度与对齐区间的关系都来自 SAM 规范，不是你说了算。
samtools.github.io

3.2 产出方式

保留你原样（略）。

3.3 split-SA 对端片段提取（重写为可证伪规则）
3.3.1 先把概念说清楚

对每个 read_id，定义它在“基因组参考”上的一组对齐记录集合 R = {r_k}（primary + supplementary；SA tag 只当线索，不直接当坐标真相）。

我们要从中提取两类片段：

flank 片段：与候选断点 r_bp 邻域共线、提供可靠基因组锚的那段

opposite 片段：与 flank 在 query 上相邻，且在参考侧不共线的那段（通常是 TE 端或远端）

这一步的动机与 SV caller 用 split-read 精确定位断点一致：split-read 和软剪切在单碱基级别提供断点信息，但前提是你能确定“哪段是锚、断点落在锚的哪一端”。DELLY/LUMPY 都把 split-read 当强证据，但它们都显式建模了断点端点与方向，而不是一句“非共线部分”。
OUP Academic
+1

同样，在 MEI/TE caller 里，软剪切和 split-read 是核心信号之一，但工具会把它们和锚定规则一起用，否则假阳性会爆炸。
genome.cshlp.org
+2
OUP Academic
+2

3.3.2 记录标准化与一致性检查（你缺的就是这个）

对每条记录 r_k，计算：

query span：[q_aln0_k, q_aln1_k)（用 3.1 的 QS/QE/softclip 规则）

ref span：[r_aln0_k, r_aln1_k)（按 reference 坐标升序）

strand_k：由 FLAG 0x10 得到（参考 SAM）
samtools.github.io

如果只能从 SA tag 获得某个补充对齐 sa_m，则用 SA 的 CIGAR 计算 q_aln0/q_aln1 的候选值，但必须过一致性门槛：

若同一 read_id 存在显式 supplementary record，与 SA 推回的 span 必须满足：

IoU(query_span) >= 0.9 且 |q_start差| <= 5 且 |q_end差| <= 5

不满足则把 SA 视为“不可靠提示”，不进入核心链路

SA 字段格式来自 SAMtags 规范，你可以引用它来证明你没乱编，但它不保证比对器行为一致，所以你必须做一致性检查。
samtools.github.io
+1

3.3.3 flank 记录选择（确定性）

给定候选断点 r_bp 及窗口 W_ref（版本号控制），定义“候选 flank 集合”：

F = { r_k | r_k.tid == locus_tid AND distance_to_window([r_aln0_k,r_aln1_k), r_bp) <= W_ref }

其中 distance_to_window 指 ref span 与 r_bp 的最近距离。

对每个 r_k ∈ F 定义锚评分：

anchor_len_k：在 r_bp ± W_anchor 内的对齐匹配列数（M/=/X 都算列，I/D 只算 gap 列）

penalty_k = NM_k + 5\*(has_large_indel_near_bp)

score_k = anchor_len_k - alpha\*penalty_k（alpha 固定，版本号给）

选择 flank：

flank = argmax score_k

tie-break（写死）：

更大 anchor_len_k

更高 MAPQ

更小 NM

primary 优先于 supplementary（FLAG 0x800 作为判别）

更小 ref_pos（POS）

3.3.4 opposite 片段提取（连接点附近非共线部分，写成具体算法）

从 R \ {flank} 中挑 mate：

mate = argmax IoU_complement，其中 IoU_complement 指它在 query 上与 flank 的“非重叠长度”最大

要求 nonoverlap_len >= L_min_opp，否则无 opposite

确定 query 上的连接点 q_junc：

若 q_aln1_flank <= q_aln0_mate + mh_max：q_junc = q_aln1_flank

若 q_aln1_mate <= q_aln0_flank + mh_max：q_junc = q_aln0_flank

若两段 overlap：认为存在微同源或重复映射

令 overlap = [max(q0), min(q1)]

若 len(overlap) <= mh_max：q_junc = midpoint(overlap)（向下取整）

否则标记 FAIL_SPLIT_OVERLAP_TOO_LARGE，split-SA 降级 POA 诊断

然后定义 opposite 片段在 query 上的截取区间（就是你要的“连接点附近”）：

opp_span = [q_junc - L_left, q_junc + L_right] ∩ (query_range \ flank_core_range)

L_left/L_right 固定并由版本号给，且必须 ≤ read_len 的合理比例

flank_core_range 定义为 flank 对齐的“高置信核心”区间（去掉两端各 trim_q）

最后，把 opp_span 对应的序列拿去走第 2 节的 TE SemiGlobalAlign，进入 P_i 枚举。

这一步为什么合理：你不是在“整段差集最大”里捞鱼，而是围绕断点邻域抽取可解释的片段，这和 split-read 断点定位的基本思路一致。
OUP Academic
+1

3.3.5 ref_side 规则（你要的确定性定义，给你）

先定义 ref_junc_pos：断点在基因组参考上的坐标（1-based 或 0-based你自行统一，但必须全链路一致）。

计算方式：

如果 flank 在 query 上位于连接点左侧（即 q_junc 接近 q_aln1_flank）：

若 flank 是 forward strand：ref_junc_pos = r_aln1_flank

若 flank 是 reverse strand：ref_junc_pos = r_aln0_flank

如果 flank 在 query 上位于连接点右侧（即 q_junc 接近 q_aln0_flank）：

若 flank 是 forward strand：ref_junc_pos = r_aln0_flank

若 flank 是 reverse strand：ref_junc_pos = r_aln1_flank

然后定义（写死）：

REF_LEFT 当且仅当 ref_junc_pos <= r_bp + eps_bp

REF_RIGHT 当且仅当 ref_junc_pos >= r_bp - eps_bp

若同时满足（只有在 eps 太大或坐标混乱时会发生）：FAIL_REF_SIDE_AMBIG

并加一个你现在缺失的可证伪门槛：

若 |ref_junc_pos - r_bp| > w_side_max：split-SA 不得进 CoreCandidates（仍可进 POA）

4. Deterministic P_i 枚举

你 4.1, 4.2, 4.3 大体可以保留；我只插一句：你用 te_start = te_start_min(j) 的代表选择是可以的，但必须把 span_start(j) 当作“模板端点不确定性”的诊断，否则你只是把不确定性藏起来。abPOA 在工程上能跑得很快，但它不替你解决“多解”这个统计问题。
OUP Academic
+1

5. 模板可识别与多解率

保留你的指标体系（Δ_i、entropy、multi_solution_frac）。我只强化一条硬规则：

start_span_best_i > s_span_max(n)：该片段不得进入 CoreCandidates，也不得影响 theta 判别（见 6.0）

理由：重复导致的多解是 TE 内部最常见的失败模式，xTea 这类工具之所以要组合多种证据，就是因为单一比对端点在重复里经常不唯一。
Nature

6. Anchor-locked（确定性一轮 EM）+ 插入方向 theta 判别（新增关键节）
   6.0 插入方向 theta 的必要性（我明确反对你 v3.3 的写死 phi）

你 v3.3 的

REF_RIGHT -> phi = te_start

REF_LEFT -> phi = te_end

除非你能给出 REF_LEFT/REF_RIGHT 与 TE 方向的严格关系，否则这就是拍脑袋。真实情况是：

插入方向为 FWD 时：左侧基因组 flank 连接 TE 的 5 端，右侧连接 TE 的 3 端

插入方向为 REV 时：左侧连接 TE 的 3 端，右侧连接 TE 的 5 端

这不是观点，是定义层面的事实：你把 TE 模板固定在一个方向时，连接端点会随插入方向翻转。你不引入 theta，就无法把左右两类证据投到同一个 TE 端点坐标系里。

6.1 CoreCandidates（严格，更新）

片段进入 CoreCandidates 必须同时满足：

ref_side 已确定（REF_LEFT/REF_RIGHT），按 3.3.5

anchor_len >= A_min

start_span_best <= s_span_max(n)

片段级模板可识别：δ_tpl >= δ_min

split-SA 片段额外要求：|ref_junc_pos - r_bp| <= w_side_max（3.3.5）

6.2 定义 phi(p, theta)（替换你写死的 phi）

令 TE_len 为模板长度，p 为某个 placement，含 te_start, te_end（均在模板 forward 坐标系里）。

定义：

若 theta = FWD：

REF_LEFT：phi = te_start

REF_RIGHT：phi = te_end

若 theta = REV：

REF_LEFT：phi = te_end

REF_RIGHT：phi = te_start

注意：这里我没有引入坐标翻转 TE_len - x，因为你当前 te_start/te_end 已经是在“模板 forward 坐标”下的端点位置，我们只需要决定“哪个端点代表连接点”。如果你后续要输出带方向的 TE 内部坐标，再另建 phi_oriented。

6.3 用 CoreCandidates 确定 theta0（确定性二选一）

对 theta ∈ {FWD, REV}，只用 CoreCandidates 的最优解 p_i1，计算：

w_i = anchor_len_i \* exp(-norm_cost_i1 / T)

c_theta = weighted_median(phi(p_i1, theta), w_i)

mad_theta = MAD(phi(p_i1, theta))（MAD 用 R 默认常数 1.4826 作为正态一致性尺度，你也可以不缩放，但必须写死）
R Project Search

选择 theta0：

theta0 = argmin_theta [ mad_theta + beta * (1 / sum_i w_i) ]

tie-break 写死：

更小 mad

更大 sum_i w_i

theta = FWD 优先（仅为确定性，不代表更真）

weighted median 的定义是“累积权重过半的分位点”，你别用均值，均值在混簇时会坏得更快。
Wikipedia

6.4 中心初始化 c0（用 theta0）

c0 = weighted_median(phi(p_i1, theta0), w_i)，仅 CoreCandidates

6.5 一轮分配（固定目标，phi 换成带 theta0 的）

对每片段 i（含非 CoreCandidates 但有 P_i）：

p*i = argmin*{p in P_i} [ norm_cost(p) + λ*|phi(p, theta0)-c0| + μ*span_start(p) ]

tie-break 你原样保留即可。

6.6 一轮更新 c1（仅 CoreCandidates）

c1 = weighted_median(phi(p_i, theta0), w_i)，仅 CoreCandidates

固定只做一轮：输出 placement 用 p_i，core 用 c1

输出时必须带上 theta0，否则 downstream 解释不了 phi 是哪个端点

7. breakpoint_core 与窗口（用 theta0 的 phi）

你原结构可以保留，只把 phi(p_i) 全部替换为 phi(p_i, theta0)。

7.1 CoreSet（收紧）

CoreSet = {i in CoreCandidates | |phi(p_i, theta0)-c1| <= w_core_gate}

ambiguous 片段不进 CoreSet

7.2 core 与窗口

te_breakpoint_core = median(phi(p_i, theta0) for i in CoreSet)

core_span = 1.4826 \* MAD(phi(p_i, theta0) for i in CoreSet)
R Project Search

breakpoint_core_window = [core-w, core+w]

w = max(w_min, gamma\*core_span)

8. POA 输入（让 split-SA “必须包含”变成真有用）

你现在的“必须包含 breakpoint_core_window”没问题，但你对 split-SA 的强制包含要加条件，否则就是你同事吐槽的那句“形式上有”。

规则改成：

主簇片段全部进入

若存在通过 3.3.4 与 5 的 split-SA opposite 片段，且其 start_span_best <= s_span_max，则每个 read_id 至少保留 1 条 split-SA opposite（取目标函数最优的那条）

非主簇片段仅在 W_eff 不足时补入，补入上限 M_extra

必须覆盖 breakpoint_core_window

这样 split-SA 进入 POA 的方式是可解释的：它不是“反正塞进去”，而是作为“每条 read 的对端信息”补足共识图的跨断点连通性，这就是 POA/共识方法的核心价值。POA 的经典表述来自 Lee 等人的 partial order graph 多序列比对框架；abPOA 是它的工程加速实现。
OUP Academic
+1

9. QC（必须可证伪）

保留你两层 QC。新增两条必须打出来的诊断字段：

theta0 与 mad_FWD, mad_REV（让人能复核方向判别是不是在胡扯）

split_sa_core_frac：进入 CoreCandidates 的 split-SA 片段比例（低则说明你“必须包含”只是在 POA 补料，别装它是核心定位贡献者）

10. 失败与回退

在你原有基础上加一条：

若 theta0 判别不稳定（例如 mad_FWD 与 mad_REV 差异 < eps_theta，且两边 sum_w 接近），标记 FAIL_THETA_UNCERTAIN：核心窗口仍可给，但不输出方向相关注释，split-SA 降权

11. 接入点（PLACER）

你原接口建议可保留。唯一要求：IAnchorLockedModule 必须输出 theta0 与 ref_junc_pos 相关诊断，否则后面没人能 debug
