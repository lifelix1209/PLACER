下面给出 **PLACER-TE-BAM v5** 的“可落地代码构建级”方案：模块划分、数据结构、单次流式架构、关键算法、阈值与降级策略、以及输出规范。你可以据此直接开仓库写实现，不需要再补概念。

---

# 0. 明确输入、依赖与非目标

## 输入

* `input.bam`：长读对齐结果（ONT/PacBio），必须包含 `SEQ`。
* `ref.fa`：参考基因组 FASTA（用于局部再比对与断点精修；不做全量重比对）。
* `te.fa`：TE 库序列（用于 TE-proxy 和 TE 家族注释）。
* 可选：`ref.fa.fai`、`ref.mmi` 或自建索引（用于局部对齐加速）。
* 可选：BAM 含 `MD` tag（用于错配密度峰加速；没有则禁用该功能）。

## 约束（写进代码与文档）

* **单次顺序遍历 BAM**：禁止全量二次扫描、禁止基于 BAM 回读随机 seek。
* 所有后续“重任务”（局部再比对 / POA / EM）只能使用在 pass 1 已抽取并缓存的片段与元信息。
* “BAM-only”定义为：不需要 FASTQ；但必须提供参考与 TE 库。

## 非目标

* 不承诺在“比对器完全隐身、且无任何 CIGAR/SA/softclip 线索、且局部统计也不异常”的情况下仍能发现插入。该信息论上不可保证。

---

# 1. 代码总体架构（工程视角）

建议三层：

1. **Stream Layer（单遍 BAM 流）**

* 负责：顺序读 BAM、窗口统计、读片段抽取、触发窗口封箱、把任务写入队列。
* 不做：POA、EM、局部对齐（只做 TE-proxy 的超轻量判别也可以在此层完成）。

2. **Task Layer（异步任务队列）**

* 任务类型：

  * `COMPONENT_BUILD`：bin/子群分解、候选落点集合整理
  * `LOCAL_ALIGN`：局部再比对（受限候选集合）
  * `ASSEMBLY`：Graph-POA / 生成 contig
  * `COLLAPSE`：结构级合并
  * `GENOTYPE`：EM(ALT/REF/NULL)+空间先验
* 任务输入必须可序列化（JSON/MsgPack + 二进制片段文件）。

3. **Output Layer**

* 生成两轨输出：

  * `scientific.vcf`（含 locus-set、Tier2 多位点、证据字段）
  * `engineering.tsv/bed`（默认不输出代表坐标；若启用则明确标注“DO_NOT_VALIDATE_BY_SINGLE_COORD”并提供集合）

实现语言建议：C++(htslib)或Python(pysam)+C++扩展。只要遵守单遍原则即可。

---

# 2. 核心数据结构（直接映射到代码类/结构体）

## 2.1 ReadSketch（在 Stream Layer 里产生，轻量）

包含后续所需的最小信息，避免存整条 read：

* `read_id`
* `chrom, pos, end, strand`
* `mapq`
* `flag`（supplementary/secondary/primary）
* `cigar_ops`（压缩存储）
* `has_SA`、`SA_summary`（仅保留断点端点坐标与方向的摘要）
* `tags_present`: {MD, NM, ...}
* `probe_fragments`: list of `ProbeFragment`

## 2.2 ProbeFragment（Gate 1 的输入）

* `source`: END5/END3/SOFTCLIP/CIGAR_INS/CIGAR_DEL/SA_BREAKPOINT/MD_PEAK(optional)
* `seq`（短片段，200–400bp，必要时更短）
* `read_offset_start, read_offset_end`
* `ref_anchor`: (chrom, pos) 近似锚点（用于后续组件归属）

## 2.3 WindowBuffer（滑动窗口缓冲）

* `window_id = (chrom, window_start)`
* `stats`：clip_bp, clip_reads, ins_events, del_events, sa_reads, te_contig_hits, low_mapq_reads
* `buffered_reads`: list of `ReadSketch`（受上限控制，超限则分层抽样保留代表性）
* `triggered`: bool
* `sealed`: bool

## 2.4 Component（事件分解后的最小单元）

* `component_id`
* `locus_set`: list of candidate loci (chrom,pos,score)
* `breakpoint_hypothesis`: orientation, left/right boundary rough
* `core_reads`: list of ReadSketch ids + required fragments pointers
* `all_evidence_summary`: counts/distributions

## 2.5 StructuralRepresentative（结构级合并产物）

* `rep_id`
* `locus_set`
* `structure_signature`：断点几何、TE 家族/方向、插入体长度区间、5’截断级别、倒置标记
* `rep_contig_seq`（POA共识或代表序列）
* `polymorphism_summary`：polyA长度分布、局部SNP/indel摘要（注释用，不进分型核心）

---

# 3. 单次流式遍历：Sensitive 模式可落地实现

## 3.1 关键原则

* 不先做全基因组 W* 再回读。
* 采用“滑动窗口 + 触发封箱 + 任务排队”。

## 3.2 流程（伪代码级）

```
for each alignment in BAM (coordinate-sorted):
    win = get_window(aln.chrom, aln.pos)
    update_stats(win.stats, aln)

    sketches = maybe_extract_readsketch(aln, mode)
    if sketches:
        add_to_buffer(win, sketches)  # with reservoir/stratified sampling

    if not win.triggered and stats_exceed_threshold(win.stats):
        win.triggered = True
        enqueue(GATE1_TASK for win.buffered_reads)

    seal_windows_behind_current_position()
```

### stats_exceed_threshold 具体化

用稳健分位数阈值而不是固定值（跨平台）：

* 每条染色体维护在线估计器（例如 t-digest / P² quantile）记录：

  * clip_bp_per_window
  * large_ins_events_per_window
  * sa_reads_per_window
* 触发条件示例：

  * `clip_bp > Q99_clip` 或
  * `sa_reads > Q99_sa` 或
  * `large_ins_events > Q99_ins` 或
  * `te_contig_hits > 0`（若联合参考）
* fast/sensitive 的区别是阈值分位更激进还是更宽松。

---

# 4. Gate 1（v5 修订版，可实现且不自杀）

Gate 1 的目标：用极低成本把绝大多数非 TE 候选丢掉，同时不依赖“只看两端”。

## 4.1 探针片段集合 P(read)

对每个 ReadSketch：

必选：

* END5, END3（各 200–400bp）

按条件加入：

* SOFTCLIP 片段（若 S≥阈值）
* CIGAR_INS：对每个 I≥I_min 的事件，取其在 read 上的邻域片段（中心±L）
* CIGAR_DEL：同理（用于识别结构异常区域）
* SA_BREAKPOINT：若有 SA，取断点两侧邻域片段

可选（仅 MD 存在）：

* MD_PEAK：从 MD 字符串中恢复错配/缺失位置，滑窗找 top-K 峰，对峰邻域取片段
  无 MD：完全禁用。

## 4.2 TE-proxy 计算（轻量）

不做全 read minimizer 扫描，只对 ProbeFragment：

* 预先为 `te.fa` 构建 `k-mer set` 或 `minimizer sketch`（k=15~19，按错误率选）。
* 对 probe 计算：

  * `hit_count`：命中 TE k-mer 数
  * `hit_density`：hit_count / probe_length
  * `family_votes`：命中到哪类 TE 的计数（用于注释与后续分解）

Gate 1 通过规则（示例）：

* 任一 probe `hit_density >= d_min` 且 `hit_count >= c_min`
* 或多 probe 累积超过阈值
* 并且排除纯低复杂度（polyA/单碱基run）主导的 probe（否则会产生大量伪命中）

通过者进入后续 Full TE scan（仍然是“稀疏”，只在触发窗口中少数 reads 上做）。

---

# 5. Full TE scan（仍然不全量，只在候选上）

对 Gate1 通过的 read：

* 只在“被触发的局部区段”做更强的 TE 匹配：

  * 可用 minimizer chaining 或局部 seed-and-extend
  * 目标：定位 TE-like 区段在 read 上的范围（start/end）与方向
* 输出：

  * `te_segments = [(read_start, read_end, family, strand, score)]`
  * `breakpoint_hints`：与 CIGAR/SA 联合推断的断点大致位置（左右）

---

# 6. Component 构建：locus-set 聚类 + 断点解释分解

在 Task Layer 中对每个触发窗口执行：

## 6.1 初分桶（binning）

按参考坐标把 reads 聚到 bins（200bp~1kb 可调），锚点优先级：

1. SA/split 端点
2. 大 I/D 事件中心
3. 软剪切端点
4. TE segment 附近的参考端点

每个 bin 保留：

* `core_reads`（抽样代表性）
* `all_evidence_summary`

## 6.2 locus-set 构建（受限候选集合）

原则：不做全基因组搜索。

候选落点来源（按可用性）：

* read 的 secondary/supplementary（若存在）
* 同 bin 内其他 reads 的落点集合并集
* TE 反向索引（TE k-mer → genome top-L 候选同源位置）

  * 这里需要一个预计算索引：把参考基因组的 minimizers 建表，但查询只返回 top-L（例如 L≤20），且只在该窗口触发时查询。

若候选集合大于 Lmax：该 bin 直接降级（Tier3候选），后续不做局部再比对。

## 6.3 断点解释分解

对同一 locus-set 内的 reads，再按“断点几何一致性”拆 component：

* 左/右断点的方向与相对位置一致
* TE segment 在 read 上的相对位置一致
* 支持插入结构的 split/ins 模式一致

拆不出来且明显混杂：降级，不组装。

输出：`Component` 列表进入后续。

---

# 7. 局部再比对（Local Re-alignment，受限搜索空间）

对每个 Component：

## 7.1 搜索空间

只允许对 component 的候选落点集合 `{locus_j}` 的每个 locus 在 ±W（例如 10kb）窗口内对齐。
禁止全基因组自由搜索。

## 7.2 输入序列

从 core_reads 的 probe/TE segment 邻域片段中重建“侧翼片段”：

* upstream flank（尽量取 1–2kb）
* downstream flank（同理）
  如果 read 不够长，就尽力取最长可用片段，并记录长度。

## 7.3 对齐器选择

实现上可以：

* 用 minimap2 的 library/API 或调用子进程（但要注意吞吐）
* 或用 edlib/ksw2 做局部半全局比对（窗口小，速度可控）

输出：每条 read 的“落点证据向量”：

* 对各 locus 的分数、边界、方向一致性
* 用于后续 placeability 与 EM。

---

# 8. Assembly：Graph-POA（组件内部，受限输出）

对每个 Component：

## 8.1 组装输入

只用 core_reads，且必须已经通过“断点解释一致性”。
把 reads 切成三段（可缺失）：

* UpFlank
* InsertSegment（TE segment 为核心）
* DownFlank

## 8.2 Graph-POA

* 上/下游侧翼分别做 POA，保留分叉信息，但输出路径数量上限为 2。
* 插入体只重建可支持一致片段；polyA 不作为分叉依据，只记录长度分布。

输出 contigs：

* `contig_k = concat(UpPath_i, InsertRep, DownPath_j)`
* 最多输出 2 条 contig（代表两种主要结构路径）

---

# 9. Contig Collapsing：结构级合并（v5 必需）

对同一 component 产出的多个 contig：

## 9.1 结构签名定义

`structure_signature` 由以下字段构成：

* 参考落点几何：左右侧翼的对齐端点区间、方向
* 插入体长度区间（粗粒度，允许 polyA 波动）
* TE 家族/方向
* 5’截断等级（按长度桶，例如 full / 5’trunc-short / trunc-long）
* 倒置或双断点结构标志

## 9.2 合并规则

若两个 contig 的结构签名等同（或在容忍阈值内），合并为一个 StructuralRepresentative：

* 代表序列：用 POA 共识或选支持度最高 contig
* 微小差异：进入 polymorphism_summary（注释）

只有当结构差异显著（截断点差异大、倒置、有额外断点）时才保留多个代表进入分型。

---

# 10. Placeability 与 Tier 判定

对每个 StructuralRepresentative：

计算 `PlaceabilityScore`：

* `Δ`：最佳 locus 与次佳 locus 的对齐分差
* `SideConsistency`：上/下游侧翼是否一致支持同一 locus
* `SupportConsistency`：支持 reads 在 locus 上的分数分布是否集中
* `CandidateSetSizePenalty`：候选集合越大惩罚越大

Tier：

* Tier1：候选集合小且分差大，侧翼一致
* Tier2：多个等价 locus（分差小或集合大），但结构解释一致
* Tier3：无法稳定落点或组件混杂

断点输出：

* 默认输出区间（跨 reads 的一致区间）
* 只有一致性足够高才输出单点

TSD：

* 仅 Tier1 且断点稳定时检测
* 并做短重复背景显著性过滤（否则 TSD=NA）

---

# 11. Genotyping：EM(ALT/REF/NULL) + 空间先验 + 结构先验（v5 核心）

对每个 Tier1（或 placeability 足够高的 Tier2，按策略）做分型。

## 11.1 组件

* ALT_struct：插入存在，结构为当前 StructuralRepresentative
* REF：无插入
* NULL：背景/其他位点/不可解释

## 11.2 似然项（必须可计算）

对每条 read 计算特征：

* `d_spatial = |x_read_primary - x0|`（x0 为候选位点中心或集合内某 locus）
* `geom_ok`：是否满足断点几何解释（split/ins 模式一致、TE segment 方位一致）
* `align_scores`：对候选 locus 的局部对齐分数与一致性
* `contig_support`：read 是否与插入代表一致（不需要逐碱基高精度，允许用种子一致性+粗对齐）

定义：

* `P(data | ALT)`：需要 geom_ok 高、contig_support 高、且在候选 locus 上对齐一致
* `P(data | REF)`：需要跨断点连续解释更优、无插入证据
* `P(data | NULL)`：geom_ok 低或对齐分散或更符合“远端来源”

## 11.3 先验（解决 nested insertion NULL 吞噬）

空间先验：

* `π_ALT ∝ exp(-d_spatial/λ) * f(geom_ok)`
* `π_REF ∝ exp(-d_spatial/λ) * f(ref_span_ok)`
* `π_NULL ∝ c + (1 - exp(-d_spatial/λ))`

结构先验：

* 如果 read 的 TE family vote 与该位点代表一致，上调 ALT 先验
* 如果 family vote 分散，上调 NULL

## 11.4 EM 更新与收敛

* 初始化：ALT/REF 按直觉计数（split/ins 读给 ALT 初值，连续跨越给 REF）
* E-step：算每条 read 的后验责任度 `γ_alt, γ_ref, γ_null`
* M-step：更新 ALT fraction（并输出 Beta-Binomial 的后验区间）
* 终止：对数似然增益 < ε 或迭代次数达上限（例如 20）

输出：

* `E[ALT], E[REF], E[NULL]`
* `AF = E[ALT]/(E[ALT]+E[REF])` 及可信区间
* `GQ`：基因型似然比转 Phred
* 判定：

  * 若 `E[ALT]` 结构一致且 AF 明确：0/1 或 1/1
  * 若 `E[NULL]` 极高且 `E[ALT]` 无结构一致：NG
  * nested insertion 场景：允许 `E[NULL]` 高但 `E[ALT]` 结构一致仍输出低GQ call（标记 HIGH_BACKGROUND）

---

# 12. 输出规范（避免误导）

## scientific.vcf（主输出）

* Tier1：单坐标
* Tier2：用多记录或集合字段表达

  * 方案 A：每个候选 locus 一条 VCF 记录，共享 `LOCUS_SET_ID`
  * 方案 B：一个主记录 + `ALT_LOCI=chr:pos,chr:pos...`（不推荐给下游工具）
* Tier3：不输出坐标记录，输出到 `evidence.tsv`（含证据ID、TE家族、读支持摘要）

INFO 字段建议：

* `TIER`
* `LOCUS_SET_ID`
* `PLACEABILITY`
* `BP_CI_L/BP_CI_R`
* `TSD_LEN/TSD_SEQ/TSD_BG_P`
* `TE_FAMILY/TE_ORIENT/TE_COV`
* `POLYA_CI`
* `EALT/EREF/ENULL`
* `AF_CI`
* `GQ`
* `FLAGS=HIGH_BACKGROUND|LOW_COMPLEXITY|HOTSPOT|...`

## engineering 输出（默认不提供代表坐标）

* 提供：

  * `LOCUS_SET_ID`
  * `CANDIDATE_LOCI_LIST`
  * `PRIMER_DESIGN_HINT`: “use locus-set list, not single coordinate”
    如果用户显式启用 `--emit-representative-locus`：
* 输出 `REP_LOCUS`，但强制附带 `REP_IS_NOT_GROUND_TRUTH=1`

---

# 13. 性能与复杂度（写进设计文档，避免“代价是谜”）

关键成本来自三处：

1. 单遍 BAM 解压与解析（线性）
2. Gate 1 TE-proxy（只对触发窗口缓冲内少数 read 的少数短片段）
3. 局部再比对与 POA（只对触发组件的 core_reads，且候选 locus 集合受限）

成本控制手段全部显式化：

* 触发阈值是分位统计而不是全量扫描后回读
* 每窗口缓冲与 core_reads 上限
* 候选 locus 集合 Lmax，超出直接降级
* POA 输出路径上限
* EM 迭代次数上限

---

# 14. 实现顺序（最短可用产品到完整体）

1. Stream + WindowBuffer + Trigger（不做 TE-proxy，先只用 CIGAR/SA/clip 触发）
2. Gate1（端部 + clip + I/D 探针）+ TE k-mer 索引
3. Component 构建 + 受限候选集合（先用 secondary/同窗并集，不做 TE反向索引）
4. 局部再比对（候选 locus ±10kb）
5. POA 组装 + 结构级合并
6. Placeability + Tier 输出
7. EM(ALT/REF/NULL) + 空间先验
8. TE反向索引（提升缺 secondary 的场景召回）
9. 可选：MD 峰加速器

---

# 15. 原理总结（每个关键选择的因果链）

* 单遍流式：把 I/O 从“2次解压BAM”降为“1次顺序解压 + 后续只处理缓存片段”，这是大文件的决定性收益。
* Gate1 探针：用少量短片段的 TE-proxy 把 full scan 限制在稀疏集合，避免扫描半个 BAM。
* 受限候选集合：局部再比对只在候选 locus 周围做，避免把工具变成 aligner wrapper。
* 结构级合并：组装可多解，但分型只接受结构假设，不让微小差异拖垮 EM。
* NULL + 空间先验：既避免 closed-world 强行分配，又避免重复区把真实局部信号全吸进 NULL。

以上是 v5 的“能写代码”的版本。按这个拆分写实现，每个模块都有明确输入输出、降级路径和性能边界，不会再出现“理论上很好但工程上写不出来”的空洞环节。



# 16. 开发计划（详细实现版）
Phase 1: 基础设施 (Stream + WindowBuffer + Trigger)

目标：构建单次流式遍历的骨架，实现内存安全的滑动窗口与自适应触发。

实现细节：

BamReader：封装 htslib，启用多线程解压 (hts_set_threads)。必须解析 BAM_FSUPPLEMENTARY (0x800)，仅过滤 Secondary。

WindowBuffer 内存管理：实现“Safe Frontier”逻辑。维护一个 safe_pos 指针，所有 end_pos < safe_pos 的窗口对象必须被物理析构并移出 Map，防止内存随染色体长度线性增长。

自适应阈值 (P² Algorithm)：使用 P² 算法（或 t-digest）在线估算 clip/SA/ins 的 Q99 分位数。避免硬编码阈值导致在不同测序深度下失效。

任务序列化：TaskQueue 必须支持将 Read 片段和元数据序列化到内存块或临时文件（若内存压力大），确保 Worker 线程不需要回读 BAM。

核心约束：

禁止 Seek：BamReader 只能顺序 sam_read1。

O(1) 查找：窗口 ID 使用 (tid << 32) | (pos / window_size) 的整数 Key，禁止字符串拼接。

Phase 2: Gate 1 (TE-proxy 初筛)

目标：在不进行全量比对的前提下，通过廉价探针快速识别 TE 信号，过滤 99% 的背景噪音。

实现细节：

ProbeFragment 提取器：

对每条 Read，提取 5' 和 3' 端部 (200bp)。

扫描 CIGAR，对 I > 50bp 或 S > 200bp 的区域提取邻域序列。

若 MD tag 存在，提取错配密集区（Top-2 peaks）。

TE 索引构建：

对 te.fa 构建 Minimizer Sketch (k=15, w=10)，允许一定容错。不追求唯一比对，只追求“命中密度”。

建立 FamilyMap：Minimizer Hash -> TE Family ID。

Proxy 评分：

计算 Probe 的 Hit Density (matches / bp)。

逻辑门：if (max_density > 0.2 || total_hits > threshold) pass;

核心约束：

稀疏性：严禁对整条 Read 做 Minimizer 扫描。只扫描提取出的 Probe 片段。

低复杂度过滤：在索引构建阶段过滤掉 PolyA/PolyT 的 k-mer，防止由 Homopolymer 造成的假阳性触发。

Phase 3: Component 构建 (Binning & Clustering)

目标：将杂乱的 Reads 整理为可解释的事件簇，并限定搜索空间。

实现细节：

Binning 策略：

以 200bp 为单位网格化。Read 归属由“强锚点”决定（SA断点 > Clip点 > TE-hit点）。

Locus-set 生成 (受限集合)：

收集 Bin 内所有 Reads 的 SA 目标位置。

收集 Secondary Alignment 位置。

(Phase 8) 查询 TE 反向索引。

合并上述位置，去重并聚类（距离 < 500bp 归一）。

一致性分解：

检查 Reads 的方向（Strand）和断点侧（Left/Right）是否矛盾。若矛盾且无法拆分，标记为 AMBIGUOUS_MIXED 并降级。

核心约束：

Lmax 上限：若单 Component 的候选落点集合 size > 20，直接标记为 Tier3，不再进入 Phase 4，防止计算爆炸。

Phase 4: 局部再比对 (Restricted Local Re-alignment)

目标：在受限的候选区域内精确定位断点，不依赖全局比对器。

实现细节：

侧翼提取：从 Core Reads 中提取 Upstream/Downstream Flank (1-2kb)。

受限对齐：

使用 edlib 或 ksw2。Target 仅为 Locus-set 里的候选位置 ±10kb 窗口。

模式：Extension Alignment (半全局)。

证据向量：为每条 Read 生成 vector<Score>，记录其对每个 Locus 的支持度。

核心约束：

禁止全局搜索：严禁加载全基因组索引进行随机比对。必须依赖 Phase 3 提供的候选列表。

Phase 5: Assembly + Collapsing (Graph-POA)

目标：生成代表性序列，并在进入分型前消除微小变异的干扰。

实现细节：

分段组装：将 Reads 切分为 Up-Flank, Insert, Down-Flank 三部分分别组装。

Graph-POA：使用 abPOA 或 spoa。保留分叉图结构，从中提取 Top-2 最优路径。

结构级合并 (Structural Identity)：

定义结构指纹：{Breakpoint_L, Breakpoint_R, TE_Family, Orientation, 5'_Trunc_Level}。

忽略 PolyA 长度差异和 SNP 差异进行合并。

产出 StructuralRepresentative。

核心约束：

深度控制：进入 POA 的 Reads 数量必须有硬上限（如 50），需采用分层抽样（Stratified Sampling）保留多样性。

Phase 6: Placeability + Tier 判定

目标：科学地量化定位的不确定性。

实现细节：

Placeability Score：计算 Score(Best) - Score(2ndBest)。

Side Consistency：检查左侧翼和右侧翼的最佳落点是否指向同一基因组位置（距离 < Gap阈值）。

Tier 规则：

Tier 1: Unique placement (Δ>30), Consistent sides.

Tier 2: Multiple placements (Δ<10), Consistent structure.

Tier 3: Inconsistent or Unmappable.

核心约束：

TSD 验证：仅在 Tier 1 且断点明确时进行 TSD 搜索。需计算 TSD 序列在局部背景中的出现频率，防止将随机微重复误判为 TSD。

Phase 7: Genotyping (EM with Spatial Priors)

目标：在复杂区域准确分型，解决 Nested Insertion 的背景吞噬问题。

实现细节：

模型组件：构建 Mixture Model P(Data) = w_alt * P(D|ALT) + w_ref * P(D|REF) + w_null * P(D|NULL)。

空间先验 (Spatial Prior)：
π<sub>ALT/REF</sub> ∝ exp(−distance/λ)
π<sub>NULL</sub> ∝ 1 − exp(−distance/λ) + base_noise

结构先验：若 Read 的 TE 家族与 Representative 不一致，强制压低其在 ALT 中的似然。

EM 迭代：运行 10-20 轮，输出最终的期望计数 E[ALT], E[REF], E[NULL]。

核心约束：

NULL Sink：NULL 组件必须存在且具备“吸纳远端/杂乱 Reads”的能力。

NG 判定：若 E[NULL] 占比过高（>50%）且 ALT 无显著聚集，输出 GQ=0 或 GT=./.。

Phase 8: 高级优化 (Optional)

目标：提升极端情况下的召回率和速度。

实现细节：

TE 反向索引：预计算 TE k-mer 在参考基因组上的出现位置（仅保留 Top-20 hits）。用于挽救没有 Secondary Alignment 的 Reads。

MD 加速：如果 BAM 头声明了 MD tag，则在 Gate 1 启用 Mismatch 密度峰检测。


# 17. 开发日志

## 2026-02-03
- 项目初始化，创建 README.md 和开发方案文档
- 确定 PLACER 核心架构: 三层架构 (Stream/Task/Output)
- 确定关键技术路线: 单遍流式、Gate1 TE-proxy、受限候选集合、EM genotyping

## 2026-02-04
### Phase 1: 基础设施 (Stream + WindowBuffer + Trigger) - 已完成

#### 已完成:
- [x] 项目结构搭建 (CMakeLists.txt, include/, src/stream/)
- [x] BamReader 类实现 - 单遍 BAM 解析, ReadSketch 结构
- [x] WindowBuffer 类实现 - 滑动窗口, 统计收集, 触发逻辑
- [x] WindowStats 类实现 - P^2 在线分位数估计算法
- [x] Trigger 类实现 - 阈值触发评估
- [x] TaskQueue 类实现 - 异步任务队列框架
- [x] Phase 1 测试用例 - 单元测试 + 集成测试

#### 核心文件:
- `include/bam_reader.h/cpp` - BAM 解析模块
- `include/window_buffer.h/cpp` - 滑动窗口缓冲
- `include/window_stats.h/cpp` - 窗口统计器
- `include/trigger.h/cpp` - 触发器
- `include/task_queue.h/cpp` - 任务队列
- `tests/test_stream.cpp` - Phase 1 测试

#### 测试数据:
- BAM: `/mnt/home1/miska/hl725/projects/tldr_optimized/test/test.bam`
- REF: `/mnt/home1/miska/hl725/projects/tldr_optimized/test/ref.fa`
- 目标区域: `chrTEST:90777-90789`

#### Phase 1 总结

**核心技术指标:**
- 单遍遍历 BAM (无随机 seek)
- P^2 算法在线分位数估算 (clip/SA/ins)
- Safe Frontier 内存管理
- O(1) 窗口查找 (整数 Key)
- 任务序列化支持

**测试结果:**
```
=== PLACER Phase 1 Tests ===
Testing WindowStats... PASS
Testing BamReader... 79 records processed
Testing WindowBuffer... 3 windows created
Testing Trigger... Score=1, triggered
Testing TaskQueue... All tasks processed
Testing integration... PASS
=== All tests passed! ===
```

#### 已知问题与后续优化:
- [x] WindowBuffer Safe Frontier 在染色体切换时需要清理 - 已实现 `flush_all_previous_chromosomes()`
- [x] TaskQueue 序列化到磁盘功能待实现 - 已实现 `TaskSerializer` 类和 `submit_serialized()` 方法
- [后续] Phase 2: Gate 1 TE-proxy 初筛

### 2026-02-04 (Phase 1 问题修复)

#### 已修复问题

1. **WindowBuffer Safe Frontier 染色体切换清理**
   - 新增 `current_chrom_tid_` 成员变量，跟踪当前活跃染色体
   - 新增 `flush_current_chromosome()` 方法，刷新前一个染色体的所有窗口
   - `add_read()` 中检测染色体切换，自动调用刷新
   - `seal_and_flush()` 简化为只接收 `safe_pos`，不再依赖 tid 数值比较
   - **修复理由**：BAM header 中染色体顺序不一定按 tid 数值排列（如 chrMT/tid=0 可能在 chr1/tid=1 之前）

2. **TaskQueue submit_serialized 实际序列化数据**
   - 新增 `SerializedTask` 类，继承 `Task` 并包含 `std::vector<ReadSketch> reads_`
   - `submit_serialized()` 现在使用 `SerializedTask` 而非 `PlaceholderTask`
   - Worker 线程可通过 `get_reads()` 访问任务数据
   - **修复理由**：原实现丢弃了 `reads` 参数，Worker 无法获取数据

3. **TaskSerializer 字节序安全序列化**
   - 新增 `endian` 命名空间，包含 `swap32()`/`swap16()` 和网络字节序转换
   - 所有整数字段（uint32_t, int32_t, uint16_t）写入前转换为网络字节序
   - bool 值显式转换为 0/1 写入（避免平台差异）
   - 新增 `load_task()` 实现反序列化
   - 新增 round-trip 测试验证
   - **修复理由**：原实现使用原生字节序，在大端/小端机器间不兼容

#### 测试结果
```
=== PLACER Phase 1 Tests ===
Testing WindowStats... PASS
Testing BamReader... 79 records processed
Testing WindowBuffer... 3 windows created
Testing Trigger... PASS
Testing TaskQueue... All tasks processed
Testing integration... PASS
=== All tests passed! ===
```

### 2026-02-03 (代码修复)

#### 已修复问题

1. **BamReader SA 解析**
   - 修复指针算术崩溃：`comma1 - strlen(comma1)` → `rname_buffer.assign()`
   - 修复 TID 查找：使用 `sam_hdr_name2tid()` 进行 header 查找
   - 修复位置基准：SA 1-based 转 0-based
   - ReadSketch 复用优化：减少内存分配

2. **Trigger 模块**
   - 修复除零风险：`max(q99 * factor, baseline)` 双重保护
   - 添加基线阈值：`min_baseline_clip = 50.0` 等
   - 合并 evaluate/check_conditions：消除冗余调用

3. **测试与编译**
   - API 同步更新：`read.read_id` → `read.qname`，`read.end` → `read.end_pos`
   - 添加 q99 字段到 WindowStats
   - 修复 WindowBuffer 构造函数与 get_stats 方法

#### 测试结果
```
=== PLACER Phase 1 Tests ===
Testing WindowStats... PASS
Testing BamReader... 79 records processed
Testing WindowBuffer... 3 windows created
Testing Trigger... Score=1, triggered
Testing TaskQueue... All tasks processed
Testing integration... PASS

=== All tests passed! ===
```
