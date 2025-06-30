## Abstract
---
全基因组测序（WGS）研究旨在识别与人类疾病相关的rare variants和traits。 
###### 传统的single-marker association analyses在rare variants的检测上缺乏能力，因此研究人员通常使用==variant-set-based analyses==。
- 然而，现有的variant-set-based analyses方法需要预先指定==genetic regions==，不直接适用于 WGS 数据，因为基因组中存在大量non-coding variants的基因间（intergenic）和内含子区域（intron regions）。
- 常用的==sliding-window method==需要预先指定fixed window sizes，这通常是先验未知的，在实践中很难指定。
	- 此外，由于遗传关联区域的大小可能在全基因组和不同的表型之间变化，这种方法也受到限制。
###### SCANG是一种动态扫描统计方法，不需要事先指定窗口大小。
- **动态检测区域**：允许检测出的稀有变异关联区域大小随基因组和表型的不同而变化。
- **错误率控制**：该方法控制全基因组的I型错误率，并考虑遗传变异之间的连锁不平衡。
###### 模拟研究
- 通过大量模拟研究，在多种场景下表明，SCANG在检测稀有变异关联方面明显优于其他多种替代方法，同时控制全基因组范围的I型错误率。
###### 实例分析
- 为展示SCANG的应用，作者分析了来自冠状动脉粥样硬化风险社区研究（ARIC）的WGS脂质数据。

## Introduction
---
#### 全基因组关联研究（GWAS）

- **应用**：过去十年中，GWAS被广泛用于解析复杂疾病和数量性状的遗传结构。
    
- **成就**：GWAS成功识别出成千上万的常见遗传变异，这些变异可能含有与复杂疾病相关的易感等位基因。
    
- **局限**：然而，这些常见变异只解释了疾病遗传力的一小部分。
- GWAS的常见分析策略是进行单个变异测试（individual variant tests），然后进行多重检验调整以识别与性状或疾病相关的变异。然而，单个变异分析不适用于rare variants的分析，因为缺乏统计能力。
###### 近年来，基于全外显子测序（WES）的研究被开展以检测编码序列中的rare variants关联。
- 然而，大多数rare variants是non-coding的，位于内含子和基因间区（introns and intergenic regions），这些区域并不被WES覆盖。
###### 越来越多的全基因组测序（WGS）关联研究正在进行，以识别与复杂性状和疾病相关的rare variants。
1. 最近有很多关于RVs分析方法的文献，这些方法共同测试多个变异在预先指定的遗传集合中的效应，此类方法包括burden tests和non-burden tests如序列核关联测试（SKAT）。
2. 这些测试通过评估多个变异在基因中的累积效应来提高统计能力。
3. 这些variant-set-based tests的一个限制是，它们要求**预先指定用于分析的genetic regions**。
4. 因此，现有的gene-based approaches不适用于分析WGS数据中的所有变异，因为在整个基因组中variant sets并不明确。
###### 滑动窗口方法被提出作为一种识别与复杂性状和疾病相关的RVs的genome regions的方式。
1. 这些程序通过**预先指定一定大小的移动窗口**进行WGS，这些窗口具有固定的跳跃长度，从每个染色体的起始位置开始，沿着基因组移动。
2. 滑动窗口方法的一个主要限制是需要预先指定一个固定的窗口大小，可能是以碱基对或number of variants为单位，但所需的大小通常在事先是未知的，且可能在整个genome和phenotypes中并不最优。
3. 实际上，滑动窗口的大小很难在实践中确定，因为与性状或疾病相关的区域的确切大小可能在基因组中有所不同，并且与表型变化相关。
4. 如果预先指定的窗口大小过大，可能会包括过多的中性变异（neutral variants），从而降低统计能力；如果窗口大小过小，则可能会排除包含关联信号（association signals）的相邻区域。
###### 因此，开发一种dynamic procedure以连续扫描基因组是非常重要的，通过允许rare-variant association testing的单位大小和位置在基因组中变化，从而灵活识别相关区域，而不需要事先指定fixed window的大小和位置。
1. 我们提出使用基于扫描统计（scan-statistic-based）的方法，通过允许不同大小的overlapping windows“向前移动”一个给定大小的窗口，每次移动少量变异，来持续扫描整个基因组，寻找包含关联信号的窗口。
2. Scan statistics是一类广泛的方法，旨在寻找时间和空间中的事件聚集，这些方法已成功应用于遗传研究。
###### 最近，likelihood-ratio-based的scan-statistic procedures被提出，用于识别大遗传区域（如基因）中的clusters of rare-disease variants。
1. 然而，这些likelihood-ratio-based的scan-statistic procedures旨在精炼基因中的疾病聚集区域，而不是用于测试整个基因组的关联。
2. 此外，当前的方法不允许covariates（如age, sex, and population structures），并且只能用于binary traits。
3. 因此，在WGS研究中，迫切需要开发强大的扫描方法，以动态testing rare-variant associations regions的大小及其在基因组中的位置。
##### 在本文中，我们开发了==SCAN the Genome（SCANG）方法==，这是一种灵活且计算效率高的scan-statistic procedure，利用variant-set-based test的p值作为每个移动窗口的扫描统计量，以检测rare-variant association regions，适用于连续性和二元性状。
- SCANG的目标是检测基因组中是否存在任何rare-variant association regions，如果存在，则识别这些关联区域的位置和大小。
###### 1. fit null model
- SCANG首先拟合一个包含协变量（例如年龄、性别和祖先主成分）的null linear or logistic model，但不包括任何遗传变异。
###### 2. set-based tests
1. SCANG对在预先指定的实际窗口范围内的所有可能候选移动窗口进行不同大小的set-based tests。
2. SCANG框架中包括三种测试：
- burden test (SCANG-B）
- 序列核关联测试（SKAT，SCANG-S）
- 高效的全局测试（SCANG-O），通过聚合Cauchy关联测试（ACAT）方法，将burden test和SKAT在不同权重选择下的信息聚合。
1. 所有set-based tests在原假设下共享相同的简化模型，因此拟合的原假设模型是相同的，只需在扫描基因组时拟合一次，这使得SCANG的计算效率非常高。
###### 3. generates an empirical threshold
1. SCANG通过Monte Carlo模拟生成一个经验threshold，以控制给定水平（例如0.05）的全基因组I型错误率。
2. p值小于该threshold的窗口被检测为genome-wisesignificant-associationregions。
3. 提供这些显著窗口的individual-window p values and the genome-wise/family-wise p values。

- 通过模拟，我们证明SCANG在广泛的研究设计中对连续性和二元性状通常比现有方法更具统计能力。
- 我们还将SCANG应用于分析动脉粥样硬化风险社区（ARIC）研究中的WGS和脂质性状。
- 通过允许估计达到全基因组显著性的变异-表型关联区域的最佳大小，SCANG检测到小、致密的低密度脂蛋白胆固醇（sdLDL-c）与位于19号染色体上NECTIN2（MIM: 600798）中的一个4,637 bp区域的稀有变异之间的显著关联，这一关联被传统的滑动窗口方法遗漏。
- 此外，SCANG还检测到比滑动窗口程序检测到的更多稀有变异与脂蛋白（a）（Lpa）之间的关联区域。
# Material and Methods
#### ==SCANG 概述==
1. SCANG是一种可以动态、适应性地检测表型与rare variant regions之间显著关联的方法。
2. 其核心思想是将“变异集合测试（variant-set-based tests）”得到的 p 值作为“扫描统计量”，在基因组范围内对不同大小的窗口进行移动和检测。  
3. 对于每一个候选窗口，SCANG 既计算set-based p value，也计算校正全基因组多重检验（包括重叠窗口和不同大小窗口）的genome-wise p value。  
##### Aggregation Tests for Multiple Variants of a Given Region
###### 1. GLM
###### 2. Burden Test
- 记 Uⱼ = ∑ᵢ (Gᵢⱼ × (yᵢ − bμᵢ)) 为第 j 个变异的得分统计量 (score statistic)
- 将窗口内所有变异的基因型聚合成一个“总变异量”或“总等位基因数”来检测该聚合评分与表型之间的关联，其统计量为Q_burden = ( ∑ⱼ wⱼ Uⱼ )²
- 当 n 较大时，Q_burden 近似服从自由度为 1 的卡方分布。
###### 3. SKAT (Sequence Kernel Association Test)
- 基于方差分量的思想来检测窗口内多变异对表型的整体效应，其统计量为Q_skat = ∑ⱼ (wⱼ Uⱼ)²
- Q_skat的分布是若干卡方分布的混合，亦可得到解析形式的 p 值。
###### 4. 权重选择及 Omnibus 测试
- 对于Burden Test和 SKAT，常用的一种加权方法是利用次要等位基因频率(MAF)在 Beta 分布下的密度函数，如 Beta(MAFⱼ; α₁, α₂)，达到“更罕见的变异给予更高权重”等不同假设。
- 由于真实疾病模型未知，不同权重、不同测试在不同场景下的表现差异很大。为了兼顾多种可能性，文中提出了一种混合（omnibus）测试 Q_omnibus，基于 Cauchy 分布的 ACAT 进行多个测试结果的合并：
1. 首先计算 SKAT 和 Burden 在不同权重设定下 (α₁, α₂) 得到的p_SKAT(α₁, α₂) 和p_burden(α₁, α₂)。
2. 使用 ACAT 将这些 p 值通过 Cauchy 方法进行组合，得到单个统计量 Q_omnibus，并进一步得到相应的 p 值。
3. 与最小 p 值法相比，ACAT 当多个测试的 p 值都相对较小（但没到最小）的情况下能提升统计效能。
4. 此外，相比另一种常见的 SKAT-O 方法，这种 ACAT 法具有更灵活（可组合不同权重）且计算更高效的优点。
##### Dynamic Detection of Rare-Variant Association Regions with Different Window Sizes by Using the Scan Statistic

- 原假设 H₀：整个基因组范围内不存在任何与表型显著关联的变异区域 (r = 0)。
- 备择假设：存在 r 个信号区域，每个区域包含至少一个与表型有显著关联的变异，且不同区域的大小可不同。
- 目标：检测是否存在这些变异-表型关联的区域，如果存在，则估计这些区域的位置与大小。

1. 在原假设下，对所有样本数据(含协变量)拟合广义线性模型，得到 bμᵢ (无遗传效应)。
2. 在基因组范围内，针对一系列可能大小和位置的窗口（这些窗口可以重叠、大小不同），分别进行负担测试、SKAT 或者 Omnibus 等集成方法计算统计量与 p 值。
3. 通过蒙特卡洛模拟（或其他方法）在全基因组层面上设定一个多重检验校正后的阈值，如家族范围 I 型错误率 (FWER) 或全基因组范围错误率 (GWER)。当某窗口 p 值小于该阈值，则判定其为全基因组显著关联区域。
4. 若全基因组的整体检验拒绝 H₀，可将所有显著窗口纳入后续定位和结果报告中，确定关联区域的具体位置与大小。
##### Threshold for Controlling the Genome-Wise (Family-Wise) Type I Error Rate
###### 为什么需要经验阈值 (Empirical Threshold)
1. 在 SCANG 框架下，对全基因组的不同位置、不同大小的滑动窗口进行检测。这些候选窗口会相互重叠，导致它们的统计量（p 值）高度相关。
2. 如果简单地采用 Bonferroni 校正来纠正多重比较，会过于保守，从而大幅降低检验的统计功效(power)。
3. 因此，作者提议利用经验阈值来控制给定水平 (α) 的全基因组（或家族范围）I 型错误率，而不是使用过于保守的 Bonferroni 方法。
###### 如何构造经验阈值

采用了基于 Monte Carlo 的模拟方法，核心思路可分为以下几步：

1. 假设在“全局原假设” (global null) 下，没有任何区域与表型存在关联，也就是说，观测到的变异在理论上都不影响表型，任何显著信号都只是随机噪声。
2. 构造“伪得分向量”（pseudo-score vectors），模拟没有任何真实关联时，各个窗口统计量的分布。
3. 计算伪集合统计量与相应的 p 值
- 对于每个伪得分向量 Û_b，像在真实数据中那样对所有移动窗口或变异集合进行测试，得到一个或多个集合统计量（比如负担测试、SKAT、或 Omnibus），进而得到一系列伪 p 值。