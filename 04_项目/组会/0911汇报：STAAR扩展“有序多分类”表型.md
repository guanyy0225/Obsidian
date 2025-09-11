
---
### 一、导入------方法概括

#### 1. 二分类表型：广义线性混合模型 (GLMM)

*   **广义 (Generalized)**: 因为表型不是标准正态分布的。对于二分类表型，是逻辑回归；对于您的有序表型，理论上可以是**有序probit模型**。
*   **线性 (Linear)**: 协变量和遗传效应仍然是线性组合。
*   **混合 (Mixed)**: 这是最关键的一点。模型中包含**两种效应**：
    *   **固定效应 (Fixed Effects)**: 我们想要直接估计其效应大小的变量，如年龄、性别、基因型 `G` 等。
    *   **随机效应 (Random Effects)**: 用于解释样本之间由于**亲缘关系**和**群体结构**导致的表型相关性。这通常通过一个基于**遗传关系矩阵 (GRM)** 的随机效应项来实现。

#### 2. 有序多分类表型：REGENIE得到的PRS解释样本亲缘关系，从而近似GLMM中的随机效应项。

优势：
- 计算效率高
- 比GLM更准确


---

### 二、具体工作流程

#### 第1步: 得到fullDat：phenotype、covariates、PC1-10

2.  **表型转换**:        ```
    *   **二分类化**: 这是为了给 `REGENIE` 准备输入。代码将6个等级的有序表型，巧妙地转换为了**5个**相关的二分类表型：
        *   `alcohol_bin_1`: 从不 vs. 喝酒 (Y ≤ 1 vs. Y > 1)
        *   `alcohol_bin_2`: 偶尔或更少 vs. 喝得更频繁 (Y ≤ 2 vs. Y > 2)
        *   `alcohol_bin_3`: 每月1-3次或更少 vs. 更频繁 (Y ≤ 3 vs. Y > 3)
        *   ... 以此类推。

3.  **协变量处理**:
    *   创建了年龄的平方项 `age2` 以捕捉非线性效应。
    *   对所有连续的协变量（`age`, `age2`, `PC1-10`）进行了**标准化 (scaling)**。

#### 第2步: 使用 REGENIE 计算个体水平预测值 (PRS)

1.  **运行 REGENIE (Step 1)**:
    *   在这一步中，`REGENIE` 会利用**全基因组**的常见变异信息，对**每一个**二分类饮酒表型分别构建一个复杂的全基因组回归模型。
    *   **关键产物**: `REGENIE` 会为**每个参与者**和**每个二分类表型**，计算出一个**个体水平的遗传预测值**。这个预测值本质上就是一个**多基因风险评分 (PRS)**。
    *   最终得到5个独立的PRS文件。
    * 注意1：这里保留了LOCO文件，LOCO=TRUE时，会给出具体的染色体。
    * 注意2：REGENIE不支持多分类表型，但是支持同时对多个二分类进行分析。

#### 第3步: ==PRS 的正交化== 

1.  **问题**: 这5个PRS分数之间**高度相关**，会引起**多重共线性**。
2.  **解决方案**: 使用**主成分分析 (PCA)**。
3.  **最终数据**: 将这5个新的、正交的PRS主成分 `prs_pc` 作为新的协变量添加到fullDat中，得到最终用于拟合零模型的数据 `data_for_null_model`。

#### 第7步: 拟合 OrdinalSTAAR Null Model

1.  **定义模型**:
    *   **结局变量 (`outcomeCol`)**: `alcohol_intake_frequency` (原始的、6个等级的有序变量)。
    *   **基础协变量 (`covCol`)**: `sex`, `age`, `age2`, `PC1-10`。
    *   **PRS协变量 (`PRSCol`)**: `prs_pc1` 到 `prs_pc5` (5个正交的PRS主成分)。
2.  **调用 `NullModel`**:
    *   代码最后调用了自己开发的 `OrdinalSTAAR::NullModel` 函数。

![[Pasted image 20250910212112.png]]



---
### 三、OrdinalSTAAR::NullModel(...) 

参考：[SurvSTAAR/R/NullModel.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/NullModel.R)

#### Part 1: 拟合有序probit模型 ordinal::clm(...）

#### Part 2: 计算残差和方差组件（[[0805潜在变量残差法]]）

*   为后续的分数检验准备必需的组件。
*   **实现**:
    1.  **计算潜在变量残差**:
        *   `eta <- X_mat %*% alpha_coefs`: 计算每个人的线性预测值 `Xβ`。
        *   `lower_b`, `upper_b`: 根据每个人的观测类别，确定其潜在误差 `ε` 被截尾的区间 `[a, b]`。
        *   `residuals <- (dnorm(lower_b) - dnorm(upper_b)) / prob_interval`: 应用**截尾正态分布的条件期望公式** `[φ(a) - φ(b)] / [Φ(b) - Φ(a)]`，计算出每个人的潜变量残差。
    2.  **计算条件方差和权重**:
        *   `var_y <- 1 + term1 - residuals^2`: 应用**截尾正态分布的条件方差公式**，计算出每个人的潜在变量的条件方差。
        *   `W_mat <- Diagonal(x = 1 / var_y)`: 创建权重矩阵 `W`，其对角线元素是条件方差的倒数。
    3.  **预计算矩阵**:
		- `X_t_W <- crossprod(X_mat, W_mat)`
        * `XWX_mat <- X_t_W %*% X_mat`
        * `XWX_inv <- solve(XWX_mat)`
        * `WX_mat <- t(X_t_W)`
        * 预先计算好分数检验方差公式 `Var(U) = G'WG - G'WX(X'WX)⁻¹X'WG` 中所有与基因型无关的、耗时的部分。

#### Part 3: 组装

1.  `base_list <- list(...)`: 创建一个包含所有核心组件的列表。
2.   根据 `LOCO` 参数，添加 `LOCO = TRUE` 或 `LOCO = FALSE` 标志，并相应地记录染色体 `chr`。这告诉下游函数这个零模型是用于全基因组分析还是LOCO分析。
3.  待完成：如果 `use_SPA = TRUE`，代码会调用 `CGF4LatentRes` 函数来额外计算用于鞍点近似的累积生成函数。
4.  `return(fit_null)`: 返回这个包含了所有信息的、可以直接被 `OrdinalSTAAR` 等函数使用的最终对象。

#### 1. 为什么 SurvSTAAR 需要用 SPA？

因为低事件率导致**鞅残差高度偏态**，破坏了标准检验的正态性假设，所以需要SPA来精确计算尾部概率。

#### 2. 为什么 OrdinalSTAAR 也需要用 SPA？

`OrdinalSTAAR` 处理的是有序分类数据，其核心统计量是基于**潜变量残差**。

*   **潜变量残差的分布特性**:
    *   在您的模型中，我们假设潜在误差 $\epsilon$ 服从标准正态分布。
    *   然而，我们实际使用的**潜变量残差 `residuals`** 并不是真正的 $\epsilon$，而是它的**条件期望 `E[ε | X, Y]`**。
    *   这个条件期望的分布**并不保证是正态的**。它的分布形状完全取决于**有序类别的分布情况**。


---

### 四、OrdinalSTAAR::Ordinal_GeneCentricCoding(...) 

参考：[SurvSTAAR/R/GeneCentricCoding.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/GeneCentricCoding.R)

![[Pasted image 20250910205817.png]]


---
### 五、OrdinalSTAAR::Ordinal_plof(...) 

参考：[SurvSTAAR/R/plof.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/plof.R)

![[Pasted image 20250911085543.png]]

---
### 六、OrdinalSTAAR::OrdinalSTAAR(...) 

参考：[SurvSTAAR/R/SurvSTAAR.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/SurvSTAAR.R)

![[Pasted image 20250911085835.png]]


---
### 七、OrdinalSTAAR::Ordinal**Burden**(...) 

参考：[SurvSTAAR/R/basicFunction.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/basicFunction.R)

**`OrdinalBurden` 和 `SurvSTAAR` 的 `Burden` 在计算公式上，其核心的代数形式是完全相同的，但它们内部变量的统计学含义和来源不同。**


1.  **单变异得分**: 首先，我们为每个变异 $j$ 计算其在零模型下的得分 $U_j$。这个 $U_j$ 衡量了该变异的基因型 $\mathbf{g}_j$ 与模型残差 $\mathbf{e}$ 之间的协方差。
    $$ U_j = \mathbf{g}_j^T \mathbf{e} $$

2.  **Burden得分为加权得分之和**: Burden检验的总得分 $U_B$ 被定义为所有单变异得分的**加权和**：
    $$ U_B = \sum_{j=1}^{p} w_j U_j = \mathbf{w}^T \mathbf{U} $$
    这在数学上与直接使用加权负担基因型得分 $B_i$ 进行检验是等价的。

3.  **方差计算与检验**: 同样地，我们需要计算 $U_B$ 的方差 $\text{Var}(U_B) = \mathbf{w}^T \mathbf{C} \mathbf{w}$（其中 $\mathbf{C}$ 是单变异得分的协方差矩阵），然后构建卡方统计量：
    $$ S_{\text{Burden}} = \frac{U_B^2}{\text{Var}(U_B)} \sim \chi^2_1 $$

真正的区别在于`Score` 向量和 `Covariance` 矩阵。

| 原材料                    | `OrdinalBurden` (来自 `Ordinal_exactScore`) | `Burden` (来自 `SurvSTAAR::exactScore`)     | 核心区别解释                                                                                                                                                                                                                                                              |
| :--------------------- | :---------------------------------------- | :---------------------------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **得分向量 `Score`**       | `Score = G' * residuals_ordinal`          | `Score = G' * residuals_survival`         | **残差的来源完全不同**。<br>- `residuals_ordinal`: 是基于**潜变量模型**和**截尾正态分布**的条件期望计算出的**潜变量残差**。它代表了每个人的“潜在饮酒倾向”与模型预测的偏差。<br>- `residuals_survival`: 是Cox比例风险模型中的**鞅残差 (Martingale Residuals)**。它代表了在观察期间，一个个体“超额”发生事件的数量（观测到的事件数 - 预期的事件数）。                                     |
| **协方差矩阵 `Covariance`** | `Cov(U) = G'WG - G'WX(X'WX)⁻¹X'WG`        | `Cov(U) = G'VG = G'(W-P'P)G - ...` (简化形式) | **方差的结构完全不同**。<br>- **`OrdinalBurden`**: 方差是基于**加权最小二乘法 (WLS)** 的标准方差公式。`W` 是潜变量条件方差的倒数。<br>- **`Burden`**: 方差是基于**Cox模型的偏似然 (Partial Likelihood)** 推导出的复杂方差结构。它不仅包含了类似于 `W` 的权重项（反映了每个时间点的风险人数），还包含一个额外的减项 `-P'P`，这个 `P` 矩阵与随时间变化的**风险集 (risk sets)** 有关，用于处理删失数据。 |

### 八、OrdinalSTAAR::OrdinalACAT(...) 

参考：[SurvSTAAR/R/basicFunction.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/basicFunction.R)

它们计算上的**唯一区别**，发生在需要调用底层 `Burden` 检验的特定情况下。在这个时候，它们会分别调用为各自统计模型（有序 vs. 生存）量身定制的、具有**完全不同统计学基础**的 `Burden` 函数。

---
### 九、OrdinalSTAAR::OrdinalSKAT(...) 

参考：[SurvSTAAR/R/basicFunction.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/basicFunction.R)

这两段代码在处理超稀有变异的复杂策略和高层算法逻辑上是完全相同的。** 它们之间唯一的、也是最关键的区别在于**它们所操作的底层统计量（Score 和 Covariance）的来源和含义不同**，并且在退化为Burden检验时，它们会调用各自模型特有的 `Burden` 函数。

---
### 十、OrdinalSTAAR::OrdinalSTAAR_O(...) 

参考：[SurvSTAAR/R/basicFunction.R at main · Cui-yd/SurvSTAAR](https://github.com/Cui-yd/SurvSTAAR/blob/main/R/basicFunction.R)

---
### 十一、plof结果

![[Pasted image 20250910214240.png]]

### 参考wenjian老师的研究发现

作者将POLMM-GENE应用于**UK Biobank 45万外显子组测序数据**，分析了**五个有序分类表型**：
1.  饮酒频率 (Alcohol intake frequency)
2.  10岁时相对身高 (Comparative height size at age 10)
3.  10岁时相对体型 (Comparative body size at age 10)
4.  睡眠类型/时型 (Morning/evening person chronotype)
5.  认知症状严重程度 (Cognitive symptoms severity)

*   **总共发现了54个显著的基因-表型关联** (p < 2.5 x 10⁻⁶)。
*   **主要发现亮点 (Table 1 & Figure 3)：**
    *   **饮酒频率：** 确认了已知基因`ADH1C`的关联，并发现了一个新的潜在关联基因`GIGYF1`。
    *   **睡眠类型：** 再次验证了生物钟核心基因`PER2`和`PER3`以及褪黑素受体基因`MTNR1B`的强关联。
    *   **10岁时相对身高：** 发现了大量与身高和生长发育相关的基因，如`ACAN`, `NPR2`, `GH1`等，其中许多是已知的人类身高相关基因，验证了方法的可靠性。同时，也发现了一些新的关联。
    *   **认知症状严重程度：** 发现了一个有趣的关联基因`MRGPRX1`，该基因与感知（如瘙痒、疼痛）有关，这为认知功能与感觉通路之间的联系提供了新的线索。