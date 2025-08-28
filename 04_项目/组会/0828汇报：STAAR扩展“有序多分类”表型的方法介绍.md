
---
### 一、预处理：REGENENIE控制sample relatedness

1.  **第一步：分解 
    *   我们拿到有序表型，比如“饮酒频率”（从“从不”到“每天”共6个等级）。
    *   我们不直接分析它，而是把它“分解”成5个相关的二分类问题。比如：
        *   问题1：是不是“从不”饮酒？ (是/否)
        *   问题2：饮酒频率是否低于“偶尔”？ (是/否)
        *   ...以此类推。
    *   这样做，是为了让数据能被`regenie`这个工具处理。

2.  **第二步：预测
    *   我们利用`regenie`去分析前面分解出的5个二分类问题，但我们的目的**不是**为了得到最终的基因关联结果。
    *   我们的核心目的是利用`regenie`的`--step 1`功能，为每个样本计算出5个不同的**多基因风险评分（PRS）**。

3.  **第三步：整合
    *   现在，每个样本都有了5个PRS分数，它们之间可能高度相关。
    *   为了提取最关键的遗传信息并避免模型过拟合，我们用**主成分分析**（PCA）对这5个PRS分数进行降维。
    *   最终，我们得到了几个关键的“**PRS主成分**”（比如`prs_pc1`, `prs_pc2`...）。

4.  **第四步：分析
    *   我们除了加入年龄、性别、基因组PC等常规协变量外，还把刚刚计算出的“**PRS主成分**”也放了进去。
    *   这样构建出的“零假设模型”就异常强大，因为它已经提前考虑并校正了群体结构。

---
### **二、fit_nullmodel 增加有序多分类表型的类别

1.  **它首先检查我们要分析的phenotype是什么类型。**
    *   **情况一：** 如果是标准的**二分类**或**连续变量**，这个函数就直接调用 `STAAR` 自己的标准流程来处理，我们完全沿用它成熟的算法。
    *   **情况二：** 如果是我们重点关注的**有序分类变量**，它就会启动我们自己设计的特殊工作流。

2.  **这个workflow是做什么的呢？**
    *   **第一步：拟合模型。** 它会先用一个专门处理有序数据的统计模型（`ordinal::clm`）进行拟合。
    *   **第二步：转换数据。** 它通过一种叫“潜在残差”的统计方法，巧妙地将“轻、中、重”这种有序的等级，**转化**成一个`STAAR`能够理解的、类似于连续变量的“残差”和对应的“权重”。
    *   **第三步：打包伪装。** 最后，它把这些新算出来的值，打包成一个和`STAAR`标准输出一模一样的对象。

---
### 三、有序多分类的workflow

有序回归模型（在代码中是用 `ordinal::clm` 函数拟合的）有一个非常直观的假设：

1.  存在一个未被观测到的、连续的潜在变量 Y*。这个Y*代表了个体真实的、潜在的疾病严重程度或性状水平。
2.  这个潜在变量 Y* 服从一个线性模型：
    `Y* = β₀ + β₁X₁ + ... + βₚXₚ + ε`
    其中，`X`是协变量（如年龄、性别），`β`是它们的效应系数，`ε`是随机误差。
3.  我们观测到的有序分类 Y (例如，疾病等级1, 2, 3) 是由 Y* 和一系列阈值（α）共同决定的。
    *   如果 `Y* ≤ α₁`，我们观测到 Y = 1
    *   如果 `α₁ < Y* ≤ α₂`，我们观测到 Y = 2
    *   如果 `α₂ < Y* ≤ α₃`，我们观测到 Y = 3
    *   ...
*   `link` 参数（如 "probit" 或 "logit"）决定了我们假设误差项 `ε` 服从什么分布。
    *   **probit**：假设 `ε` 服从标准正态分布 N(0, 1)。
    *   **logit**：假设 `ε` 服从标准逻辑斯谛分布。

**计算步骤**：
1.  **拟合模型**：首先，我们用 `clm` 拟合模型，得到协变量的效应 `β` 和阈值 `α`。
2.  **确定残差范围**：对于一个特定的个体，我们知道他的协变量值 `X` 和他最终被观察到的表型分类 `Y=k`。
    *   模型的预测部分是 `η = Xβ`。
    *   根据 `Y=k`，我们知道他的潜在变量 `Y*` 落在 `(α_{k-1}, α_k]` 这个区间内。
    *   因为 `ε = Y* - η`，所以我们同样知道了他的**真实残差 `ε` 必须落在 `(α_{k-1} - η, α_k - η]` 这个区间内**。
3.  **计算条件期望 (即“残差”)**: 使用截断分布的公式来计算 $E[\epsilon_i | a_i < \epsilon_i \le b_i]$。
    *   对于 **probit** 链接 (正态分布)，公式为：
        $Residual_i = \frac{\phi(a_i) - \phi(b_i)}{\Phi(b_i) - \Phi(a_i)}$
    其中 $\phi$ 是标准正态分布的PDF，$\Phi$ 是CDF。
4.  **残差的条件期望是`ε`的最佳统计估计**：在统计理论中，条件期望具有一个非常重要的性质：它是**最小均方误差 (Minimum Mean Square Error, MMSE)** 估计量。
5.  **计算条件方差 (用于“权重”)**: 同样使用截断分布的公式计算 $Var(\epsilon_i | a_i < \epsilon_i \le b_i]$。
    *   对于 **probit** 链接，公式为：
        $Variance_i = 1 + \frac{a_i\phi(a_i) - b_i\phi(b_i)}{\Phi(b_i) - \Phi(a_i)} - (Residual_i)^2$

#### 代码实现

**第 1 步：拟合累积链接模型

```R
clm_obj <- ordinal::clm(formula = formula, data = data, link = "probit", ...)
```

*   CLM正是基于我们前面提到的“潜在变量”理论。它会估计出两组重要的参数：
	1.  **协变量的效应 (`clm_obj$beta`)：** 比如，年龄每增加一岁，那个“潜在的疾病严重程度”会增加多少。
	2.  **阈值 (`clm_obj$alpha`)：** 模型会估计出划分不同等级的“切点”。例如，当潜在严重程度低于`alpha1`时，表现为“轻度”；在`alpha1`和`alpha2`之间时，表现为“中度”，以此类推。
*   我们选择`link = "probit"`是因为它假设潜在变量的误差项服从**标准正态分布**，这是后续计算残差的数学基础。

**第 2 步：计算每个个体的“预测区间”

```R
    # 1. 协变量的线性预测值 (eta)
    eta <- as.vector(X[, names(clm_obj$beta), drop = FALSE] %*% clm_obj$beta)
    
    # 2. 确定每个样本的潜在变量必须落在哪两个阈值之间
    thresholds <- c(-Inf, clm_obj$alpha, Inf)
    y_ordinal_numeric <- as.numeric(clm_obj$y)
    lower_bounds_eta <- thresholds[y_ordinal_numeric]     # 区间下界
    upper_bounds_eta <- thresholds[y_ordinal_numeric + 1] # 区间上界
    ```

**第 3 步：计算残差的条件期望

我们无法直接观测到残差，但我们可以计算出它的“条件期望值”，也就是在已知它落在一个特定区间的前提下，对它的最佳猜测。

```R
# 1. 将潜在变量的区间，转换成残差的区间
lower_bounds_eps <- lower_bounds_eta - eta
upper_bounds_eps <- upper_bounds_eta - eta

# 2. 利用正态分布的数学性质，计算条件期望
phi_a <- dnorm(lower_bounds_eps) # 正态分布密度函数 PDF
phi_b <- dnorm(upper_bounds_eps)
Phi_a <- pnorm(lower_bounds_eps) # 正态分布累积函数 CDF
Phi_b <- pnorm(upper_bounds_eps)

prob_in_interval <- Phi_b - Phi_a # 残差落在这个区间的概率

# 这就是截尾正态分布的期望公式
residuals <- (phi_a - phi_b) / prob_in_interval 
y_numeric <- residuals
```

*   `residuals`的计算公式，是截尾正态分布的期望值。它的直观意义是：“**已知一个服从标准正态分布的随机变量（残差）落在了区间 `[a, b]` 内，那么它的期望值是多少？**”
*   这个计算出的`residuals`就是我们想要的，它是一个连续值，代表了除去所有已知协变量影响后，每个样本剩余的、需要用基因来解释的表型信息。

**第 4 步：计算权重

```R
# 计算截尾正态分布的方差
term1 <- (lower_bounds_eps * phi_a - upper_bounds_eps * phi_b) / prob_in_interval
var_y <- 1 + term1 - residuals^2

# 权重是方差的倒数
weights <- 1 / var_y
```

---

### 四、实际遇到的问题

#### 样本方差为0

```R
--- Ordinal phenotype detected. Using custom ordinal null model workflow. ---
Part 1: Fitting cumulative link model with ordinal::clm...
Part 1: Ordinal model fitting successful.
Part 2: Converting 'clm' object to a GMMAT/STAAR-compatible null model...
Part 2.1: Calculating conditional expectation of latent residuals...
Part 3: Assembling the final 'glmmkin' object...
--- Ordinal null model fitting complete. ---
警告信息:
1: Numerical instability detected in weight calculation.
  - 137565 out of 484058 samples (28.4191%) had near-zero variance for latent residuals.
  - This can happen with extreme covariate values (e.g., PRS) leading to perfect prediction.
  - Their weights have been regularized. If the percentage is high, consider checking for covariate separation. 
2: In fit_null_model(formula = null_model_formula, data = data_for_null_model,  :
  Non-finite values detected in weights. Replacing with 1.
```

**解决：

```R
var_y[var_y < 1e-8] <- 1e-8
```

**关联分析结果
"D:\desktop\multiclass\STAARpipelinePheWAS\ADH1C.xlsx"
