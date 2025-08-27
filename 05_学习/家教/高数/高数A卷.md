
---

## 2019年

### **一、单项选择题**

**1. 曲线 $\sin(xy) + \ln(x - y) = y$ 在点 $(1, 0)$ 处的切线方程是 ( A ).**

> [!NOTE]
> *   **思路:** 使用隐函数求导法求出切线的斜率 $dy/dx$，然后利用点斜式方程。

*   **步骤:**
    1.  对原方程两边关于 $x$ 求导：
        $$
        \cos(xy) \cdot \left(y + x \frac{dy}{dx}\right) + \frac{1}{x-y} \cdot \left(1 - \frac{dy}{dx}\right) = \frac{dy}{dx}
        $$
    2.  将点 $(x, y) = (1, 0)$ 代入上式：
        $$
        \cos(1 \cdot 0) \cdot \left(0 + 1 \cdot \frac{dy}{dx}\right) + \frac{1}{1 - 0} \cdot \left(1 - \frac{dy}{dx}\right) = \frac{dy}{dx}
        $$
        $$
        \cos(0) \cdot \left(\frac{dy}{dx}\right) + 1 \cdot \left(1 - \frac{dy}{dx}\right) = \frac{dy}{dx}
        $$
        $$
        1 \cdot \frac{dy}{dx} + 1 - \frac{dy}{dx} = \frac{dy}{dx}
        $$
        $$
        1 = \frac{dy}{dx}
        $$
    3.  切线的斜率 $m = 1$。
    4.  使用点斜式方程 $y - y₀ = m(x - x₀)$：
        $y - 0 = 1 \cdot (x - 1)$
        $y = x - 1$
*   **答案:** **(A) $y = x - 1$**














**2. 曲线 $y = x² + x$ ($x < 0$) 上曲率为 $\sqrt{2} / 2$ 的点的坐标是 ( B ).**

*   **思路:** 利用曲率公式 $K = \frac{|y''|}{(1 + (y')²)^{3/2}}$，求解 $x$。

*   **步骤:**
    1.  求一阶和二阶导数：
        $y' = 2x + 1$
        $y'' = 2$
    2.  代入曲率公式：
        $K = \frac{|2|}{(1 + (2x + 1)²)^{3/2}}$
    3.  令 $K = \sqrt{2} / 2$：
        $$
        \frac{\sqrt{2}}{2} = \frac{2}{(1 + (2x + 1)²)^{3/2}}
        $$
        $$
        \sqrt{2} \cdot (1 + (2x + 1)²)^{3/2} = 4
        $$
        $$
        (1 + (2x + 1)²)^{3/2} = \frac{4}{\sqrt{2}} = 2\sqrt{2}
        $$
    4.  两边取 $(2/3)$ 次方：
        $$
        1 + (2x + 1)² = (2\sqrt{2})^{2/3} = (8^{1/2})^{2/3} = 8^{1/3} = 2
        $$
        $(2x + 1)² = 1$
    5.  解得 $2x + 1 = 1$ 或 $2x + 1 = -1$。
        $x = 0$ 或 $x = -1$。
    6.  根据题意 $x < 0$，所以 $x = -1$。
    7.  当 $x = -1$ 时，$y = (-1)² + (-1) = 1 - 1 = 0$。
*   **答案:** **(B) (-1, 0)**

---
#### 1. 显函数形式 $y = f(x)$

这是最常用的情况。如果曲线由函数 $y = f(x)$ 给出，其曲率 $K$ 在点 $(x, y)$ 处为：
$$
K = \frac{|y''|}{(1 + (y')^2)^{3/2}}
$$
或者写成 $f(x)$ 的形式：
$$
K = \frac{|f''(x)|}{(1 + [f'(x)]^2)^{3/2}}
$$
*   $y'$ 和 $y''$ (或 $f'(x)$ 和 $f''(x)$) 分别是函数对 $x$ 的一阶和二阶导数。
*   公式中的绝对值 $|y''|$ 确保曲率 $K$ 是一个非负值。

**曲率半径 (Radius of Curvature)**
曲率半径 $R$ 是曲率的倒数，它表示在该点最接近曲线的圆（称为“密切圆”）的半径。
$$
R = \frac{1}{K} = \frac{(1 + (y')^2)^{3/2}}{|y''|}
$$
---

#### 2. 参数方程形式 $\{ x = x(t), y = y(t) \}$

如果曲线由参数方程给出，其曲率 $K$ 在参数 $t$ 对应的点处为：
$$
K = \frac{|x'y'' - x''y'|}{((x')^2 + (y')^2)^{3/2}}
$$
*   $x'$, $x''$ 分别是 $x(t)$ 对 $t$ 的一阶和二阶导数。
*   $y'$, $y''$ 分别是 $y(t)$ 对 $t$ 的一阶和二阶导数。

**这个公式更通用**。实际上，显函数形式可以看作是参数方程 $x = t, y = f(t)$ 的一个特例。将 $x' = 1$, $x'' = 0$, $y' = f'(t)$, $y'' = f''(t)$ 代入参数方程的曲率公式，就能推导出显函数形式的公式。

---











**3. 极限 $\lim_{n \to \infty} \ln\left[\left(1 + \frac{1}{n}\right)^2\left(1 + \frac{2}{n}\right)^2\cdots\left(1 + \frac{n}{n}\right)^2\right]^{1/n}$ = ( B ).**

*   **思路:** 将极限表达式化为定积分的形式。
*   **步骤:**
    1.  设 $L$ 为所求极限。利用对数性质化简：
        $L = \lim_{n \to \infty} \frac{1}{n} \ln\left[\left(1 + \frac{1}{n}\right)^2\left(1 + \frac{2}{n}\right)^2\cdots\left(1 + \frac{n}{n}\right)^2\right]$
        $L = \lim_{n \to \infty} \frac{1}{n} \left[\ln\left(1 + \frac{1}{n}\right)^2 + \ln\left(1 + \frac{2}{n}\right)^2 + \dots + \ln\left(1 + \frac{n}{n}\right)^2\right]$
        $L = \lim_{n \to \infty} \frac{1}{n} \left[2\ln\left(1 + \frac{1}{n}\right) + 2\ln\left(1 + \frac{2}{n}\right) + \dots + 2\ln\left(1 + \frac{n}{n}\right)\right]$
        $L = \lim_{n \to \infty} \frac{2}{n} \sum_{k=1}^{n} \ln\left(1 + \frac{k}{n}\right)$
    2.  根据定积分的定义，这可以看作是对函数 $f(x) = 2\ln(1+x)$ 在区间 $[0, 1]$ 上的黎曼和：
        $L = \int_{0}^{1} 2\ln(1+x) dx$
    3.  为了匹配选项，进行变量代换。令 $u = 1 + x$，则 $du = dx$。当 $x = 0$ 时, $u = 1$；当 $x = 1$ 时, $u = 2$。
        $L = \int_{1}^{2} 2\ln(u) du$
    4.  将变量 u 写回 x：
        $L = 2 \int_{1}^{2} \ln(x) dx$
*   **答案:** **(B) $2\int_{1}^{2} \ln x dx$**









**4. 设 C 为任意常数, 则微分方程 $\frac{dy}{dx} + \frac{x}{y} = 0$ 的通解是 ( A ).**

> [!NOTE]
> *   **思路:** 这是一个可分离变量的微分方程。

*   **步骤:**
    1.  分离变量：
        $\frac{dy}{dx} = -\frac{x}{y}$
        $y \, dy = -x \, dx$
    2.  两边积分：
        $\int y \, dy = \int -x \, dx$
        $\frac{y^2}{2} = -\frac{x^2}{2} + C_1$ ($C_1$ 是积分常数)
    3.  整理方程：
        $\frac{y^2}{2} + \frac{x^2}{2} = C_1$
        $x^2 + y^2 = 2C_1$
    4.  令 $C^2 = 2C_1$ (因为 $2C_1$ 是一个任意正的常数，可以用 $C^2$ 表示)。
        $x^2 + y^2 = C^2$
*   **答案:** **(A) $x^2 + y^2 = C^2$**

**解题路线图：如何识别方程类型**

当你遇到一个一阶微分方程时，按以下顺序进行检查：

1.  **真的是不可分离吗？**
    *   尝试代数变形，看是否能分离。如果不行，进入下一步。

2.  **是一阶线性方程 (First-Order Linear) 吗？**
    *   标准形式：$\frac{dy}{dx} + P(x)y = Q(x)$
    *   这是最常见、最重要的类型之一。

3.  **是齐次方程 (Homogeneous) 吗？**
    *   标准形式：$\frac{dy}{dx} = F\left(\frac{y}{x}\right)$
    *   标志是方程中的 $x$ 和 $y$ 可以组合成 $y/x$ 的形式。

---

**不可分离方程的解法**

**1. 一阶线性方程 (First-Order Linear)

*   **标准形式**: $y' + P(x)y = Q(x)$
*   **核心方法**: **积分因子法 (Integrating Factor)**

**求解步骤:**
1.  **化为标准形式**，找出 $P(x)$ 和 $Q(x)$。
2.  **计算积分因子 $\mu(x)$**:
    $\mu(x) = e^{\int P(x)dx}$
    (计算积分时，暂时忽略常数 C)
3.  **将标准方程两边同乘以 $\mu(x)$**:
    $\mu(x)y' + \mu(x)P(x)y = \mu(x)Q(x)$
    神奇之处在于，方程的左边现在正好是 $(\mu(x)y)$ 的导数，即 $(\mu(x)y)'$。
4.  **方程变为**: $(\mu(x)y)' = \mu(x)Q(x)$
5.  **两边对 $x$ 积分**:
    $\int (\mu(x)y)' dx = \int \mu(x)Q(x) dx$
    $\mu(x)y = \int \mu(x)Q(x) dx + C$
6.  **解出 y**:
    $y = \frac{1}{\mu(x)} \left[ \int \mu(x)Q(x) dx + C \right]$

**示例**: 求解 $xy' - 2y = x^2$
1.  标准形式: $y' - \frac{2}{x}y = x$。这里 $P(x) = -2/x$, $Q(x) = x$。
2.  积分因子: $\mu(x) = e^{\int(-2/x)dx} = e^{-2\ln|x|} = e^{\ln(x^{-2})} = x^{-2} = 1/x^2$。
3.  乘以 $1/x^2$: $\frac{1}{x^2}y' - \frac{2}{x^3}y = \frac{1}{x}$。
4.  左边是 $(\frac{1}{x^2}y)'$。所以 $(\frac{1}{x^2}y)' = \frac{1}{x}$。
5.  积分: $\frac{1}{x^2}y = \int \frac{1}{x} dx = \ln|x| + C$。
6.  解出 y: $y = x^2(\ln|x| + C)$。


**2. 齐次方程 (Homogeneous)

*   **标准形式**: $\frac{dy}{dx} = F\left(\frac{y}{x}\right)$
*   **核心方法**: **变量代换 $v = y/x$**

**求解步骤:**
1.  **将方程变形**，确认它可以写成 $F(y/x)$ 的形式。
2.  **进行代换**: 令 $v = y/x$，则 $y = vx$。
3.  **求 $\frac{dy}{dx}$**: 对 $y = vx$ 使用乘法法则求导：$\frac{dy}{dx} = \frac{d}{dx}(vx) = v + x\frac{dv}{dx}$。
4.  **代入原方程**:
    $v + x\frac{dv}{dx} = F(v)$
5.  **分离变量**:
    $x\frac{dv}{dx} = F(v) - v$
    $\frac{dv}{F(v) - v} = \frac{dx}{x}$
    现在这个关于 $v$ 和 $x$ 的方程是**可分离的**了。
6.  **积分求解** $v$。
7.  **换回 $y$**: 将 $v = y/x$ 代回去，得到最终解。

**示例**: 求解 $\frac{dy}{dx} = \frac{x+y}{x}$
1.  变形: $\frac{dy}{dx} = 1 + \frac{y}{x}$。这是 $F(y/x)$ 的形式，其中 $F(v) = 1+v$。
2.  代换 $y = vx$， $\frac{dy}{dx} = v + x\frac{dv}{dx}$。
3.  代入: $v + x\frac{dv}{dx} = 1 + v$。
4.  分离变量: $x\frac{dv}{dx} = 1 \implies dv = \frac{dx}{x}$。
5.  积分: $\int dv = \int \frac{dx}{x} \implies v = \ln|x| + C$。
6.  换回 y: $y/x = \ln|x| + C \implies y = x(\ln|x| + C)$。





---

### **二、填空题**

**1. 极限 $\lim_{x \to 0} \frac{e^x - e^{\sin x}}{x - \sin x}$ = __1__**

> [!NOTE]
> *   **思路:** 这是一个 "0/0" 型的极限，可以使用洛必达法则、等价无穷小或泰勒展开）。使用一个巧妙的代换会更简单。

**方法1 (观察):**
令 $f(t) = e^t$。极限可以看作 $\lim_{x \to 0} \frac{f(x) - f(\sin x)}{x - \sin x}$。
由于当 $x \to 0$ 时, $\sin x \to 0$。
根据导数的定义 $f'(a) = \lim_{b \to a} \frac{f(b) - f(a)}{b - a}$，该极限等于 $f'(0)$。
$f'(t) = e^t$，所以 $f'(0) = e^0 = 1$。






**方法2 (洛必达):**
$\lim_{x \to 0} \frac{e^x - e^{\sin x} \cos(x)}{1 - \cos(x)}$
利用等价无穷小 $1-\cos(x) \sim \frac{1}{2}x^2$。
$\lim_{x \to 0} \frac{e^x - e^{\sin x} \cos(x)}{\frac{1}{2}x^2}$
再用一次洛必达：
$\lim_{x \to 0} \frac{e^x - e^{\sin x}\cos^2(x) + e^{\sin x}\sin(x)}{x} = \lim_{x \to 0} \left(\frac{e^x - e^{\sin x}\cos^2(x)}{x} + e^{\sin x}\frac{\sin x}{x}\right)$
第二项是 $e^0 \cdot 1 = 1$。
第一项 $\lim_{x \to 0} \frac{e^x - e^{\sin x}\cos^2(x)}{x}$ 还是 "0/0" 型，再用一次洛必达：
$\lim_{x \to 0} \frac{e^x - (e^{\sin x}\cos^3(x) - e^{\sin x}2\cos(x)(-\sin(x)))}{1} = e^0 - (e^0 \cdot 1 - 0) = 1-1 = 0$.
所以结果是 $0+1=1$。**方法1和2的结果是正确的。**







**方法3 (泰勒展开):**

我们将对分子和分母分别进行泰勒展开（在 $x=0$ 附近），直到找到它们各自的**首个非零项**。然后用这两个首项的比值来确定极限。

**第一步：展开分母**

分母 $x - \sin x$ 。
我们需要 $\sin x$ 的泰勒展开式：
$$
\sin x = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \dots = x - \frac{x^3}{6} + O(x^5)
$$
代入分母：
$$
x - \sin x = x - \left(x - \frac{x^3}{6} + O(x^5)\right) = \frac{x^3}{6} + O(x^5)
$$
**关键信息**：分母的首个非零项是 $\frac{x^3}{6}$。这意味着，为了得到一个确定的、非零的极限，我们也必须将分子展开到 $x^3$ 项。

**第二步：展开分子**

分子是 $e^x - e^{\sin x}$。

1.  **展开 $e^x$:**
    $$
    e^x = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + O(x^4) = 1 + x + \frac{x^2}{2} + \frac{x^3}{6} + O(x^4)
    $$

2.  **展开 $e^{\sin x}$:**
    这是一个复合函数的展开。我们使用 $e^u$ 的展开式，其中 $u = \sin x$。
    $$
    e^u = 1 + u + \frac{u^2}{2!} + \frac{u^3}{3!} + O(u^4)
    $$
    现在，我们将 $\sin x$ 的展开式代入 $u$。
    $$
    u = \sin x = x - \frac{x^3}{6} + O(x^5)
    $$
    接着，我们计算 $u^2, u^3, \dots$ 并保留到 $x^3$ 项：
    *   $u = x - \frac{x^3}{6} + O(x^5)$
    *   $u^2 = \left(x - \frac{x^3}{6} + O(x^5)\right)^2 = x^2 - 2 \cdot x \cdot \frac{x^3}{6} + O(x^{10}) = x^2 + O(x^4)$
    *   $u^3 = \left(x - \frac{x^3}{6} + O(x^5)\right)^3 = x^3 + O(x^5)$

    现在，将这些代回到 $e^u$ 的展开式中：====
    $$
    e^{\sin x} = 1 + \left(x - \frac{x^3}{6}\right) + \frac{1}{2}\left(x^2 + O(x^4)\right) + \frac{1}{6}\left(x^3 + O(x^5)\right) + O(x^4)
    $$
    整理并按 $x$ 的幂次排序：
    $$
    e^{\sin x} = 1 + x + \frac{x^2}{2} + \left(-\frac{x^3}{6} + \frac{x^3}{6}\right) + O(x^4)
    $$
    **注意这里！** 来自 $u$ 的 $-\frac{x^3}{6}$ 项和来自 $\frac{u^3}{6}$ 的 $+\frac{x^3}{6}$ 项**恰好相互抵消**了。这是一个非常关键的细节，之前的错误计算很可能忽略了这一点。
    所以，$e^{\sin x}$ 的展开式是：
    $$
    e^{\sin x} = 1 + x + \frac{x^2}{2} + O(x^4)
    $$

3.  **计算分子的差：**
    $$
    e^x - e^{\sin x} = \left(1 + x + \frac{x^2}{2} + \frac{x^3}{6} + O(x^4)\right) - \left(1 + x + \frac{x^2}{2} + O(x^4)\right)
    $$
    $$
    e^x - e^{\sin x} = \frac{x^3}{6} + O(x^4)
    $$
    分子的首个非零项是 $\frac{x^3}{6}$。

**第三步：计算极限**

现在我们将分子和分母的展开结果代回原极限表达式：
$$
\lim_{x \to 0} \frac{e^x - e^{\sin x}}{x - \sin x} = \lim_{x \to 0} \frac{\frac{x^3}{6} + O(x^4)}{\frac{x^3}{6} + O(x^5)}
$$
当 $x \to 0$ 时，所有高阶无穷小项 $O(x)$ 和 $O(x^2)$ 都趋近于0。
$$
= \frac{\frac{1}{6}}{\frac{1}{6}} = 1
$$

**方法4：等价无穷小**

这种方法可以看作是泰勒展开的快捷方式，能有效避免复杂的代数计算。

1.  对分子进行变形：
    $$
    e^x - e^{\sin x} = e^{\sin x} \left( \frac{e^x}{e^{\sin x}} - 1 \right) = e^{\sin x} (e^{x - \sin x} - 1)
    $$
2.  代回极限表达式：
    $$
    \lim_{x \to 0} \frac{e^{\sin x} (e^{x - \sin x} - 1)}{x - \sin x}
    $$
3.  将极限拆分：
    $$
    = \left( \lim_{x \to 0} e^{\sin x} \right) \cdot \left( \lim_{x \to 0} \frac{e^{x - \sin x} - 1}{x - \sin x} \right)
    $$
4.  分别计算两个极限：
    *   第一个极限：当 $x \to 0$ 时，$\sin x \to 0$，所以 $\lim_{x \to 0} e^{\sin x} = e^0 = 1$。
    *   第二个极限：这是一个标准形式 $\lim_{u \to 0} \frac{e^u - 1}{u} = 1$。在这里，令 $u = x - \sin x$。当 $x \to 0$ 时，显然 $u \to 0$。所以这个极限的值也是 1。

5.  最终结果：
    $$
    1 \cdot 1 = 1
    $$




**2. 曲线 $y = x³ - 3x² + 6x - 5$ 的拐点为 __(1, -1)__**

> [!NOTE]
> *   **思路:** 拐点出现在二阶导数等于零或不存在的地方，并且二阶导数在该点两侧变号。

*   **步骤:**
    1.  求一阶导数: $y' = 3x² - 6x + 6$
    2.  求二阶导数: $y'' = 6x - 6$
    3.  令 $y'' = 0$，得 $6x - 6 = 0$，解得 $x = 1$。
    4.  检查 $x = 1$ 两侧 $y''$ 的符号：
        当 $x < 1$ 时, $y'' < 0$ (曲线凹)
        当 $x > 1$ 时, $y'' > 0$ (曲线凸)
        符号发生改变，因此 $x = 1$ 是拐点的横坐标。
    5.  计算 y 坐标: $y(1) = 1³ - 3(1)² + 6(1) - 5 = 1 - 3 + 6 - 5 = -1$。
*   **答案:** **(1, -1)**

---
##### 驻点

对于一个可导函数 $y = f(x)$，如果在其定义域内的某一点 $x_0$ 满足：
$f'(x_0) = 0$
那么，点 $(x_0, f(x_0))$ 就是函数 $f(x)$ 的一个**驻点**。

**核心思想**：在驻点处，函数的变化率（即导数）为零。这意味着函数在这一点“暂时停止”了上升或下降的趋势。

---

### **驻点的三种类型**

驻点虽然都满足 $f'(x) = 0$，但它们周围的函数行为可以完全不同。主要有以下三种类型：

#### 1. **局部极大值点 

*   **特征**: 在这一点，函数值比它紧邻的左右两边的任何点的值都**大**。
*   **导数行为**: 在该点左侧，$f'(x) > 0$ (函数上升)；在该点右侧，$f'(x) < 0$ (函数下降)。导数符号由正变负。
*   **二阶导数判别法**: 如果 $f''(x_0) < 0$，则该驻点是局部极大值点（开口向下）。
*   **图像**: 看起来像一座**山峰**。

#### 2. **局部极小值点

*   **特征**: 在这一点，函数值比它紧邻的左右两边的任何点的值都**小**。
*   **导数行为**: 在该点左侧，$f'(x) < 0$ (函数下降)；在该点右侧，$f'(x) > 0$ (函数上升)。导数符号由负变正。
*   **二阶导数判别法**: 如果 $f''(x_0) > 0$，则该驻点是局部极小值点（开口向上）。
*   **图像**: 看起来像一个**山谷**。

#### 3. **鞍点

*   **特征**: 在这一点，函数既不是局部最大值也不是局部最小值。
*   **导数行为**: 导数在该点两侧**符号不变**。
    *   例如，在 $f(x) = x^3$ 的 $x=0$ 点，左侧 $f'(x) > 0$，右侧 $f'(x) > 0$。函数一直在上升，只是在 $x=0$ 这一点“歇了口气”，切线变平了。
*   **二阶导数判别法**: 如果 $f''(x_0) = 0$，则无法判断，需要进一步分析（如观察更高阶导数或 $f'(x)$ 在该点两侧的符号）。如果 $f''(x)$ 在该点两侧变号，则该点是拐点。
*   **图像**: 看起来像一个**水平的 “S” 形**或者马鞍的中心部分。

---
### **驻点 vs. 临界点

*   **驻点**: **$f'(x) = 0$** 的点。
*   **临界点**: **$f'(x) = 0$ 或 $f'(x)$ 不存在**的点。

也就是说，**所有的驻点都是临界点，但临界点不一定是驻点**。

**一个 $f'(x)$ 不存在的临界点例子**:
函数 $f(x) = |x|$。在 $x=0$ 处，函数图像有一个尖点，导数不存在。因此，$x=0$ 是一个临界点，但它不是驻点（因为切线不是水平的）。这个点也是一个局部极小值点。






**3. 不定积分 $\int \frac{5x - 1}{x^2 - x - 2} dx$ = __$3 \ln|x - 2| + 2 \ln|x + 1| + C$__**

*   **思路:** 使用部分分式分解法。
*   **步骤:**
    1.  因式分解分母: $x^2 - x - 2 = (x - 2)(x + 1)$
    2.  设 $\frac{5x - 1}{(x - 2)(x + 1)} = \frac{A}{x - 2} + \frac{B}{x + 1}$
    3.  通分得到: $5x - 1 = A(x + 1) + B(x - 2)$
    4.  解出 A 和 B:
        令 $x = 2$: $5(2) - 1 = A(2 + 1) \implies 9 = 3A \implies A = 3$
        令 $x = -1$: $5(-1) - 1 = B(-1 - 2) \implies -6 = -3B \implies B = 2$
    5.  积分: $\int \left(\frac{3}{x - 2} + \frac{2}{x + 1}\right) dx = 3\int \frac{1}{x-2} dx + 2\int \frac{1}{x+1} dx = 3 \ln|x - 2| + 2 \ln|x + 1| + C$

> [!NOTE]
> #### 1. **换元积分法**

*   **第一类换元法 (凑微分法)**:
    *   **核心思想**: 观察被积函数，将其凑成 $\int f(g(x)) g'(x) dx$ 的形式，然后令 $u = g(x)$，则 $du = g'(x)dx$，积分变为 $\int f(u) du$。
    *   **适用场景**: 被积函数中，一部分是另一部分的导数（或只差一个常数倍）。
    *   **示例**:
        *   $\int 2x e^{x^2} dx$。
        * 观察到 $(x^2)' = 2x$，令 $u=x^2$，$du=2xdx$。
        * 积分变为 $\int e^u du = e^u + C = e^{x^2} + C$。
        
        *   $\int \tan(x) dx = \int \frac{\sin(x)}{\cos(x)} dx$。
        * 观察到 $(\cos(x))' = -\sin(x)$，令 $u=\cos(x)$，$du=-sin(x)dx$。
        * 积分变为 $\int \frac{-1}{u} du = -\ln|u| + C = -\ln|\cos(x)| + C = \ln|\sec(x)| + C$。







*   **第二类换元法**:
    *   **核心思想**: 当被积函数形式复杂，特别是含有根式时，通过一个巧妙的变量代换 $x = g(t)$ 来简化积分。
    *   **适用场景**:
        *   **三角换元**: 处理含有 $\sqrt{a^2 - x^2}$, $\sqrt{a^2 + x^2}$, $\sqrt{x^2 - a^2}$ 的积分。
            *   令 $x = a\sin(t)$ 来消除 $\sqrt{a^2 - x^2}$。
            *   令 $x = a\tan(t)$ 来消除 $\sqrt{a^2 + x^2}$。
            *   令 $x = a\sec(t)$ 来消除 $\sqrt{x^2 - a^2}$。
        *   **根式换元**: 处理含有 $\sqrt[n]{ax+b}$ 等形式的积分。
            *   令 $t = \sqrt[n]{ax+b}$，反解出 $x$。

**类型一：处理 $\sqrt{a^2 - x^2}$ -> 令 $x = a \sin(t)$**

**示例：计算不定积分 $\int \frac{dx}{x^2 \sqrt{9 - x^2}}$**

1.  **识别类型并进行代换**:
    *   被积函数中含有 $\sqrt{9 - x^2} = \sqrt{3^2 - x^2}$，所以 $a = 3$。
    *   令 $x = 3 \sin(t)$。
    *   为了能换元回来，我们限制 $t$ 的范围，通常取 $t \in [-\pi/2, \pi/2]$。在这个范围内，$cos(t) \ge 0$。
    *   求 $dx$: $dx = 3 \cos(t) dt$。

2.  **替换积分中的每一部分**:
    *   $x^2 = (3 \sin(t))^2 = 9 \sin^2(t)$
    *   $\sqrt{9 - x^2} = \sqrt{9 - 9 \sin^2(t)} = \sqrt{9(1 - \sin^2(t))} = \sqrt{9 \cos^2(t)} = 3|\cos(t)|$。
    *   因为我们限制了 $t \in [-\pi/2, \pi/2]$，所以 $\cos(t) \ge 0$，因此 $3|\cos(t)| = 3 \cos(t)$。
    *   $dx = 3 \cos(t) dt$

3.  **代入并化简积分**:
    $\int \frac{dx}{x^2 \sqrt{9 - x^2}} = \int \frac{3 \cos(t) dt}{9 \sin^2(t) \cdot 3 \cos(t)}$
    *   $3 \cos(t)$ 在分子分母中可以约掉。
    $= \int \frac{dt}{9 \sin^2(t)} = \frac{1}{9} \int \csc^2(t) dt$

4.  **计算三角函数积分**:
    *   我们知道 $\int \csc^2(t) dt = -\cot(t) + C$。
    *   所以，$\frac{1}{9} \int \csc^2(t) dt = - \frac{1}{9} \cot(t) + C$。

5.  **将结果换回原始变量 x**:
    *   我们需要根据最初的代换 $x = 3 \sin(t)$ 来表达 $\cot(t)$。
    *   从 $x = 3 \sin(t)$ 可得 $\sin(t) = x/3$。
    *   我们可以画一个直角三角形来帮助理解：
        *   角 $t$ 的对边是 $x$，斜边是 $3$。
        *   根据勾股定理，邻边是 $\sqrt{3^2 - x^2} = \sqrt{9 - x^2}$。
    *   $\cot(t) = \frac{\text{邻边}}{\text{对边}} = \frac{\sqrt{9 - x^2}}{x}$。
    *   将此代回结果：$- \frac{1}{9} \cdot \frac{\sqrt{9 - x^2}}{x} + C$

**最终答案**: $-\frac{\sqrt{9 - x^2}}{9x} + C$

---

**类型二：处理 $\sqrt{a^2 + x^2}$ -> 令 $x = a \tan(t)$**

**利用的恒等式**: $1 + \tan^2(t) = \sec^2(t)$

**示例：计算不定积分 $\int \frac{dx}{\sqrt{4 + x^2}}$**

1.  **识别类型并进行代换**:
    *   被积函数中含有 $\sqrt{4 + x^2} = \sqrt{2^2 + x^2}$，所以 $a = 2$。
    *   令 $x = 2 \tan(t)$。
    *   限制 $t \in (-\pi/2, \pi/2)$，在这个范围内 $\sec(t) > 0$。
    *   求 $dx$: $dx = 2 \sec^2(t) dt$。

2.  **替换积分中的每一部分**:
    *   $\sqrt{4 + x^2} = \sqrt{4 + 4 \tan^2(t)} = \sqrt{4(1 + \tan^2(t))} = \sqrt{4 \sec^2(t)} = 2|\sec(t)|$。
    *   因为 $t \in (-\pi/2, \pi/2)$，所以 $\sec(t) > 0$，因此 $2|\sec(t)| = 2 \sec(t)$。
    *   $dx = 2 \sec^2(t) dt$

3.  **代入并化简积分**:
    $\int \frac{dx}{\sqrt{4 + x^2}} = \int \frac{2 \sec^2(t) dt}{2 \sec(t)}$
    *   约掉 $2 \sec(t)$。
    $= \int \sec(t) dt$

4.  **计算三角函数积分**:
    *   这是一个标准积分：$\int \sec(t) dt = \ln|\sec(t) + \tan(t)| + C$。

5.  **将结果换回原始变量 x**:
    *   从 $x = 2 \tan(t)$ 可得 $\tan(t) = x/2$。
    *   画一个直角三角形：
        *   角 $t$ 的对边是 $x$，邻边是 $2$。
        *   斜边是 $\sqrt{x^2 + 2^2} = \sqrt{x^2 + 4}$。
    *   $\sec(t) = \frac{\text{斜边}}{\text{邻边}} = \frac{\sqrt{x^2 + 4}}{2}$。
    *   将 $\tan(t)$ 和 $\sec(t)$ 代回结果：
        $\ln\left| \frac{\sqrt{x^2 + 4}}{2} + \frac{x}{2} \right| + C$
        $= \ln\left| \frac{x + \sqrt{x^2 + 4}}{2} \right| + C$
        $= \ln|x + \sqrt{x^2 + 4}| - \ln(2) + C$。
    *   因为 $-\ln(2)$ 是一个常数，可以被吸收到积分常数 $C$ 中。

**最终答案**: $\ln|x + \sqrt{x^2 + 4}| + C$ (这也是 $\text{arsinh}(x/2)$ 的表达式)

---

**类型三：处理 $\sqrt{x^2 - a^2}$ -> 令 $x = a \sec(t)$**

**利用的恒等式**: $\sec^2(t) - 1 = \tan^2(t)$

**示例：计算不定积分 $\int \frac{\sqrt{x^2 - 25}}{x} dx$ (假设 $x > 5$)**

1.  **识别类型并进行代换**:
    *   被积函数中含有 $\sqrt{x^2 - 25} = \sqrt{x^2 - 5^2}$，所以 $a = 5$。
    *   令 $x = 5 \sec(t)$。
    *   为了方便，限制 $t \in [0, \pi/2)$。在这个范围内，$\tan(t) \ge 0$ 且 $\sec(t) \ge 1$。
    *   求 $dx$: $dx = 5 \sec(t) \tan(t) dt$。

2.  **替换积分中的每一部分**:
    *   $\sqrt{x^2 - 25} = \sqrt{25 \sec^2(t) - 25} = \sqrt{25(\sec^2(t) - 1)} = \sqrt{25 \tan^2(t)} = 5|\tan(t)|$。
    *   因为 $t \in [0, \pi/2)$，所以 $\tan(t) \ge 0$，因此 $5|\tan(t)| = 5 \tan(t)$。
    *   $x = 5 \sec(t)$
    *   $dx = 5 \sec(t) \tan(t) dt$

3.  **代入并化简积分**:
    $\int \frac{\sqrt{x^2 - 25}}{x} dx = \int \frac{5 \tan(t)}{5 \sec(t)} \cdot (5 \sec(t) \tan(t) dt)$
    *   约掉 $5 \sec(t)$。
    $= \int 5 \tan^2(t) dt$

4.  **计算三角函数积分**:
    *   利用 $\tan^2(t) = \sec^2(t) - 1$。
    $= 5 \int (\sec^2(t) - 1) dt = 5 \left[\int \sec^2(t) dt - \int 1 dt\right]$
    $= 5 [\tan(t) - t] + C$

5.  **将结果换回原始变量 x**:
    *   从 $x = 5 \sec(t)$ 可得 $\sec(t) = x/5$。
    *   画一个直角三角形：
        *   $\sec(t) = \frac{\text{斜边}}{\text{邻边}}$，所以斜边是 $x$，邻边是 $5$。
        *   对边是 $\sqrt{x^2 - 5^2} = \sqrt{x^2 - 25}$。
    *   $\tan(t) = \frac{\text{对边}}{\text{邻边}} = \frac{\sqrt{x^2 - 25}}{5}$。
    *   $t$ 怎么表示？从 $\sec(t) = x/5$ 可得 $t = \text{arcsec}(x/5)$。
    *   代回结果：$5 \left[ \frac{\sqrt{x^2 - 25}}{5} - \text{arcsec}(x/5) \right] + C$
        $= \sqrt{x^2 - 25} - 5 \text{arcsec}(x/5) + C$

**最终答案**: $\sqrt{x^2 - 25} - 5 \text{arcsec}(x/5) + C$




> [!NOTE]
> #### 2. 分部积分法 

*   **核心思想**: 乘积求导法则 $(uv)' = u'v + uv'$ 的逆运算。
*   **公式**: $\int u \, dv = uv - \int v \, du$
*   **适用场景**: 被积函数是两种不同类型函数（如幂函数、指数函数、对数函数、三角函数）的乘积。
*   **选择 u 和 dv 的技巧 (LIATE 法则)**:
    按照 **L**ogarithmic (对数), **I**nverse Trig (反三角), **A**lgebraic (多项式/幂函数), **T**rigonometric (三角), **E**xponential (指数) 的顺序，排在前面的函数优先选作 `u`。
*   **示例**:
    *   $\int x e^x dx$。
    * 根据法则，$u = x$ (代数), $dv = e^x dx$ (指数)。
    * 则 $du = dx$, $v = e^x$。
    * $\int x e^x dx = xe^x - \int e^x dx = xe^x - e^x + C$。

    *   $\int \ln(x) dx$。
    * 看似不是乘积，但可以看作 $\int \ln(x) \cdot 1 dx$。
    * 令 $u = \ln(x)$, $dv = 1 dx$。
    * 则 $du = (1/x)dx$, $v = x$。
    * $\int \ln(x) dx = x\ln(x) - \int x \cdot \frac{1}{x} dx = x\ln(x) - \int 1 dx = x\ln(x) - x + C$。




> [!NOTE]
> #### 3. **有理函数积分**

*   **核心思想**: 任何有理函数（两个多项式的商）都可以通过==多项式长除法和部分分式分解==，转化为基本函数和简单分式的积分。

##### **例1：分子次数 ≥ 分母次数 (需要长除法)**

**问题：求解不定积分 $\int \frac{x^3 + x}{x - 1} dx$**

1.  **化为真分式**:
    *   分子次数 (3) 大于分母次数 (1)，所以需要用多项式长除法。
    ```
          x² + x + 2
        ____________
    x-1 | x³ + 0x² + x + 0
         -(x³ - x²)
        _________
              x² + x
             -(x² - x)
            _________
                   2x + 0
                  -(2x - 2)
                 ________
                        2 
    ```
    *   除法结果是：商式为 $x^2 + x + 2$，余式为 $2$。
    *   所以，$\frac{x^3 + x}{x - 1} = x^2 + x + 2 + \frac{2}{x - 1}$。

2.  **积分**:
    *   $\int \frac{x^3 + x}{x - 1} dx = \int \left(x^2 + x + 2 + \frac{2}{x - 1}\right) dx$
    *   利用线性法则拆分：
        $= \int x^2 dx + \int x dx + \int 2 dx + \int \frac{2}{x - 1} dx$
    *   分别积分：
        $= \frac{x^3}{3} + \frac{x^2}{2} + 2x + 2 \ln|x - 1| + C$

**最终答案**: $\frac{x^3}{3} + \frac{x^2}{2} + 2x + 2 \ln|x - 1| + C$

---

##### **例2：分母含不重复的一次因式和不可约二次因式**

**问题：求解不定积分 $\int \frac{1}{x(x^2 + 4)} dx$**

1.  **部分分式分解**:
    *   $\frac{1}{x(x^2 + 4)} = \frac{A}{x} + \frac{Bx + C}{x^2 + 4}$
    *   两边同乘以 $x(x^2 + 4)$，得到恒等式：
        $1 = A(x^2 + 4) + (Bx + C)x$
    *   展开右边：$1 = Ax^2 + 4A + Bx^2 + Cx$
    *   合并同类项：$1 = (A + B)x^2 + Cx + 4A$
    *   比较系数：
        *   $x^2$ 项系数: $A + B = 0$
        *   $x$ 项系数: $C = 0$
        *   常数项: $4A = 1 \implies A = 1/4$
    *   从 $A + B = 0$ 可得 $B = -A = -1/4$。
    *   展开式为：$\frac{1/4}{x} + \frac{-1/4 \cdot x}{x^2 + 4}$

2.  **积分**:
    *   $\int \frac{1}{x(x^2 + 4)} dx = \int \left( \frac{1/4}{x} - \frac{x}{4(x^2 + 4)} \right) dx$
    *   $= \frac{1}{4} \int \frac{1}{x} dx - \frac{1}{4} \int \frac{x}{x^2 + 4} dx$
    *   第一个积分很简单：$\frac{1}{4} \ln|x|$。
    *   第二个积分 $\int \frac{x}{x^2 + 4} dx$，使用凑微分法：
        *   令 $u = x^2 + 4$，则 $du = 2x dx$，所以 $x dx = du/2$。
        *   $\int \frac{1}{u} \frac{du}{2} = \frac{1}{2} \int \frac{1}{u} du = \frac{1}{2} \ln|u| = \frac{1}{2} \ln(x^2 + 4)$。
    *   将结果合并：
        $= \frac{1}{4} \ln|x| - \frac{1}{4} \cdot \frac{1}{2} \ln(x^2 + 4) + C$
        $= \frac{1}{4} \ln|x| - \frac{1}{8} \ln(x^2 + 4) + C$

**最终答案**: $\frac{1}{4} \ln|x| - \frac{1}{8} \ln(x^2 + 4) + C$

---

##### **例3：分母含重复的一次因式**

**问题：求解不定积分 $\int \frac{x + 4}{(x + 2)^2} dx$**

1.  **部分分式展开**:
    *   $\frac{x + 4}{(x + 2)^2} = \frac{A}{x + 2} + \frac{B}{(x + 2)^2}$
    *   两边同乘以 $(x + 2)^2$：
        $x + 4 = A(x + 2) + B$
    *   令 $x = -2$：$-2 + 4 = A(0) + B \implies B = 2$
    *   比较 $x$ 的系数，直接得到 $A = 1$。
    *   展开式为：$\frac{1}{x + 2} + \frac{2}{(x + 2)^2}$

2.  **积分**:
    *   $\int \frac{x + 4}{(x + 2)^2} dx = \int \frac{1}{x + 2} dx + \int \frac{2}{(x + 2)^2} dx$
    *   $= \ln|x + 2| + 2 \int (x + 2)^{-2} dx$
    *   $= \ln|x + 2| + 2 \cdot \frac{(x + 2)^{-1}}{-1} + C$
    *   $= \ln|x + 2| - \frac{2}{x + 2} + C$

**最终答案**: $\ln|x + 2| - \frac{2}{x + 2} + C$








**4. 以曲线 $y = e^x + x$ ($x \in [0, 1]$) 为曲边绕 x 轴旋转一周所得立体的体积 V = __$\pi(\frac{e^2}{2} + \frac{11}{6})$__**

> [!NOTE]
> *   **思路:** 使用圆盘法求旋转体体积，公式为 $V = \pi \int_{a}^{b} y^2 dx$。

*   **步骤:**
    1.  $V = \pi \int_{0}^{1} (e^x + x)^2 dx$
    2.  展开被积函数: $(e^x + x)^2 = (e^x)^2 + 2xe^x + x^2 = e^{2x} + 2xe^x + x^2$
    3.  分部积分 $\int 2xe^x dx$：令 $u=2x$, $dv=e^x dx$，则 $du=2dx$, $v=e^x$。
        $\int 2xe^x dx = 2xe^x - \int 2e^x dx = 2xe^x - 2e^x$
    4.  计算定积分：
        $V = \pi \left[ \int_{0}^{1} e^{2x} dx + \int_{0}^{1} 2xe^x dx + \int_{0}^{1} x^2 dx \right]$
        $V = \pi \left[ \left[\frac{1}{2}e^{2x}\right]_{0}^{1} + [2xe^x - 2e^x]_{0}^{1} + \left[\frac{1}{3}x^3\right]_{0}^{1} \right]$
        $V = \pi \left[ \left(\frac{1}{2}e^2 - \frac{1}{2}e^0\right) + \left((2e - 2e) - (0 - 2e^0)\right) + \left(\frac{1}{3} - 0\right) \right]$
        $V = \pi \left[ \left(\frac{e^2}{2} - \frac{1}{2}\right) + (0 - (-2)) + \frac{1}{3} \right]$
        $V = \pi \left[ \frac{e^2}{2} - \frac{1}{2} + 2 + \frac{1}{3} \right]$
        $V = \pi \left[ \frac{e^2}{2} + \frac{3}{2} + \frac{1}{3} \right] = \pi \left[ \frac{e^2}{2} + \frac{9}{6} + \frac{2}{6} \right] = \pi\left(\frac{e^2}{2} + \frac{11}{6}\right)$


---

##### **1. 圆盘法** 

*   **适用场景**:
    *   旋转区域**紧贴**着旋转轴，中间没有空隙。
    *   积分微元（一个薄片）垂直于旋转轴。
*   **公式**:
    *   **绕 x 轴旋转**: $V = \int_{a}^{b} \pi [f(x)]^2 dx$
    *   **绕 y 轴旋转**: $V = \int_{c}^{d} \pi [g(y)]^2 dy$ (曲线表示为 $x = g(y)$)

---
##### **2. 圆环法** 

*   **适用场景**:
    *   旋转区域**不紧贴**旋转轴，被两条曲线所包围。
    *   积分微元（一个薄片）仍然垂直于旋转轴。
*   **公式**:
    *   **绕 x 轴旋转** (区域由 $y=f(x)$ 和 $y=g(x)$ 包围，$f(x) \ge g(x)$)：
        $V = \int_{a}^{b} \pi [ (f(x))^2 - (g(x))^2 ] dx$
    *   **绕 y 轴旋转** (区域由 $x=f(y)$ 和 $x=g(y)$ 包围，$f(y) \ge g(y)$)：
        $V = \int_{c}^{d} \pi [ (f(y))^2 - (g(y))^2 ] dy$

**示例**:
求由 $y = x$ 和 $y = x^2$ 围成的区域绕 x 轴旋转一周所得的体积。
*   交点是 (0,0) 和 (1,1)。在 $[0,1]$ 上, $x \ge x^2$。
*   外半径 $R(x) = x$。
*   内半径 $r(x) = x^2$。
*   $V = \int_{0}^{1} \pi [ (x)^2 - (x^2)^2 ] dx = \pi \int_{0}^{1} (x^2 - x^4) dx = \pi \left[\frac{x^3}{3} - \frac{x^5}{5}\right]_{0}^{1} = \pi\left(\frac{1}{3} - \frac{1}{5}\right) = \frac{2\pi}{15}$。

---








**5. 若 $f(x) = \sqrt{1 - x^2} + x^3 \int_{0}^{1} f(x) dx$, 则 $\int_{0}^{1} f(x) dx$ = __$\pi/3$__**

> [!NOTE]
> *   **思路:** $\int_{0}^{1} f(x) dx$ 是一个常数，可以设其为 C，然后解出 C。

*   **步骤:**
    1.  令 $C = \int_{0}^{1} f(x) dx$。
    2.  则 $f(x) = \sqrt{1 - x^2} + C x^3$。
    3.  将 f(x) 的表达式代入 C 的定义中：
        $C = \int_{0}^{1} (\sqrt{1 - x^2} + C x^3) dx$
        $C = \int_{0}^{1} \sqrt{1 - x^2} dx + \int_{0}^{1} C x^3 dx$
    4.  计算两个积分：
        $\int_{0}^{1} \sqrt{1 - x^2} dx$ 表示半径为1的圆在第一象限的面积，即 $\frac{1}{4}\pi(1)^2 = \frac{\pi}{4}$。
        $\int_{0}^{1} C x^3 dx = C \left[\frac{x^4}{4}\right]_{0}^{1} = C \left(\frac{1}{4}\right)$
    5.  代回方程解 C：
        $C = \frac{\pi}{4} + \frac{C}{4}$
        $C - \frac{C}{4} = \frac{\pi}{4}$
        $\frac{3}{4}C = \frac{\pi}{4}$
        $C = \frac{\pi}{3}$








**6. 星形线 $\{ x = \cos^3 t, y = \sin^3 t \}$ 在第一象限 ($0 \le t \le \pi/2$) 的弧长为 __3/2__**

> [!NOTE]
> *   **思路:** 使用参数方程的弧长公式
> $$
> L = \int_{a}^{b} \sqrt{\left(\frac{dx}{dt}\right)^2 + \left(\frac{dy}{dt}\right)^2} \, dt
> $$

*   **步骤:**
    1.  求导数：
        $\frac{dx}{dt} = 3\cos^2 t \cdot (-\sin t) = -3\cos^2 t \sin t$
        $\frac{dy}{dt} = 3\sin^2 t \cdot (\cos t)$
    2.  计算 $(\frac{dx}{dt})^2 + (\frac{dy}{dt})^2$：
        $(-3\cos^2 t \sin t)^2 + (3\sin^2 t \cos t)^2$
        $= 9\cos^4 t \sin^2 t + 9\sin^4 t \cos^2 t$
        $= 9\sin^2 t \cos^2 t (\cos^2 t + \sin^2 t)$
        $= 9\sin^2 t \cos^2 t$
    3.  开方：$\sqrt{(\frac{dx}{dt})^2 + (\frac{dy}{dt})^2} = \sqrt{9\sin^2 t \cos^2 t} = 3|\sin t \cos t|$。
        在第一象限 $0 \le t \le \pi/2$，$\sin t$ 和 $\cos t$ 均为非负，所以 $|\sin t \cos t| = \sin t \cos t$。
    4.  计算积分：
        $L = \int_{0}^{\pi/2} 3\sin t \cos t dt$
        利用 $\sin(2t) = 2\sin t \cos t$：
        $L = \frac{3}{2} \int_{0}^{\pi/2} \sin(2t) dt$
        $L = \frac{3}{2} \left[-\frac{\cos(2t)}{2}\right]_{0}^{\pi/2}$
        $L = \frac{3}{2} \left[ \left(-\frac{\cos(\pi)}{2}\right) - \left(-\frac{\cos(0)}{2}\right) \right]$
        $L = \frac{3}{2} \left[ \left(-\frac{-1}{2}\right) - \left(-\frac{1}{2}\right) \right]$
        $L = \frac{3}{2} \left[ \frac{1}{2} + \frac{1}{2} \right] = \frac{3}{2} \cdot 1 = \frac{3}{2}$

---









### **三、计算题**

**1. 已知 $\int_{0}^{x} f(t)dt = xf(ux)$，且 $f(x) = e^x$，求极限 $\lim_{x \to 0^+} u$ .**

**解题思路:**
首先将已知的 $f(x) = e^x$ 代入给定的方程中，得到一个包含 $x$ 和 $u$ 的关系式。然后从这个关系式中解出 $u$，最后对 $u$ 求极限。

**详细步骤:**

1.  **代入 f(x):**
    将 $f(x) = e^x$ 代入方程 $\int_{0}^{x} f(t)dt = xf(ux)$。
    *   左边：$\int_{0}^{x} e^t dt = [e^t]_{0}^{x} = e^x - e^0 = e^x - 1$
    *   右边：$xf(ux) = x e^{ux}$
    *   所以，我们得到方程：$e^x - 1 = x e^{ux}$

2.  **解出 u:**
    *   $e^{ux} = \frac{e^x - 1}{x}$
    *   对两边取自然对数：
        $ux = \ln\left(\frac{e^x - 1}{x}\right)$
    *   解得 $u$：
        $u = \frac{1}{x} \ln\left(\frac{e^x - 1}{x}\right)$

3.  **求极限:**
    $\lim_{x \to 0^+} u = \lim_{x \to 0^+} \frac{\ln\left(\frac{e^x - 1}{x}\right)}{x}$
    *   这是一个 "0/0" 型的极限，因为当 $x \to 0^+$ 时：
        *   分子 $\ln\left(\frac{e^x - 1}{x}\right) \to \ln(1) \to 0$ (使用了重要极限 $\lim_{x \to 0} \frac{e^x - 1}{x} = 1$)
        *   分母 $x \to 0$


#### **方法一：洛必达法则

1.  **求分母的导数:**
    $$
    \frac{d}{dx}(x) = 1
    $$

2.  **求分子的导数 (这是关键步骤):**
    我们需要对 $f(x) = \ln\left(\frac{e^x - 1}{x}\right)$ 求导。

    为了计算方便，我们可以先利用对数性质将函数拆分：
    $$
    f(x) = \ln(e^x - 1) - \ln(x)
    $$
    现在分别对这两项求导：
    *   $\frac{d}{dx}[\ln(e^x - 1)] = \frac{1}{e^x - 1} \cdot \frac{d}{dx}(e^x - 1) = \frac{e^x}{e^x - 1}$
    *   $\frac{d}{dx}[\ln(x)] = \frac{1}{x}$

    所以，分子的导数是：
    $$
    f'(x) = \frac{e^x}{e^x - 1} - \frac{1}{x}
    $$

3.  **应用法则并计算新极限:**
    根据洛必达法则，原极限 $L$ 等于：
    $$
    L = \lim_{x \to 0^+} \frac{\frac{e^x}{e^x - 1} - \frac{1}{x}}{1} = \lim_{x \to 0^+} \left(\frac{e^x}{e^x - 1} - \frac{1}{x}\right)
    $$

**第二次洛必达
4.  **通分:**
    $$
    L = \lim_{x \to 0^+} \frac{x e^x - (e^x - 1)}{x(e^x - 1)} = \lim_{x \to 0^+} \frac{x e^x - e^x + 1}{x e^x - x}
    $$

*   **新分子的导数:**
	$\frac{d}{dx}(x e^x - e^x + 1) = (1 \cdot e^x + x \cdot e^x) - e^x + 0 = e^x + xe^x - e^x = xe^x$
*   **新分母的导数:**
	$\frac{d}{dx}(x e^x - x) = (1 \cdot e^x + x \cdot e^x) - 1 = e^x + xe^x - 1$

5.  **计算最终的极限:**
    $$
    L = \lim_{x \to 0^+} \frac{xe^x}{e^x + xe^x - 1}
    $$

**第三次洛必达**

6.  **第三次对分子和分母求导:**
*   **分子的导数:**
	$\frac{d}{dx}(xe^x) = 1 \cdot e^x + x \cdot e^x = e^x + xe^x$
*   **分母的导数:**
	$\frac{d}{dx}(e^x + xe^x - 1) = e^x + (e^x + xe^x) - 0 = 2e^x + xe^x$

7.  **计算极限:**
    $$
    L = \lim_{x \to 0^+} \frac{e^x + xe^x}{2e^x + xe^x}
    $$
    现在，这个极限可以直接通过代入 $x=0$ 来计算了，因为它不再是不定型。
    $$
    L = \frac{e^0 + 0 \cdot e^0}{2e^0 + 0 \cdot e^0} = \frac{1 + 0}{2 + 0} = \frac{1}{2}
    $$


#### **方法二：泰勒展开**
*   当 $x \to 0$ 时，$e^x \approx 1 + x + \frac{x^2}{2}$。
*   所以 $\frac{e^x - 1}{x} \approx \frac{x + x^2/2}{x} = 1 + \frac{x}{2}$。
*   代入 $u$ 的表达式中：$u = \frac{\ln\left(1 + \frac{x}{2}\right)}{x}$
*   当 $y \to 0$ 时，$\ln(1+y) \approx y$。令 $y = x/2$：
	$u \approx \frac{1}{x} \cdot \frac{x}{2} = \frac{1}{2}$
*   因此，$\lim_{x \to 0^+} u = \frac{1}{2}$。



#### **方法三：导数定义

**核心思想：** 将极限表达式凑成 $f'(a) = \lim_{x \to a} \frac{f(x) - f(a)}{x - a}$ 的形式。

1.  **定义一个函数**
    让我们定义一个函数 $g(x)$ 如下：
    $$
    g(x) = \ln\left(\frac{e^x - 1}{x}\right)
    $$
    这个函数在 $x=0$ 处没有定义。但是，我们可以计算它在 $x \to 0$ 时的极限：
    $$
    \lim_{x \to 0} g(x) = \lim_{x \to 0} \ln\left(\frac{e^x - 1}{x}\right) = \ln\left(\lim_{x \to 0} \frac{e^x - 1}{x}\right) = \ln(1) = 0
    $$
    所以，我们可以定义一个在 $x=0$ 处连续的新函数 $f(x)$：
    $$
    f(x) =
    \begin{cases}
    \ln\left(\frac{e^x - 1}{x}\right) & \text{if } x \neq 0 \\
    0 & \text{if } x = 0
    \end{cases}
    $$

2.  **将原极限改写为导数形式**
    现在，我们要求的极限可以被重写为：
    $$
    L = \lim_{x \to 0^+} \frac{f(x)}{x} = \lim_{x \to 0^+} \frac{f(x) - 0}{x - 0} = \lim_{x \to 0^+} \frac{f(x) - f(0)}{x - 0}
    $$
    这正是函数 $f(x)$ 在 $x=0$ 处的**右导数**的定义，即 $f'_+(0)$。

3.  **计算导数**
    根据求导法则，如果 $f'(x)$ 在 $x=0$ 处连续，那么 $f'(0) = \lim_{x \to 0} f'(x)$。我们来计算 $f'(x)$ (当 $x \neq 0$ 时)：
    $$
    f'(x) = \frac{e^x}{e^x - 1} - \frac{1}{x}
    $$
    通分得到：
    $$
    f'(x) = \frac{x e^x - (e^x - 1)}{x(e^x - 1)} = \frac{x e^x - e^x + 1}{x e^x - x}
    $$

4.  **求解导数的极限**
    *   **使用洛必达法则 (对 $f'(x)$):**
        $$
        \lim_{x \to 0} \frac{x e^x - e^x + 1}{x e^x - x} \overset{L'H}{=} \lim_{x \to 0} \frac{(e^x + xe^x) - e^x}{(e^x + xe^x) - 1} = \lim_{x \to 0} \frac{xe^x}{e^x + xe^x - 1}
        $$
        这还是 "0/0" 型，再用一次：
        $$
        \overset{L'H}{=} \lim_{x \to 0} \frac{e^x + xe^x}{e^x + (e^x + xe^x)} = \frac{e^0 + 0}{e^0 + e^0 + 0} = \frac{1}{2}
        $$
    
    所以，我们得到 $f'(0) = 1/2$。

**结论：** 原极限等于 $f'(0)$，也就是 **1/2**。



#### **方法四：利用夹逼定理

这种方法非常严谨，但需要构造合适的不等式，技巧性较强。

**核心思想：** 找到两个函数 $g(x)$ 和 $h(x)$，使得 $g(x) \le f(x) \le h(x)$，并且 $\lim_{x \to 0} g(x) = \lim_{x \to 0} h(x) = L$。

1.  **构造不等式**
    我们需要一个关于 $\ln(1+u)$ 的不等式。对于 $u > 0$，有一个著名且有用的不等式：
    $$
    \frac{u}{1+u} \le \ln(1+u) \le u
    $$

2.  **变量代换**
    我们的表达式是 $\frac{1}{x} \ln\left(\frac{e^x - 1}{x}\right)$。
    我们令 $u = \frac{e^x - 1}{x} - 1$。
    当 $x \to 0^+$ 时，$u = (\frac{1+x+x^2/2+... - 1}{x}) - 1 = (1+x/2+...) - 1 = x/2+... > 0$。
    所以我们可以使用上面的不等式。
    原表达式中的 $\ln\left(\frac{e^x-1}{x}\right)$ 可以写成 $\ln(1+u)$。

3.  **应用不等式**
    将 $\ln(1+u)$ 替换为它的上下界，并同除以 $x$ ($x > 0$，不等号方向不变)：
    $$
    \frac{1}{x} \left( \frac{u}{1+u} \right) \le \frac{\ln(1+u)}{x} \le \frac{u}{x}
    $$
    代回 $u$ 的表达式，我们得到一个下界和一个上界。

4.  **计算上界的极限**
    $$
    \lim_{x \to 0^+} \frac{u}{x} = \lim_{x \to 0^+} \frac{\frac{e^x - 1}{x} - 1}{x} = \lim_{x \to 0^+} \frac{e^x - 1 - x}{x^2}
    $$
    这是一个 "0/0" 型，用两次洛必达法则：
    $$
    \overset{L'H}{=} \lim_{x \to 0^+} \frac{e^x - 1}{2x} \overset{L'H}{=} \lim_{x \to 0^+} \frac{e^x}{2} = \frac{1}{2}
    $$
    所以，我们找到了上界的极限是 **1/2**。

5.  **计算下界的极限**
    $$
    \lim_{x \to 0^+} \frac{1}{x} \left( \frac{u}{1+u} \right) = \left( \lim_{x \to 0^+} \frac{u}{x} \right) \cdot \left( \lim_{x \to 0^+} \frac{1}{1+u} \right)
    $$
    我们已经知道 $\lim_{x \to 0^+} \frac{u}{x} = \frac{1}{2}$。
    并且当 $x \to 0^+$ 时，$u \to 0$，所以 $\lim_{x \to 0^+} \frac{1}{1+u} = \frac{1}{1+0} = 1$。
    因此，下界的极限是 $\frac{1}{2} \cdot 1 = \frac{1}{2}$。

6.  **得出结论**
    因为下界和上界的极限都等于 **1/2**，根据夹逼定理，原极限也必须等于 **1/2**。


| 方法        | 核心思想                                          | 优点        | 缺点                    |
| :-------- | :-------------------------------------------- | :-------- | :-------------------- |
| **洛必达法则** | 处理不定型 $\frac{0}{0}$ 或 $\frac{\infty}{\infty}$ | 机械化，直接    | 可能需要多次求导，计算繁琐         |
| **泰勒展开**  | 用多项式逼近函数                                      | 功能强大，揭示本质 | 需要记住常用展开式，代数运算易错      |
| **导数定义**  | 将极限看作某函数的导数                                   | 思路巧妙，形式优雅 | 需要构造合适的函数，可能最终还是要用洛必达 |
| **夹逼定理**  | 用两个极限相同的函数夹住目标                                | 逻辑严谨，非常根本 | 技巧性强，找到合适的不等式是难点      |

---













**2. 设 $f(x)$ 的一个原函数为 $\ln(x + \sqrt{1+x^2})$，求 $\int xf'(x)dx$.**

**解题思路:**
首先根据原函数的定义求出 $f(x)$，然后使用==分部积分法==计算所求的不定积分。

**详细步骤:**

1.  **求 f(x):**
    设 $F(x) = \ln(x + \sqrt{1+x^2})$。根据定义，$f(x) = F'(x)$。
    $f(x) = \frac{d}{dx} \left[\ln(x + \sqrt{1+x^2})\right]$
    $f(x) = \frac{1}{x + \sqrt{1+x^2}} \cdot \frac{d}{dx}\left[x + \sqrt{1+x^2}\right]$
    $f(x) = \frac{1}{x + \sqrt{1+x^2}} \cdot \left[1 + \frac{1}{2\sqrt{1+x^2}} \cdot 2x\right]$
    $f(x) = \frac{1}{x + \sqrt{1+x^2}} \cdot \left[1 + \frac{x}{\sqrt{1+x^2}}\right]$
    $f(x) = \frac{1}{x + \sqrt{1+x^2}} \cdot \left[\frac{\sqrt{1+x^2} + x}{\sqrt{1+x^2}}\right]$
    $f(x) = \frac{1}{\sqrt{1+x^2}}$

2.  **使用分部积分法:**
    计算 $\int xf'(x)dx$。
    *   令 $u = x$， $dv = f'(x)dx$。
    *   则 $du = dx$， $v = \int f'(x)dx = f(x)$。
    *   根据分部积分公式 $\int u \, dv = uv - \int v \, du$：
        $\int xf'(x)dx = x f(x) - \int f(x) dx$
    *   我们知道 $f(x) = \frac{1}{\sqrt{1+x^2}}$。
    *   我们还知道 $\int f(x) dx$ 就是 $f(x)$ 的原函数，题目已给出为 $F(x) = \ln(x + \sqrt{1+x^2})$。
    *   将这些代入公式：
        $\int xf'(x)dx = x \cdot \frac{1}{\sqrt{1+x^2}} - \ln(x + \sqrt{1+x^2}) + C$



**尝试使用第一类换元法

1.  **选择换元变量 `u`:**
    一个自然的选择是令 $u = f(x)$。
    那么 $du = f'(x)dx$。


2.  **进行换元:**
    $$
    \int x \underbrace{f'(x)dx}_{du}
    $$
    我们用 $du$ 替换了 $f'(x)dx$，但积分式里还剩下一个 $x$。我们必须把 $x$ 也用 $u$ 来表示。

3.  **反解出 `x`:**
    $$
    \sqrt{1+x^2} = \frac{1}{u}
    $$
    $$
    1 + x^2 = \frac{1}{u^2}
    $$
    $$
    x^2 = \frac{1}{u^2} - 1 = \frac{1-u^2}{u^2}
    $$
    $$
    x = \frac{\sqrt{1-u^2}}{u} \quad (\text{假设 x > 0})
    $$

4.  **得到新的积分:**
    将 $x = \frac{\sqrt{1-u^2}}{u}$ 代入换元后的积分式，我们得到一个全新的、只关于 $u$ 的积分：
    $$
    \int \frac{\sqrt{1-u^2}}{u} du
    $$

5.  **解决这个新积分:**
    这个新积分比原来的问题复杂多了！它需要用**第二次换元**，也就是**三角换元法**来解决。
    *   令 $u = \sin\theta$，则 $du = \cos\theta \, d\theta$。
    *   积分变为 $\int \frac{\sqrt{1-\sin^2\theta}}{\sin\theta} \cos\theta \, d\theta = \int \frac{\cos\theta}{\sin\theta} \cos\theta \, d\theta = \int \frac{\cos^2\theta}{\sin\theta} d\theta$
    *   $= \int \frac{1-\sin^2\theta}{\sin\theta} d\theta = \int (\csc\theta - \sin\theta) d\theta$
    *   $= \ln|\csc\theta - \cot\theta| + \cos\theta + C$

6.  **换元回来:**
    这还没完，我们得把 $\theta$ 换回 $u$，再把 $u$ 换回 $x$。
    *   从 $u = \sin\theta$ 得到 $\cos\theta = \sqrt{1-u^2}$，$\csc\theta=1/u$，$\cot\theta = \sqrt{1-u^2}/u$。
    *   代入得到：$\ln\left|\frac{1-\sqrt{1-u^2}}{u}\right| + \sqrt{1-u^2} + C$
    *   再代入 $u = \frac{1}{\sqrt{1+x^2}}$ ... 最终也能得到正确答案，但过程极其复杂。

对于 $\int x f'(x) dx$ 这种形式的积分，**分部积分法**几乎永远是首选的、最高效的方法。

**最终答案:** $\frac{x}{\sqrt{1+x^2}} - \ln(x + \sqrt{1+x^2}) + C$

---













~~**3. 计算由函数 $f(x) = x\cos^8 x$，直线 $x = 2\pi$ 以及两个坐标轴在第一象限所围图形面积。**~~

**解题思路:**
面积就是定积分 $A = \int_{0}^{2\pi} x\cos^8 x \, dx$。这个积分直接计算很困难，可以使用定积分的对称性性质来简化。

**详细步骤:**

1.  **建立面积积分:**
    $A = \int_{0}^{2\pi} x\cos^8 x \, dx$

2.  **使用积分性质:**
    $A = \int_{0}^{2\pi} (2\pi - x)\cos^8(2\pi - x) \, dx$
    因为 $\cos(2\pi - x) = \cos(x)$，所以 $\cos^8(2\pi - x) = \cos^8 x$。
    $A = \int_{0}^{2\pi} (2\pi - x)\cos^8 x \, dx = \int_{0}^{2\pi} 2\pi\cos^8 x \, dx - \int_{0}^{2\pi} x\cos^8 x \, dx$
    $A = 2\pi \int_{0}^{2\pi} \cos^8 x \, dx - A$
    $2A = 2\pi \int_{0}^{2\pi} \cos^8 x \, dx \implies A = \pi \int_{0}^{2\pi} \cos^8 x \, dx$

3.  **计算简化后的积分:**
    *   函数 $\cos^8 x$ 的周期是 $\pi$，所以 $\int_{0}^{2\pi} \cos^8 x \, dx = 2 \int_{0}^{\pi} \cos^8 x \, dx$。
    *   又因为 $\cos^8(\pi - x) = (-\cos x)^8 = \cos^8 x$，
    *   所以 $\int_{0}^{\pi} \cos^8 x \, dx = 2 \int_{0}^{\pi/2} \cos^8 x \, dx$。
    *   综上，$A = \pi \cdot (2 \cdot 2 \int_{0}^{\pi/2} \cos^8 x \, dx) = 4\pi \int_{0}^{\pi/2} \cos^8 x \, dx$。

4.  **使用瓦力斯 (Wallis) 积分公式:**
    对于偶数 n, $\int_{0}^{\pi/2} \cos^n x \, dx = \frac{(n-1)!!}{n!!} \frac{\pi}{2}$。
    *   这里 $n = 8$。
    *   $\int_{0}^{\pi/2} \cos^8 x \, dx = \frac{7 \cdot 5 \cdot 3 \cdot 1}{8 \cdot 6 \cdot 4 \cdot 2} \frac{\pi}{2} = \frac{105}{384} \frac{\pi}{2} = \frac{35}{128} \frac{\pi}{2} = \frac{35\pi}{256}$

5.  **计算最终面积 A:**
    $A = 4\pi \cdot \frac{35\pi}{256} = \frac{140\pi^2}{256} = \frac{35\pi^2}{64}$

**最终答案:** $\frac{35\pi^2}{64}$

---
















**4. 已知 $\int_{0}^{+\infty} \frac{\sin x}{x} dx = \frac{\pi}{2}$，求 $\int_{0}^{+\infty} \frac{\sin^2 x}{x^2} dx$.**

> [!NOTE]
> **解题思路:**
> 所求积分的形式很适合使用分部积分法，通过分部积分，可以将其转化为已知的狄利克雷积分 $\int_{0}^{+\infty} \frac{\sin x}{x} dx$。
> 

**详细步骤:**

1.  **设所求积分为 I:**
    $I = \int_{0}^{+\infty} \frac{\sin^2 x}{x^2} dx$

2.  **使用分部积分法:**
    令 $u = \sin^2 x$，$dv = \frac{1}{x^2} dx$。
    *   $du = 2\sin x \cos x \, dx = \sin(2x) \, dx$
    *   $v = \int \frac{1}{x^2} dx = -\frac{1}{x}$
    *   $I = \left[\sin^2 x \cdot \left(-\frac{1}{x}\right)\right]_{0}^{+\infty} - \int_{0}^{+\infty} \left(-\frac{1}{x}\right) \sin(2x) \, dx$

3.  **计算边界项:**
    *   $\lim_{x \to +\infty} -\frac{\sin^2 x}{x} = 0$ (夹逼定理)。
    *   $\lim_{x \to 0^+} -\frac{\sin^2 x}{x} = \lim_{x \to 0^+} -\left(\frac{\sin x}{x}\right)^2 x = -1^2 \cdot 0 = 0$。
    *   所以边界项为 $0 - 0 = 0$。

4.  **化简剩余的积分:**
    $I = \int_{0}^{+\infty} \frac{\sin(2x)}{x} dx$

5.  **变量代换:**
    *   令 $t = 2x$，则 $x = t/2$，$dx = dt/2$。
    *   $I = \int_{0}^{+\infty} \frac{\sin(t)}{t/2} \frac{dt}{2} = \int_{0}^{+\infty} \frac{\sin(t)}{t} dt$

6.  **利用已知条件:**
    $\int_{0}^{+\infty} \frac{\sin(t)}{t} dt = \frac{\pi}{2}$

**最终答案:** $\frac{\pi}{2}$

---



















**5. 求微分方程 $y'' - 3y' - 10y = 0$ 满足初始条件 $y(0) = 2$ 和 $y'(0) = 10$ 的特解。**

**解题思路:**
这是一个二阶常系数齐次线性微分方程。

**详细步骤:**

1.  **求解特征方程:**
    $r^2 - 3r - 10 = 0$
    $(r - 5)(r + 2) = 0$
    特征根为 $r_1 = 5$, $r_2 = -2$。

2.  **写出通解:**
    $y(x) = C_1 e^{5x} + C_2 e^{-2x}$

3.  **利用初始条件求特解:**
    $y'(x) = 5C_1 e^{5x} - 2C_2 e^{-2x}$
    *   `y(0) = 2`: $C_1 e^0 + C_2 e^0 = C_1 + C_2 = 2$
    *   `y'(0) = 10`: $5C_1 e^0 - 2C_2 e^0 = 5C_1 - 2C_2 = 10$
    解方程组:
    $C_1 + C_2 = 2 \implies C_2 = 2 - C_1$
    $5C_1 - 2(2 - C_1) = 10 \implies 5C_1 - 4 + 2C_1 = 10 \implies 7C_1 = 14 \implies C_1 = 2$
    $C_2 = 2 - 2 = 0$

4.  **写出特解:**
    $y(x) = 2e^{5x} + 0 \cdot e^{-2x} = 2e^{5x}$

**最终答案:** $y = 2e^{5x}$

---

### **四、(8 分) 设函数 $\phi(x)$ 连续，$f(x) = \int_{0}^{1} \phi(xt)dt$，且 $\lim_{x \to 0} \frac{\phi(x)}{x} = 2019$，求 $f(x)$ 和 $f'(x)$。**

**详细步骤:**

1.  **对 f(x) 进行变量代换:**
    令 $u = xt$，则 $dt = \frac{1}{x}du$。
    `f(x) = \int_{0}^{x} \phi(u) \frac{1}{x}du = \frac{1}{x} \int_{0}^{x} \phi(u)du`

2.  **利用极限条件确定 $\phi(x)$ 的形式:**
    $\lim_{x \to 0} \frac{\phi(x)}{x} = 2019$。因为极限存在且分母趋于0，所以分子 $\lim_{x \to 0} \phi(x) = 0$。又因为 $\phi(x)$ 连续，所以 $\phi(0) = 0$。
    该极限是导数的定义：$\phi'(0) = \lim_{x \to 0} \frac{\phi(x) - \phi(0)}{x - 0} = 2019$。
    最简单的满足此条件的连续函数形式是 $\phi(x) = 2019x$。

3.  **求 f(x):**
    将 $\phi(t) = 2019t$ 代入 $f(x) = \frac{1}{x} \int_{0}^{x} \phi(t)dt$：
    $f(x) = \frac{1}{x} \int_{0}^{x} 2019t \, dt = \frac{2019}{x} \left[\frac{t^2}{2}\right]_{0}^{x} = \frac{2019}{x} \frac{x^2}{2} = \frac{2019x}{2}$

4.  **求 f'(x):**
    $f'(x) = \frac{d}{dx} \left(\frac{2019x}{2}\right) = \frac{2019}{2}$

**最终答案:**
$f(x) = \frac{2019x}{2}$
$f'(x) = \frac{2019}{2}$

---

### **五、(9 分) 求微分方程 $x^2 y' + xy = y^2$ 满足初始条件 $y(1) = 1$ 的特解。**

**解题思路:**
这是一个**伯努利方程**。标准形式为 $y' + P(x)y = Q(x)y^n$。

**详细步骤:**

1.  **化为标准形式:**
    $y' + \frac{1}{x}y = \frac{1}{x^2}y^2$
    这是一个 $n=2$ 的伯努利方程。

2.  **进行变量代换:**
    令 $z = y^{1-n} = y^{1-2} = y^{-1}$。
    则 $z' = -y^{-2}y'$。

3.  **转换原方程:**
    原方程两边同除以 $y^2$：
    $y^{-2}y' + \frac{1}{x}y^{-1} = \frac{1}{x^2}$
    代入 $z$ 和 $z'$：
    $-z' + \frac{1}{x}z = \frac{1}{x^2} \implies z' - \frac{1}{x}z = -\frac{1}{x^2}$

4.  **求解关于 z 的一阶线性微分方程:**
    积分因子 $\mu(x) = e^{\int(-1/x)dx} = e^{-\ln x} = x^{-1} = \frac{1}{x}$。
    $z = \frac{1}{\mu(x)} \left[ \int \mu(x)q(x)dx + C \right] = x \left[ \int \frac{1}{x} \left(-\frac{1}{x^2}\right)dx + C \right]$
    $z = x \left[ \int -x^{-3}dx + C \right] = x \left[ \frac{1}{2}x^{-2} + C \right] = \frac{1}{2x} + Cx$

5.  **换回 y 并使用初始条件:**
    $z = 1/y \implies \frac{1}{y} = \frac{1}{2x} + Cx$
    代入 $y(1)=1$：$1 = \frac{1}{2} + C \implies C = \frac{1}{2}$。

6.  **写出特解:**
    $\frac{1}{y} = \frac{1}{2x} + \frac{x}{2} = \frac{1+x^2}{2x}$
    $y = \frac{2x}{1+x^2}$

**最终答案:** $y = \frac{2x}{1 + x^2}$

---

### **六、(8 分) 设函数 $f(x)$ 在 $[0,1]$ 上具有二阶导数且 $f''(x) \le 0$，$F(x) = \int_{0}^{1} f(y)|x - y|dy$，其中 $0 \le x \le 1$，证明 (1) $f(x) = \frac{1}{2}F''(x)$；(2) $\int_{0}^{1} f(x^2)dx \le f(1/3)$。**

**详细证明:**

**(1) 证明 $f(x) = \frac{1}{2}F''(x)$**

1.  **处理绝对值:**
    $F(x) = \int_{0}^{x} f(y)(x - y)dy + \int_{x}^{1} f(y)(y - x)dy$
    $F(x) = x\int_{0}^{x} f(y)dy - \int_{0}^{x} yf(y)dy + \int_{x}^{1} yf(y)dy - x\int_{x}^{1} f(y)dy$

2.  **求一阶导数 F'(x):**
    $F'(x) = \left(\int_{0}^{x} f(y)dy + xf(x)\right) - (xf(x)) + (-xf(x)) - \left(\int_{x}^{1} f(y)dy + x(-f(x))\right)$
    $F'(x) = \int_{0}^{x} f(y)dy - \int_{x}^{1} f(y)dy$

3.  **求二阶导数 F''(x):**
    $F''(x) = \frac{d}{dx} \left(\int_{0}^{x} f(y)dy\right) - \frac{d}{dx} \left(\int_{x}^{1} f(y)dy\right)$
    $F''(x) = f(x) - (-f(x)) = 2f(x)$
    因此，$f(x) = \frac{1}{2}F''(x)$。**证毕。**

---

**(2) 证明 $\int_{0}^{1} f(x^2)dx \le f(1/3)$**

1.  **利用凹函数性质:**
    $f''(x) \le 0$ 意味着 $f(x)$ 在 $[0, 1]$ 上是凹函数。

2.  **应用琴生不等式的积分形式 (Jensen's Inequality):**
    对于一个凹函数 $f$，有 $\frac{1}{b-a}\int_{a}^{b} f(g(x))dx \le f\left(\frac{1}{b-a}\int_{a}^{b} g(x)dx \right)$。
    在本题中，$a=0, b=1$，凹函数是 $f$，内层函数是 $g(x) = x^2$。

3.  **代入不等式:**
    $\int_{0}^{1} f(x^2)dx \le f\left(\int_{0}^{1} x^2 dx\right)$

4.  **计算右侧积分:**
    $\int_{0}^{1} x^2 dx = \left[\frac{x^3}{3}\right]_{0}^{1} = \frac{1}{3}$

5.  **结论:**
    代入计算结果，得 $\int_{0}^{1} f(x^2)dx \le f\left(\frac{1}{3}\right)$。**证毕。**