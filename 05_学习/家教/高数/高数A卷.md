
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

**1. 一阶线性方程

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


**2. 齐次方程 

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





**3. 伯努利方程 

*   **标准形式**: $y' + P(x)y = Q(x)y^n$  ($n \neq 0, 1$)
    *   **识别技巧**：它看起来**几乎**就是一个一阶线性方程，只是在右边的 $Q(x)$ 后面多了一个 $y$ 的幂次项 $y^n$。
    *   当 $n=0$ 或 $n=1$ 时，它本身就是个线性方程，无需特殊方法。
*   **核心方法**: **通过变量代换 $z = y^{1-n}$，将其巧妙地转化为一个一阶线性方程**。

#### **求解步骤:**

1.  **化为标准形式**，并两边同除以 $y^n$ (假设 $y \neq 0$)：
    $$
    y^{-n}y' + P(x)y^{1-n} = Q(x)
    $$
2.  **进行变量代换**：
    令 $z = y^{1-n}$。
3.  **求 $z$ 的导数 $z'$**：
    使用链式法则对 $z$ 求导：
    $\frac{dz}{dx} = (1-n)y^{-n} \frac{dy}{dx}$
    整理可得：$y^{-n}y' = \frac{1}{1-n}z'$。
4.  **代入方程**：
    将第 2 步和第 3 步的结果代入第 1 步的方程中：
    $$
    \frac{1}{1-n}z' + P(x)z = Q(x)
    $$
5.  **整理成标准线性形式**：
    方程两边同乘以 $(1-n)$：
    $$
    z' + (1-n)P(x)z = (1-n)Q(x)
    $$
    现在，这已经是一个关于变量 $z$ 和 $x$ 的**标准一阶线性方程**了！
6.  **使用积分因子法**求解出 $z$ 的通解。
7.  **换回 $y$**：最后一步，将 $z = y^{1-n}$ 代回去，解出 $y$，得到原方程的最终解。

**示例**: 求解 $y' + \frac{1}{x}y = x y^2$

1.  **识别**: 这是伯努利方程，其中 $P(x) = 1/x$, $Q(x) = x$, $n=2$。
2.  两边同除以 $y^2$： $y^{-2}y' + \frac{1}{x}y^{-1} = x$。
3.  **代换**: 令 $z = y^{1-2} = y^{-1}$。则 $z' = -y^{-2}y'$，即 $y^{-2}y' = -z'$。
4.  **代入**: $-z' + \frac{1}{x}z = x$。
5.  **整理**: $z' - \frac{1}{x}z = -x$。这是一个关于 $z$ 的一阶线性方程。
6.  **求解 $z$**:
    *   积分因子 $\mu(x) = e^{\int(-1/x)dx} = e^{-\ln|x|} = 1/|x|$。取 $x>0$，则 $\mu(x) = 1/x$。
    *   $z = \frac{1}{1/x} \left[ \int \frac{1}{x}(-x) dx + C \right] = x \left[ \int -1 dx + C \right] = x(-x+C) = -x^2 + Cx$。
7.  **换回 y**: $z=1/y \implies \frac{1}{y} = -x^2+Cx \implies y = \frac{1}{C x - x^2}$。

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

#### **驻点的三种类型**

驻点虽然都满足 $f'(x) = 0$，但它们周围的函数行为可以完全不同。主要有以下三种类型：

##### 1. **局部极大值点 

*   **特征**: 在这一点，函数值比它紧邻的左右两边的任何点的值都**大**。
*   **导数行为**: 在该点左侧，$f'(x) > 0$ (函数上升)；在该点右侧，$f'(x) < 0$ (函数下降)。导数符号由正变负。
*   **二阶导数判别法**: 如果 $f''(x_0) < 0$，则该驻点是局部极大值点（开口向下）。
*   **图像**: 看起来像一座**山峰**。

##### 2. **局部极小值点

*   **特征**: 在这一点，函数值比它紧邻的左右两边的任何点的值都**小**。
*   **导数行为**: 在该点左侧，$f'(x) < 0$ (函数下降)；在该点右侧，$f'(x) > 0$ (函数上升)。导数符号由负变正。
*   **二阶导数判别法**: 如果 $f''(x_0) > 0$，则该驻点是局部极小值点（开口向上）。
*   **图像**: 看起来像一个**山谷**。

##### 3. **鞍点

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
    *   $\lim_{x \to +\infty} -\frac{\sin^2 x}{x} = 0$。
    *   $\lim_{x \to 0} -\frac{\sin^2 x}{x} = \lim_{x \to 0} -\left(\frac{\sin x}{x}\right)^2 x = -1^2 \cdot 0 = 0$。
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

#### **补充：二阶齐次线性微分方程**的通用方法

这类方程的形式为：
$$
y'' + P(x)y' + Q(x)y = 0
$$
其中 $P(x)$ 和 $Q(x)$ 是关于 $x$ 的函数。

根据系数 $P(x)$ 和 $Q(x)$ 的不同，求解方法主要分为两大类：

1.  **常系数**：$P(x)$ 和 $Q(x)$ 都是常数。这是最简单、最常见的情况。
2.  **变系数**：$P(x)$ 或 $Q(x)$ 不是常数。这类方程通常更难解，只有少数特定类型有通用解法。

---

**第一类：常系数二阶齐次线性微分方程

这是最重要的类型，形式为：
$$
ay'' + by' + cy = 0 \quad (a, b, c \text{ 是常数})
$$

**核心解法：特征方程法**

**步骤 1：写出特征方程**
将微分方程中的 $y''$ 换成 $r^2$，$y'$ 换成 $r$，$y$ 换成 1：
$$
ar^2 + br + c = 0
$$

**步骤 2：求解特征根 $r$**
这是一个标准的一元二次方程，根据判别式 $\Delta = b^2 - 4ac$ 的值，分三种情况讨论根。

**步骤 3：根据根的形式写出通解**

| **判别式** | **特征根 $r_1, r_2$** | **通解 $y(x)$** |
| :--- | :--- | :--- |
| $\Delta > 0$ | 两个不等的实根 $r_1, r_2$ | $y = C_1 e^{r_1 x} + C_2 e^{r_2 x}$ |
| $\Delta = 0$ | 两个相等的实根 $r_1 = r_2 = r$ | $y = (C_1 + C_2 x) e^{rx}$ |
| $\Delta < 0$ | 一对共轭复根 $\alpha \pm i\beta$ | $y = e^{\alpha x}(C_1 \cos(\beta x) + C_2 \sin(\beta x))$ |

**示例:** 求解 $y'' - 6y' + 9y = 0$
1.  **特征方程**: $r^2 - 6r + 9 = 0$
2.  **求解根**: $(r-3)^2 = 0 \implies r_1 = r_2 = 3$ (重根)
3.  **通解**: $y = (C_1 + C_2 x)e^{3x}$

---










### **四、(8 分) 设函数 $\phi(x)$ 连续，$f(x) = \int_{0}^{1} \phi(xt)dt$，且 $\lim_{x \to 0} \frac{\phi(x)}{x} = 2019$，求 $f(x)$ 和 $f'(x)$。**

**详细步骤:**

1.  **对 f(x) 进行变量代换:**
    令 $u = xt$，则 $dt = \frac{1}{x}du$。
    $f(x) = \int_{0}^{x} \phi(u) \frac{1}{x}du = \frac{1}{x} \int_{0}^{x} \phi(u)du$

2.  **利用极限条件确定 $\phi(x)$ 的形式:**
    $\lim_{x \to 0} \frac{\phi(x)}{x} = 2019$。因为极限存在且分母趋于0，所以分子 $\lim_{x \to 0} \phi(x) = 0$。又因为 $\phi(x)$ 连续，所以 $\phi(0) = 0$。
    该极限是导数的定义：$\phi'(0) = \lim_{x \to 0} \frac{\phi(x) - \phi(0)}{x - 0} = 2019$。
    满足此条件的连续函数形式是 $\phi(x) = 2019x$。

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

**方法一：伯努利方程**
标准形式为 $y' + P(x)y = Q(x)y^n$。

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


**方法二：齐次方程**

**第一步：变形，凑出 $\frac{y}{x}$ 的形式**

$$
\frac{dy}{dx} = \frac{y^2}{x^2} - \frac{xy}{x^2}
$$
$$
\frac{dy}{dx} = \left(\frac{y}{x}\right)^2 - \frac{y}{x}
$$
**第二步：变量代换**

根据齐次方程的解法，我们进行变量代换：
*   令 $v = \frac{y}{x}$，则 $y = vx$。
*   对 $y = vx$ 求导，得到 $\frac{dy}{dx} = v + x \frac{dv}{dx}$。

**第三步：代入方程，分离变量**

将上述结果代入我们转化好的齐次方程 $\frac{dy}{dx} = \left(\frac{y}{x}\right)^2 - \frac{y}{x}$ 中：
$$
v + x \frac{dv}{dx} = v^2 - v
$$
现在，整理这个关于 $v$ 和 $x$ 的新方程，并分离变量：
$$
x \frac{dv}{dx} = v^2 - 2v
$$
$$
\frac{dv}{v^2 - 2v} = \frac{dx}{x}
$$

**第四步：两边积分**

现在我们对两边进行积分。左边的积分需要使用**部分分式分解**。
$$
\int \frac{1}{v(v-2)} dv = \int \frac{1}{x} dx
$$
对左边的被积函数进行分解：
$$
\frac{1}{v(v-2)} = \frac{A}{v} + \frac{B}{v-2}
$$
通分得到 $1 = A(v-2) + Bv$。
*   令 $v=0$，得 $1 = -2A \implies A = -1/2$。
*   令 $v=2$，得 $1 = 2B \implies B = 1/2$。

所以，左边的积分变为：
$$
\int \left(\frac{1/2}{v-2} - \frac{1/2}{v}\right) dv = \frac{1}{2} \int \left(\frac{1}{v-2} - \frac{1}{v}\right) dv
$$
$$
\frac{1}{2} (\ln|v-2| - \ln|v|) = \ln|x| + C_1
$$
利用对数性质合并：
$$
\frac{1}{2} \ln\left|\frac{v-2}{v}\right| = \ln|x| + C_1
$$
$$
\ln\left|\frac{v-2}{v}\right| = 2\ln|x| + 2C_1 = \ln(x^2) + C_2
$$
两边取指数：
$$
\left|\frac{v-2}{v}\right| = e^{\ln(x^2) + C_2} = e^{\ln(x^2)} \cdot e^{C_2} = C_3 x^2 \quad (C_3 = e^{C_2} > 0)
$$
去掉绝对值，让常数 $C$ 可以取正负值：
$$
\frac{v-2}{v} = C x^2
$$

**第五步：换回 $y$**

将 $v = \frac{y}{x}$ 代回：
$$
\frac{\frac{y}{x} - 2}{\frac{y}{x}} = C x^2
$$
$$
1 - \frac{2x}{y} = C x^2
$$
$$
1 - C x^2 = \frac{2x}{y}
$$
$$
y = \frac{2x}{1 - C x^2}
$$
这就是方程的通解。

**第七步：使用初始条件求特解**

代入初始条件 $y(1) = 1$：
$$
1 = \frac{2(1)}{1 - C (1)^2} = \frac{2}{1-C}
$$
$$
1 - C = 2 \implies C = -1
$$

**第八步：写出最终特解**

将 $C = -1$ 代入通解：
$$
y = \frac{2x}{1 - (-1)x^2} = \frac{2x}{1+x^2}
$$


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

2.  **应用琴生不等式的积分形式:**
    对于一个凹函数 $f$，有 $\frac{1}{b-a}\int_{a}^{b} f(g(x))dx \le f\left(\frac{1}{b-a}\int_{a}^{b} g(x)dx \right)$。
    在本题中，$a=0, b=1$，凹函数是 $f$，内层函数是 $g(x) = x^2$。

3.  **代入不等式:**
    $\int_{0}^{1} f(x^2)dx \le f\left(\int_{0}^{1} x^2 dx\right)$

4.  **计算右侧积分:**
    $\int_{0}^{1} x^2 dx = \left[\frac{x^3}{3}\right]_{0}^{1} = \frac{1}{3}$

5.  **结论:**
    代入计算结果，得 $\int_{0}^{1} f(x^2)dx \le f\left(\frac{1}{3}\right)$。**证毕。**

# 2025（选填部分）


### **一、单项选择题 (每小题 4 分, 共 20 分)**

**1. 曲线 $x=2+\sin t, y=\cos t, z=(1+t)^2$ 在对应于 $t=0$ 的点处的切线方程为 [A]**

*   **思路:**
    1.  求出 $t=0$ 时点的坐标 $(x_0, y_0, z_0)$。
    2.  求出曲线在 $t=0$ 处的切向量 $(x'(0), y'(0), z'(0))$。
    3.  写出空间曲线切线的**点向式方程**。

*   **步骤:**
    1.  当 $t=0$ 时：
        $x_0 = 2 + \sin(0) = 2$
        $y_0 = \cos(0) = 1$
        $z_0 = (1+0)^2 = 1$
        所以切点坐标为 $(2, 1, 1)$。
    2.  求导数：
        $x'(t) = \cos t$
        $y'(t) = -\sin t$
        $z'(t) = 2(1+t)$
    3.  当 $t=0$ 时，切向量为：
        $x'(0) = \cos(0) = 1$
        $y'(0) = -\sin(0) = 0$
        $z'(0) = 2(1+0) = 2$
        切向量为 $\vec{v} = (1, 0, 2)$。
    4.  切线方程为：
        $$
        \frac{x-x_0}{x'(0)} = \frac{y-y_0}{y'(0)} = \frac{z-z_0}{z'(0)}
        $$
        $$
        \frac{x-2}{1} = \frac{y-1}{0} = \frac{z-1}{2}
        $$










**2. 直线 $\frac{x-2}{3} = \frac{y+1}{-1} = \frac{z-2}{2}$ 与平面 $2x-4y+z=7$ 的夹角 $\theta$ 的正弦值 $\sin\theta$ 为 [ C ]**

*   **思路:**
    直线与平面的夹角 $\theta$ 的正弦值等于直线方向向量 $\vec{s}$ 与平面法向量 $\vec{n}$ 夹角 $\phi$ 的余弦值的绝对值，即 $\sin\theta = |\cos\phi| = \frac{|\vec{s} \cdot \vec{n}|}{||\vec{s}|| \cdot ||\vec{n}||}$。

*   **步骤:**
    1.  直线方向向量：$\vec{s} = (3, -1, 2)$
    2.  平面法向量：$\vec{n} = (2, -4, 1)$
    3.  计算点积：$|\vec{s} \cdot \vec{n}| = |(3)(2) + (-1)(-4) + (2)(1)| = |6 + 4 + 2| = 12$
    4.  计算模长：
        $||\vec{s}|| = \sqrt{3^2 + (-1)^2 + 2^2} = \sqrt{9+1+4} = \sqrt{14}$
        $||\vec{n}|| = \sqrt{2^2 + (-4)^2 + 1^2} = \sqrt{4+16+1} = \sqrt{21}$
    5.  计算 $\sin\theta$：
        $$
        \sin\theta = \frac{12}{\sqrt{14} \cdot \sqrt{21}} = \frac{12}{\sqrt{2 \cdot 7} \cdot \sqrt{3 \cdot 7}} = \frac{12}{7\sqrt{6}} = \frac{12\sqrt{6}}{7 \cdot 6} = \frac{2\sqrt{6}}{7}
        $$

### **1. 如何从直线方程确定方向向量**

直线的方程通常有多种形式，这里我们主要看**对称式（点向式）方程**，因为这是题目中给出的形式。

#### **直线的对称式方程**

其标准形式为：
$$
\frac{x - x_0}{a} = \frac{y - y_0}{b} = \frac{z - z_0}{c}
$$

**如何解读这个方程：**

*   **分子**告诉我们直线经过的一个**点**。
    *   通过令分子为零，我们可以得到点的坐标：
        *   $x - x_0 = 0 \implies x = x_0$
        *   $y - y_0 = 0 \implies y = y_0$
        *   $z - z_0 = 0 \implies z = z_0$
    *   所以，这条直线经过点 $P_0(x_0, y_0, z_0)$。

*   **分母**直接给出了直线的**方向向量**。
    *   这个向量描述了直线在空间中的朝向。
    *   方向向量 $\vec{s} = (a, b, c)$。


---

### **2. 如何从平面方程确定法向量**

平面的方程通常以**一般式**给出。

#### **平面的一般式方程**

其标准形式为：
$$
Ax + By + Cz + D = 0
$$

方程中 $x, y, z$ 的**系数** $A, B, C$ 直接构成了该平面的一个**法向量**。

**如果一个向量与平面上的任意向量都垂直，那么这个向量就是该平面的法向量。**
向量 $\vec{n} = (A, B, C)$ 正好满足这个条件。

**推导过程

1.  **在平面上任取两个点**
    假设我们有一个由方程 $Ax + By + Cz + D = 0$ 定义的平面。
    让我们在这个平面上随便选取两个不同的点：
    *   点 $P_0 = (x_0, y_0, z_0)$
    *   点 $P_1 = (x_1, y_1, z_1)$

2.  **利用点在平面上的性质**
    因为这两个点都在平面上，所以它们的坐标都满足平面方程：
    *   对于 $P_0$: $Ax_0 + By_0 + Cz_0 + D = 0$  --- (1)
    *   对于 $P_1$: $Ax_1 + By_1 + Cz_1 + D = 0$  --- (2)

3.  **构造一个位于平面内的向量**
    现在，我们构造一个从 $P_0$ 指向 $P_1$ 的向量 $\vec{P_0 P_1}$。这个向量的起点和终点都在平面上，所以**整个向量 $\vec{P_0 P_1}$ 完全位于这个平面内**。

    该向量的坐标表示为：
    $$
    \vec{P_0 P_1} = (x_1 - x_0, y_1 - y_0, z_1 - z_0)
    $$

4.  **证明系数向量 $\vec{n}$ 与平面内向量 $\vec{P_0 P_1}$ 垂直**
    要证明两个向量垂直，我们只需要证明它们的**点积（内积）等于零**。

    让我们来计算向量 $\vec{n} = (A, B, C)$ 和向量 $\vec{P_0 P_1}$ 的点积：
    $$
    \vec{n} \cdot \vec{P_0 P_1} = (A, B, C) \cdot (x_1 - x_0, y_1 - y_0, z_1 - z_0)
    $$
    $$
    = A(x_1 - x_0) + B(y_1 - y_0) + C(z_1 - z_0)
    $$
    展开这个表达式：
    $$
    = (Ax_1 + By_1 + Cz_1) - (Ax_0 + By_0 + Cz_0)
    $$
    现在，我们可以用第 2 步得到的方程 (1) 和 (2) 来替换这两个括号里的内容：
    *   从 (2) 式可知： $Ax_1 + By_1 + Cz_1 = -D$
    *   从 (1) 式可知： $Ax_0 + By_0 + Cz_0 = -D$

    代入点积的计算结果：
    $$
    \vec{n} \cdot \vec{P_0 P_1} = (-D) - (-D) = -D + D = 0
    $$









**3. 极限 $\lim_{(x,y)\to(0,0)} \frac{\sqrt{1+x^2y^2}-1}{(x^2+y^2)e^{\sin(xy^2)}}$ = [ A ]**

*   **思路:**
    这是一个 "0/0" 型极限。我们可以使用等价无穷小来简化计算。

*   **步骤:**
    1.  **分子等价无穷小:** 
        在这里，令 $u = x^2y^2$。当 $(x,y) \to (0,0)$ 时，$u \to 0$。
        所以，分子 $\sqrt{1+x^2y^2}-1 \sim \frac{1}{2}x^2y^2$。
    2.  **分母中的指数部分:** 
		当 $(x,y) \to (0,0)$ 时，$xy^2 \to 0$，所以 $\sin(xy^2) \to 0$。
        因此，分母中的 $e^{\sin(xy^2)} \to e^0 = 1$。
    3.  **代入等价无穷小:**
        $$
        \text{原极限} = \lim_{(x,y)\to(0,0)} \frac{\frac{1}{2}x^2y^2}{(x^2+y^2) \cdot 1} = \frac{1}{2} \lim_{(x,y)\to(0,0)} \frac{x^2y^2}{x^2+y^2}
        $$
    4.  **计算剩余极限 (使用极坐标、夹逼定理):**
        令 $x = r\cos\theta, y=r\sin\theta$。当 $(x,y) \to (0,0)$ 时，$r \to 0$。
        $$
        \frac{x^2y^2}{x^2+y^2} = \frac{(r^2\cos^2\theta)(r^2\sin^2\theta)}{r^2} = r^2\cos^2\theta\sin^2\theta
        $$
        因为 $0 \le \cos^2\theta\sin^2\theta \le 1$，所以：
        $$
        0 \le r^2\cos^2\theta\sin^2\theta \le r^2
        $$
        当 $r \to 0$ 时，$r^2 \to 0$。根据夹逼定理，$\lim_{r\to 0} r^2\cos^2\theta\sin^2\theta = 0$。
    5.  **最终结果:**
        $$
        L = \frac{1}{2} \times 0 = 0
        $$











**4. 设积分区域 $D: x^2+y^2 \le 1, y \ge 0$，则二重积分 $\iint_D \sqrt{1-x^2-y^2} \,dxdy$ = [ B ]**

*   **思路:**
    该积分的几何意义是求一个球心在原点、半径为1的半球体（$z=\sqrt{1-x^2-y^2}$）在 xy 平面上半圆区域 $D$ 上方的体积。这相当于整个球体积的 $1/4$。

*   **步骤 (几何法):**
    1.  球的体积公式为 $V_{sphere} = \frac{4}{3}\pi r^3$。
    2.  半径 $r=1$。
    3.  所求体积为 $V = \frac{1}{4} V_{sphere} = \frac{1}{4} \cdot \frac{4}{3}\pi (1)^3 = \frac{\pi}{3}$。

*   **步骤 (极坐标法):**
    1.  **坐标转换:**
        $x = r\cos\theta, y = r\sin\theta$, $dxdy = r\,drd\theta$
        $\sqrt{1-x^2-y^2} = \sqrt{1-r^2}$
    2.  **积分区域:**
        $0 \le r \le 1$, $0 \le \theta \le \pi$
    3.  **计算积分:**
        $$
        I = \int_{0}^{\pi} d\theta \int_{0}^{1} \sqrt{1-r^2} \cdot r \, dr = \pi \left[ -\frac{1}{2} \cdot \frac{2}{3}(1-r^2)^{3/2} \right]_0^1 = \pi \left( 0 - (-\frac{1}{3}) \right) = \frac{\pi}{3}
        $$


以下是求解二重积分最主要的几种方法：

---

### 1. 利用直角坐标系计算（化为二次积分）

这是最基本、最直接的方法。核心思想是将二重积分 $\iint_D f(x,y) \,dA$ 转化为两次独立的定积分（即二次积分）。

#### 类型一：X-型区域 (上下边界是函数)

如果积分区域 $D$ 可以被描述为：
$$
D = \{(x,y) \mid a \le x \le b, \quad g_1(x) \le y \le g_2(x)\}
$$
*   **几何特征**：区域的左右边界是两条竖直线 $x=a, x=b$，上下边界是两条函数曲线 $y=g_1(x)$ 和 $y=g_2(x)$。
*   **积分次序**：**先对 $y$ 积分，再对 $x$ 积分**。
*   **公式**：
    $$
    \iint_D f(x,y) \,dA = \int_a^b \left[ \int_{g_1(x)}^{g_2(x)} f(x,y) \,dy \right] dx
    $$
    *   **计算技巧**：在计算内层关于 $y$ 的积分时，要将 $x$ 视为**常数**。

#### 类型二：Y-型区域 (左右边界是函数)

如果积分区域 $D$ 可以被描述为：
$$
D = \{(x,y) \mid c \le y \le d, \quad h_1(y) \le x \le h_2(y)\}
$$
*   **几何特征**：区域的上下边界是两条水平线 $y=c, y=d$，左右边界是两条函数曲线 $x=h_1(y)$ 和 $x=h_2(y)$。
*   **积分次序**：**先对 $x$ 积分，再对 $y$ 积分**。
*   **公式**：
    $$
    \iint_D f(x,y) \,dA = \int_c^d \left[ \int_{h_1(y)}^{h_2(y)} f(x,y) \,dx \right] dy
    $$
    *   **计算技巧**：在计算内层关于 $x$ 的积分时，要将 $y$ 视为**常数**。

**选择哪种类型？**
*   **画图是关键！** 首先画出积分区域。
*   观察区域形状，判断是描述成 X-型区域更简单，还是 Y-型区域更简单。
*   有时一个区域需要被**分割**成几个不同的小区域才能用这种方法计算。
*   有时一种积分次序非常难算，而**交换积分次序**后会变得很简单。

---

### 2. 利用极坐标系计算

*   **适用场景**：
    *   积分区域是**圆形、扇形、环形**或其一部分。
    *   被积函数中含有 $x^2+y^2$、$\frac{y}{x}$ 等形式。

*   **核心转换公式**：
    1.  **坐标代换**:
        $x = r\cos\theta$
        $y = r\sin\theta$
    2.  **面积微元代换
        $dxdy \rightarrow r \,dr d\theta$  
    3.  **被积函数代换**:
        $f(x,y) \rightarrow f(r\cos\theta, r\sin\theta)$
    4.  **积分区域转换**:
        将直角坐标系下的区域 $D$ 转化为极坐标系下的区域 $D^*$。

*   **公式**：
    $$
    \iint_D f(x,y) \,dxdy = \iint_{D^*} f(r\cos\theta, r\sin\theta) \,r \,dr d\theta
    $$

---

### 3. 利用对称性简化计算

如果积分区域或被积函数具有对称性，可以大大减少计算量。

#### 区域对称性 + 函数奇偶性

设积分区域 $D$ 关于 **y 轴** 对称。
*   如果 $f(x,y)$ 是关于 $x$ 的**奇函数** (即 $f(-x, y) = -f(x, y)$)，那么：
    $$
    \iint_D f(x,y) \,dA = 0
    $$
*   如果 $f(x,y)$ 是关于 $x$ 的**偶函数** (即 $f(-x, y) = f(x, y)$)，那么：
    $$
    \iint_D f(x,y) \,dA = 2 \iint_{D_1} f(x,y) \,dA
    $$
    其中 $D_1$ 是 $D$ 在 y 轴右侧的一半区域。

**关于 x 轴对称的情况同理。**

**示例：** 计算 $\iint_D x \,dA$，其中 $D$ 是单位圆盘 $x^2+y^2 \le 1$。
*   区域 $D$ 关于 y 轴对称。
*   被积函数 $f(x,y) = x$ 是关于 $x$ 的奇函数。
*   所以，$\iint_D x \,dA = 0$。

---

### 4. 利用二重积分的几何意义

有时，一个复杂的二重积分可能具有非常简单的几何或物理意义，可以直接写出答案而无需计算。

*   **几何意义**：
    *   如果 $f(x,y) = 1$，则 $\iint_D 1 \,dA = \text{Area}(D)$，即积分值就是区域 $D$ 的**面积**。
    *   如果 $f(x,y) \ge 0$，则 $\iint_D f(x,y) \,dA$ 代表了以 $D$ 为底，曲面 $z=f(x,y)$ 为顶的**曲顶柱体的体积**。


**示例：** 计算 $\iint_D \sqrt{1-x^2-y^2} \,dA$，其中 $D$ 是单位圆盘 $x^2+y^2 \le 1$。
*   这个积分计算的是以单位圆盘为底，以半球面 $z=\sqrt{1-x^2-y^2}$ 为顶的物体的体积。
*   这正是一个**半球体**的体积。
*   半球体积公式为 $\frac{2}{3}\pi r^3$。半径 $r=1$，所以积分值是 $\frac{2}{3}\pi$。

### 总结

| **方法**      | **适用场景**                   | **关键点**                         |
| :---------- | :------------------------- | :------------------------------ |
| **直角坐标**    | 矩形、梯形、三角形等由直线或简单函数围成的区域    | **画图**，确定积分上下限和积分次序             |
| **极坐标**     | 圆形、扇形、环形区域；被积函数含 $x^2+y^2$ | 坐标代换，面积微元 $dxdy \to rdrd\theta$ |
| **对称性**     | 区域对称，函数有奇偶性                | 判断奇偶性，快速得出 0 或简化计算              |
| **几何/物理意义** | 被积函数是 1，或代表简单几何体/物理量       | 识别几何形状或物理模型，直接用公式               |







**5. 方程 $16x^2 + 4y^2 - z^2 = 64$ 对应的曲面示意图为 [ D ]**

*   **思路:**
    将方程化为标准形式，然后判断曲面类型。

*   **步骤:**
    1.  **化为标准形式:**
        方程两边同除以 64：
        $$
        \frac{16x^2}{64} + \frac{4y^2}{64} - \frac{z^2}{64} = 1 \implies \frac{x^2}{4} + \frac{y^2}{16} - \frac{z^2}{64} = 1
        $$
    2.  **判断类型:**
        这是**单叶双曲面的标准方程 ($\frac{x^2}{a^2} + \frac{y^2}{b^2} - \frac{z^2}{c^2} = 1$)。它的对称轴是与负号项对应的 z 轴。
    3.  **对比图像:**
        (A) 双曲抛物面（马鞍面）。
        (B) 抛物柱面。
        (C) 椭圆抛物面。
        (D) 单叶双曲面。


| 曲面名称 | 标准方程形式 ($z$ 为特殊轴) | 关键特征 |
| :--- | :--- | :--- |
| **椭球面** | $\frac{x^2}{a^2}+\frac{y^2}{b^2}+\frac{z^2}{c^2}=1$ | 3个正二次项, =1 (封闭曲面) |
| **椭圆抛物面 (碗)** | $\frac{x^2}{a^2}+\frac{y^2}{b^2}=z$ | 2个同号二次项, 1个一次项 |
| **双曲抛物面 (马鞍)** | $\frac{x^2}{a^2}-\frac{y^2}{b^2}=z$ | 2个异号二次项, 1个一次项 |
| **单叶双曲面** | $\frac{x^2}{a^2}+\frac{y^2}{b^2}-\frac{z^2}{c^2}=1$ | 2正1负二次项, =1 (一个整体) |
| **双叶双曲面** | $-\frac{x^2}{a^2}-\frac{y^2}{b^2}+\frac{z^2}{c^2}=1$ | 1正2负二次项, =1 (两片分离) |
| **椭圆锥面** | $\frac{x^2}{a^2}+\frac{y^2}{b^2}-\frac{z^2}{c^2}=0$ | 3个二次项, =0 |
| **柱面** | 方程只含两个变量 | 沿缺失变量轴平移基准线 |

---



### **二、填空题 (每小题 4 分, 共 24 分)**

**1. 设单位向量 $\vec{a}, \vec{b}, \vec{c}$ 满足条件 $\vec{a}+\vec{b}+\vec{c}=\vec{0}$，则 $\vec{a}\cdot\vec{b}+\vec{b}\cdot\vec{c}+\vec{c}\cdot\vec{a} = \underline{-\frac{3}{2}}$**

*   **思路:** 对向量和式 $\vec{a}+\vec{b}+\vec{c}=\vec{0}$ 进行平方（与自身点积）。
*   **步骤:**
    $(\vec{a}+\vec{b}+\vec{c}) \cdot (\vec{a}+\vec{b}+\vec{c}) = \vec{0} \cdot \vec{0} = 0$
    $|\vec{a}|^2 + |\vec{b}|^2 + |\vec{c}|^2 + 2(\vec{a}\cdot\vec{b}+\vec{b}\cdot\vec{c}+\vec{c}\cdot\vec{a}) = 0$
    $1^2 + 1^2 + 1^2 + 2(\dots) = 0$
    $3 + 2(\dots) = 0 \implies \vec{a}\cdot\vec{b}+\vec{b}\cdot\vec{c}+\vec{c}\cdot\vec{a} = -\frac{3}{2}$



**2. 设 $\vec{OA} = \vec{i}+\vec{j}+3\vec{k}$, $\vec{OB} = \vec{i}+2\vec{k}$，则 $\triangle OAB$ 的面积 $S_{\triangle OAB} = \underline{\frac{\sqrt{6}}{2}}$**

*   **思路:** 由两个向量构成的三角形的面积等于这两个向量叉积所得向量的模长的一半。
*   **步骤:**
    1.  $\vec{OA} = (1, 1, 3)$, $\vec{OB} = (1, 0, 2)$
    2.  计算叉积 $\vec{OA} \times \vec{OB}$：
        $$
        \begin{vmatrix} \vec{i} & \vec{j} & \vec{k} \\ 1 & 1 & 3 \\ 1 & 0 & 2 \end{vmatrix} = 2\vec{i} + \vec{j} - \vec{k} = (2, 1, -1)
        $$
    3.  $S_{\triangle OAB} = \frac{1}{2} ||\vec{OA} \times \vec{OB}|| = \frac{1}{2} \sqrt{2^2 + 1^2 + (-1)^2} = \frac{\sqrt{6}}{2}$


为什么三角形的面积会和两个向量的“叉积”有关呢？

这背后其实联系了三个关键的几何概念：
1.  **三角形的面积公式** 
2.  **向量点积与夹角的关系**
3.  **向量叉积的定义与几何意义**

### 1. 从最基础的三角形面积公式说起

$$
S_{\triangle OAB} = \frac{1}{2} ||\vec{OA}|| \cdot ||\vec{OB}|| \sin\theta
$$
### 2. 向量叉积的几何意义

**向量叉积**这个运算，它的**大小（模长）在几何上的定义就是：
$$
||\vec{OA} \times \vec{OB}|| = ||\vec{OA}|| \cdot ||\vec{OB}|| \sin\theta
$$
### 3. 如何计算叉积？—— 行列式方法

前面我们理解了“为什么”可以用叉积来算面积。现在的问题就是“怎么算”叉积。

对于两个三维向量 $\vec{a} = (a_1, a_2, a_3)$ 和 $\vec{b} = (b_1, b_2, b_3)$，它们的叉积 $\vec{a} \times \vec{b}$ 的计算规则可以用一个 3x3 的行列式来方便地记忆和计算：
$$
\vec{a} \times \vec{b} = \begin{vmatrix} \vec{i} & \vec{j} & \vec{k} \\ a_1 & a_2 & a_3 \\ b_1 & b_2 & b_3 \end{vmatrix}
$$
其中 $\vec{i}, \vec{j}, \vec{k}$ 是 x, y, z 方向的单位向量。

*   展开这个行列式，就得到了叉积的结果向量。







**3. 已知曲面 $z=7-x^2+y^2$ 上点 $M$ 处的切平面平行于平面 $2x-4y-z+5=0$，则点 $M$ 的坐标为 $\underline{(-1, -2, 10)}$**

*   **思路:**
    两平面平行，则它们的法向量平行。
*   **步骤:**
    1.  ==曲面 $f(x,y) = 7-x^2+y^2$ 在点 $M(x_0, y_0, z_0)$ 处的法向量为 $\vec{n}_M = (f_x, f_y, -1) = (-2x_0, 2y_0, -1)$。==
    2.  给定平面的法向量为 $\vec{n}_P = (2, -4, -1)$。
    3.  因为 $\vec{n}_M$ 与 $\vec{n}_P$ 平行，所以 $\vec{n}_M = k\vec{n}_P$。
        $(-2x_0, 2y_0, -1) = k(2, -4, -1)$
    4.  比较 $z$ 分量: $-1 = -k \implies k=1$。
    5.  比较 $x, y$ 分量:
        $-2x_0 = 2k = 2 \implies x_0 = -1$
        $2y_0 = -4k = -4 \implies y_0 = -2$
    6.  点 $M$ 在曲面上，代入求 $z_0$：
        $z_0 = 7 - (-1)^2 + (-2)^2 = 7 - 1 + 4 = 10$


### 方法一

**核心思想：** 任何曲面都可以看作一个三维空间中标量场 $F(x,y,z)$ 的**等值面**（Level Surface）。而这个标量场的**梯度向量 $\nabla F$ 在任意一点都与该点所在的等值面垂直**。因此，梯度向量就是法向量。

**步骤分解：**

1.  **将显式函数 $z = f(x,y)$ 转化为隐式函数 $F(x,y,z)=0$ 的形式。**
    我们只需要把所有项移到等式的一边即可。
    对于题目中的曲面 $z = 7 - x^2 + y^2$，我们可以写成：
    $$
    7 - x^2 + y^2 - z = 0
    $$
    我们定义一个三元函数 $F(x,y,z) = 7 - x^2 + y^2 - z$。那么，我们的曲面就是 $F(x,y,z)=0$ 这个**等值面**。

2.  **计算梯度向量 $\nabla F$。**
    梯度向量的定义是：
    $$
    \nabla F = \left( \frac{\partial F}{\partial x}, \frac{\partial F}{\partial y}, \frac{\partial F}{\partial z} \right)
    $$
    现在我们来计算 $F(x,y,z) = 7 - x^2 + y^2 - z$ 的偏导数：
    *   $\frac{\partial F}{\partial x} = -2x$
    *   $\frac{\partial F}{\partial y} = 2y$
    *   $\frac{\partial F}{\partial z} = -1$

    所以，梯度向量是：
    $$
    \nabla F = (-2x, 2y, -1)
    $$

3.  **得出结论。**
    梯度向量 $\nabla F$ 就是曲面在点 $(x,y,z)$ 处的法向量 $\vec{n}$。
    因此，法向量 $\vec{n} = (-2x, 2y, -1)$。

    在题目中的特定点 $M(x_0, y_0, z_0)$，法向量就是：
    $$
    \vec{n}_M = (-2x_0, 2y_0, -1)
    $$

---

### 方法二

这种方法能从几何上直观地理解为什么是这个形式。

**核心思想：** 曲面上某一点的**法向量**，必然与该点上**所有切线**都垂直。我们只需要在该点找到两条**不平行**的切线，计算它们的**叉积**，得到的结果就必然是法向量（因为叉积的结果向量同时垂直于原来的两个向量）。

**步骤分解：**

1.  **在曲面 $z=f(x,y)$ 上构造两条特殊的曲线。**
    我们考虑过曲面上一点 $P_0(x_0, y_0, z_0)$ 的两条最简单的曲线：
    *   **曲线1 ($C_1$)**: 保持 $y$ 不变 ($y=y_0$)，让 $x$ 变化。这条曲线可以参数化为：
        $\vec{r}_1(t) = (t, y_0, f(t, y_0))$。
        我们关心的是 $t=x_0$ 这一点。
    *   **曲线2 ($C_2$)**: 保持 $x$ 不变 ($x=x_0$)，让 $y$ 变化。这条曲线可以参数化为：
        $\vec{r}_2(t) = (x_0, t, f(x_0, t))$。
        我们关心的是 $t=y_0$ 这一点。

2.  **求这两条曲线在点 $P_0$ 处的切向量。**
    *   **$C_1$ 的切向量 $\vec{v}_1$**: 对 $\vec{r}_1(t)$ 求导：
        $\vec{r}_1'(t) = (1, 0, \frac{\partial f}{\partial x}(t, y_0))$。
        在 $t=x_0$ 处，切向量为 $\vec{v}_1 = (1, 0, f_x(x_0, y_0))$。
    *   **$C_2$ 的切向量 $\vec{v}_2$**: 对 $\vec{r}_2(t)$ 求导：
        $\vec{r}_2'(t) = (0, 1, \frac{\partial f}{\partial y}(x_0, t))$。
        在 $t=y_0$ 处，切向量为 $\vec{v}_2 = (0, 1, f_y(x_0, y_0))$。

    这两个向量 $\vec{v}_1$ 和 $\vec{v}_2$ 构成了点 $P_0$ 处的**切平面**。

3.  **计算两个切向量的叉积，得到法向量 $\vec{n}$。**
    $$
    \vec{n} = \vec{v}_1 \times \vec{v}_2 = \begin{vmatrix} \vec{i} & \vec{j} & \vec{k} \\ 1 & 0 & f_x \\ 0 & 1 & f_y \end{vmatrix}
    $$
    展开行列式：
    *   $\vec{i}$ 分量: $(0 \cdot f_y - f_x \cdot 1) = -f_x$
    *   $\vec{j}$ 分量: $-(1 \cdot f_y - f_x \cdot 0) = -f_y$
    *   $\vec{k}$ 分量: $(1 \cdot 1 - 0 \cdot 0) = 1$

    所以，我们得到一个法向量 $\vec{n} = (-f_x, -f_y, 1)$。

**这个结果和我们的公式 $\vec{n} = (f_x, f_y, -1)$ 有什么关系？**
向量 $(-f_x, -f_y, 1)$ 和 $(f_x, f_y, -1)$ 只是方向完全相反，大小相等。它们都是平面的法向量（法向量的方向可以朝上也可以朝下）。所以这两个公式是等价的。

### 总结

| 方法       | 适用函数              | 核心思想      | 最终公式                                                                                                                 |
| :------- | :---------------- | :-------- | :------------------------------------------------------------------------------------------------------------------- |
| **梯度法**  | 隐式函数 $F(x,y,z)=C$ | 梯度垂直于等值面  | $\vec{n} = \nabla F = (\frac{\partial F}{\partial x}, \frac{\partial F}{\partial y}, \frac{\partial F}{\partial z})$ |
| **快捷公式** | 显式函数 $z=f(x,y)$   | 是梯度法的一种特例 | $\vec{n} = (f_x, f_y, -1)$ 或 $\vec{n} = (-f_x, -f_y, 1)$                                                             |

对于 $z=f(x,y)$ 形式的曲面，直接使用快捷公式 $\vec{n} = (f_x, f_y, -1)$ 是最方便的，而它的理论基础就是**梯度法**。





**4. 设二元函数 $f(x,y)$ 具有一阶连续偏导数，且 $f(0,0)=0$ 以及梯度 $\text{grad}f(x,y) = \frac{2y}{1+x^2y^2}\vec{i} + \frac{2x}{1+x^2y^2}\vec{j}$，那么 $f(x,y) = \underline{2\arctan(xy)}$**

*   **思路:**
    根据梯度定义，我们有 $\frac{\partial f}{\partial x}$ 和 $\frac{\partial f}{\partial y}$。通过积分可以求出 $f(x,y)$。
*   **步骤:**
    1.  $\frac{\partial f}{\partial x} = \frac{2y}{1+(xy)^2}$
    2.  对 $x$ 积分：
        $f(x,y) = \int \frac{2y}{1+(xy)^2} dx$
        令 $u=xy$, $du = y\,dx$
        $f(x,y) = \int \frac{2}{1+u^2} du = 2\arctan(u) + h(y) = 2\arctan(xy) + h(y)$
    3.  对结果求 $y$ 的偏导数：
        $\frac{\partial f}{\partial y} = \frac{\partial}{\partial y}(2\arctan(xy)+h(y)) = \frac{2x}{1+x^2y^2} + h'(y)$
    4.  令其等于已知的 $\frac{\partial f}{\partial y} = \frac{2x}{1+x^2y^2}$，得到 $h'(y)=0$，所以 $h(y)=C$。
    5.  $f(x,y) = 2\arctan(xy) + C$。
    6.  使用初始条件 $f(0,0)=0$：
        $f(0,0) = 2\arctan(0) + C = 0 \implies C=0$
    7.  最终 $f(x,y) = 2\arctan(xy)$。











**5. 设 $a>1$，交换二次积分 $\int_{0}^{1} dx \int_{0}^{x^2} f(x,y) dy + \int_{1}^{a} dx \int_{0}^{\frac{a-x}{a-1}} f(x,y) dy$ 的积分次序为 $\underline{\int_{0}^{1} dy \int_{\sqrt{y}}^{a-(a-1)y} f(x,y) dx}$**

*   **思路:**
    1.  画出两个积分对应的积分区域 $D_1$ 和 $D_2$。
    2.  合并区域 $D = D_1 \cup D_2$。
    3.  将合并后的区域 $D$ 描述为 X-简单区域（先积 $x$ 后积 $y$）。

*   **步骤:**
    1.  **区域 $D_1$:** 由 $\int_{0}^{1} dx \int_{0}^{x^2} f(x,y) dy$ 给出。
        $0 \le x \le 1$
        $0 \le y \le x^2$
        这是由 $x$ 轴，$x=1$ 和抛物线 $y=x^2$ 围成的区域。

    2.  **区域 $D_2$:** 由 $\int_{1}^{a} dx \int_{0}^{\frac{a-x}{a-1}} f(x,y) dy$ 给出。
        $1 \le x \le a$
        $0 \le y \le \frac{a-x}{a-1}$
        下边界是 $y=0$ ($x$ 轴)，上边界是 $y = \frac{a-x}{a-1}$。这是一条直线，当 $x=1$ 时 $y=1$；当 $x=a$ 时 $y=0$。所以这条直线连接点 $(1,1)$ 和 $(a,0)$。
        $D_2$ 是由 $x=1$, $y=0$ 和直线 $y = \frac{a-x}{a-1}$ 围成的三角形区域。

    3.  **合并区域 $D$:**
        $D_1$ 和 $D_2$ 在直线 $x=1$ 处无缝拼接，合并成一个完整的区域。
        这个区域的左边界是抛物线 $x=\sqrt{y}$，右边界是直线 $y = \frac{a-x}{a-1}$。

    4.  **交换积分次序:**
        我们需要将区域 $D$ 描述成：$y$ 从一个常数到另一个常数， $x$ 从一条曲线到另一条曲线。
        *   **y 的范围:** 区域的最低点是 $y=0$，最高点是 $x=1$ 时的交点，即 $y=1^2=1$。所以 $0 \le y \le 1$。
        *   **x 的范围:** 对于固定的 $y$， $x$ 的左边界是抛物线 $y=x^2$，即 $x=\sqrt{y}$。右边界是直线 $y = \frac{a-x}{a-1}$，需要反解出 $x$：
            $(a-1)y = a-x \implies x = a - (a-1)y$。
            所以 $\sqrt{y} \le x \le a-(a-1)y$。

    5.  **新的积分表达式:**
        $$
        \int_{0}^{1} dy \int_{\sqrt{y}}^{a-(a-1)y} f(x,y) dx
        $$


# 2022（选填部分）

### **一、填空题 (每题 4 分, 共 20 分)**

**1. $z = \sqrt{\log_a(x^2+y^2)} \quad (a>0)$ 的定义域为 D = _______**


*   **解答:**
    1.  $\log_a(x^2+y^2) \ge 0$
    2.  $x^2+y^2 > 0$
    
    对第一个条件，需要分情况讨论 $a$ 的取值：
    *   **当 $a>1$ 时:** $\log_a(x^2+y^2) \ge 0 \implies x^2+y^2 \ge a^0 = 1$。
    *   **当 $0<a<1$ 时:** $\log_a(x^2+y^2) \ge 0 \implies 0 < x^2+y^2 \le a^0 = 1$。(第二个条件 $x^2+y^2>0$ 已经包含在内)
    
    所以定义域为：
    *   **若 $a>1$, D = $\{ (x,y) \mid x^2+y^2 \ge 1 \}$**
    *   **若 $0<a<1$, D = $\{ (x,y) \mid 0 < x^2+y^2 \le 1 \}$**






**2. 设曲线 $L$ 的参数方程表示为 $\begin{cases} x = \phi(t) \\ y = \psi(t) \end{cases} \quad (\alpha \le t \le \beta)$，则弧长元素 $ds = \underline{\sqrt{[\phi'(t)]^2 + [\psi'(t)]^2} dt}$**
*(题目给的范围是 $a \le x \le b$，但参数方程通常用参数范围 $\alpha \le t \le \beta$ 来定义，这里按标准参数方程弧长元素填写)*

*   **思路:** 弧长元素 $ds$ 是由 $dx$ 和 $dy$ 构成的微小直角三角形的斜边，即 $ds = \sqrt{(dx)^2 + (dy)^2}$。在参数方程下，$dx = \phi'(t)dt$, $dy = \psi'(t)dt$。
*   **解答:**
    $ds = \sqrt{(\phi'(t)dt)^2 + (\psi'(t)dt)^2} = \sqrt{([\phi'(t)]^2 + [\psi'(t)]^2)(dt)^2} = \sqrt{[\phi'(t)]^2 + [\psi'(t)]^2} dt$

**3. 设 $I = \int_{0}^{2} dx \int_{x}^{2x} f(x,y) dy$，交换积分次序后，$I = \underline{\int_{0}^{2} dy \int_{y/2}^{y} f(x,y) dx + \int_{2}^{4} dy \int_{y/2}^{2} f(x,y) dx}$**

*   **思路:** 画出积分区域，然后将 X-型区域重新描述为 Y-型区域。
*   **解答:**
    1.  **积分区域 D:**
        $0 \le x \le 2$
        $x \le y \le 2x$
        这个区域是由直线 $y=x$, $y=2x$ 和 $x=2$ 围成的三角形。三个顶点是 $(0,0)$, $(2,2)$ 和 $(2,4)$。
    2.  **交换次序 (Y-型):**
        需要将区域在 $y=2$ 处分割成两部分 $D_1$ 和 $D_2$。
        *   **区域 $D_1$ (下半部分):** $0 \le y \le 2$。此时左边界是 $y=2x \implies x=y/2$，右边界是 $y=x \implies x=y$。
        *   **区域 $D_2$ (上半部分):** $2 \le y \le 4$。此时左边界是 $y=2x \implies x=y/2$，右边界是 $x=2$。
    3.  **新的积分表达式:**
        $$
        I = \int_{0}^{2} dy \int_{y/2}^{y} f(x,y) dx + \int_{2}^{4} dy \int_{y/2}^{2} f(x,y) dx
        $$

**4. 若级数 $\sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n^p}$ 发散，则 $p \le 0$**

*   **思路:**
    这是一个交错级数。
    1.  当 $p>0$ 时，通项 $\frac{1}{n^p} \to 0$ 且单调递减，根据莱布尼茨判别法，级数收敛。
    2.  当 $p \le 0$ 时，令 $p=-q$ ($q \ge 0$)，通项为 $(-1)^{n-1}n^q$。当 $n \to \infty$ 时，通项的绝对值 $|(-1)^{n-1}n^q| = n^q$ 不趋于 0。根据级数收敛的必要条件，级数发散。
*   **解答:** $p \le 0$

**5. 设 $L$ 为取正向的圆周 $x^2+y^2=4$，则曲线积分 $\oint_L y(ye^x+1)dx + (2ye^x-x)dy = \underline{-8\pi}$**

*   **思路:**
    被积函数形式复杂，且积分路径是封闭曲线，首选**格林公式**。
    $\oint_L Pdx + Qdy = \iint_D (\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}) dxdy$

*   **解答:**
    1.  $P(x,y) = y(ye^x+1) = y^2e^x + y$
    2.  $Q(x,y) = 2ye^x - x$
    3.  计算偏导数：
        $\frac{\partial Q}{\partial x} = 2ye^x - 1$
        $\frac{\partial P}{\partial y} = 2ye^x + 1$
    4.  $\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y} = (2ye^x - 1) - (2ye^x + 1) = -2$
    5.  应用格林公式：
        $$
        \oint_L \dots = \iint_D (-2) dxdy
        $$
        其中 $D$ 是由圆周 $x^2+y^2=4$ 围成的圆盘区域。
    6.  $\iint_D (-2) dxdy = -2 \iint_D dxdy = -2 \times (\text{Area of D})$
    7.  圆的面积是 $\pi r^2 = \pi (2^2) = 4\pi$。
    8.  最终结果: $-2 \times 4\pi = -8\pi$。

---

### **二、选择题 (每题 4 分, 共 20 分)**

**1. 二元函数 $z = f(x,y)$ 在 $(x_0,y_0)$ 处可微的充分条件是 ( )**
*(题目问的是充分条件，需要找一个能推出可微的选项)*

*   **(A) $f(x,y)$ 在 $(x_0,y_0)$ 处连续;**
    *   (错误) 连续是可微的必要条件，不是充分条件。例如 $f(x,y)=\sqrt{x^2+y^2}$ 在(0,0)连续但不可微。
*   **(B) $f'_x(x,y), f'_y(x,y)$ 在 $(x_0,y_0)$ 的某邻域内存在;**
    *   (错误) 仅仅存在偏导数不能保证连续，更不能保证可微。
*   **(C) $\Delta z - f'_x(x_0,y_0)\Delta x - f'_y(x_0,y_0)\Delta y$ 当 $\sqrt{(\Delta x)^2+(\Delta y)^2} \to 0$ 时，是无穷小;**
    *   (错误) 这是可微的**定义**。题目问的是充分**条件**。根据定义，$\Delta z - f'_x\Delta x - f'_y\Delta y = o(\rho)$，其中 $\rho=\sqrt{(\Delta x)^2+(\Delta y)^2}$。选项 C 的表述是 $\Delta z - f'_x\Delta x - f'_y\Delta y = o(1)$，比定义弱。
*   **通常的充分条件是: (D) $f'_x(x,y)$ 和 $f'_y(x,y)$ 在 $(x_0, y_0)$ 的某邻域内存在且在该点连续。**
    *   由于图片中 D 选项不完整，但根据标准教材，这才是可微的充分条件。
    
    **根据可见的选项，没有一个是标准的可微充分条件。但如果必须选，最接近的是(B)的加强版。然而，最有可能的情况是选项(D)被遮挡了，内容是“偏导数连续”。**
    *重新审视 (C)*: (C) 的表述 $\Delta z - f'_x(x_0,y_0)\Delta x - f'_y(x_0,y_0)\Delta y$ 是无穷小，就是 $\lim_{\rho\to 0} (\Delta z - f'_x\Delta x - f'_y\Delta y) = 0$。而可微的定义是 $\lim_{\rho\to 0} \frac{\Delta z - f'_x\Delta x - f'_y\Delta y}{\rho} = 0$。后者是更高阶的无穷小，条件更强。所以 C 也不是正确答案。
    
    **结论：这道题很可能正确答案是未显示的选项D，其内容为“函数$f(x,y)$的偏导数$f'_x, f'_y$在点$(x_0,y_0)$的某邻域内存在，且在点$(x_0,y_0)$处连续”。**


好的，这是对试卷第二页内容的解答。

### **二、选择题 (续)**

**(续上一页选择题第1题)**

*   **(D) $\lim_{\substack{\Delta x \to 0 \\ \Delta y \to 0}} \frac{\Delta z - f'_x(x_0,y_0)\Delta x - f'_y(x_0,y_0)\Delta y}{\sqrt{(\Delta x)^2 + (\Delta y)^2}} = 0$**
    *   **分析:** 这个表达式正是二元函数在一点处**可微的定义**。题目问的是**充分条件**，而定义本身是**充要条件**。但在很多选择题的语境下，如果没有更强的充分条件（如偏导数连续）可选，定义本身也会被视为一个选项。在上一页的选项(A)(B)(C)都错误的情况下，(D)是正确的描述。

    **结论:** 结合上一页的分析，选项(D)是可微的定义，因此是正确的。

**2. 设 $f(x,y)$ 在曲线弧 $L$ 上有定义且连续，$L$ 的参数方程为 $\begin{cases} x=\phi(t) \\ y=\psi(t) \end{cases}$ ($\alpha \le t \le \beta$)，其中 $\phi(t), \psi(t)$ 在 $[\alpha, \beta]$ 上具有一阶连续导数，且 $[\phi'(t)]^2+[\psi'(t)]^2 \neq 0$，则曲线积分 $\int_L f(x,y)ds = (\text{ C })$**

*   **思路:** 这是第一类曲线积分（对弧长的积分）的计算公式。需要将所有变量都用参数 $t$ 来表示。
*   **转换关系:**
    *   $x \to \phi(t)$
    *   $y \to \psi(t)$
    *   $f(x,y) \to f(\phi(t), \psi(t))$
    *   弧长微元 $ds = \sqrt{[x'(t)]^2 + [y'(t)]^2} dt = \sqrt{[\phi'(t)]^2 + [\psi'(t)]^2} dt$
*   **选项分析:**
    (A) $\int_\alpha^\beta f(\phi(t),\psi(t))dt$ (错误，缺少了弧长微元中的根号部分)
    (B) (图片部分遮挡，但格式类似A或D，不含根号)
    (C) $\int_\alpha^\beta f(\phi(t),\psi(t))\sqrt{[\phi'(t)]^2+[\psi'(t)]^2}dt$ (**正确**，完全符合公式)
    (D) $\int_\alpha^\beta f(\phi(t),\psi(t))dt$ (这个选项似乎与A重复，可能是抄录或印刷问题)

**3. 设有限闭区域 D 由分段光滑曲线 L 所围成，L 取正向，函数 $P(x,y), Q(x,y)$ 在 D 上具有一阶连续偏导数，则 $\oint_L Pdx+Qdy = (\text{ D })$**

*   **思路:** 这是**格林公式** (Green's Theorem) 的标准形式。
*   **格林公式:** $\oint_L Pdx+Qdy = \iint_D \left(\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right)dxdy$
*   **选项分析:**
    (A) $\iint_D (\frac{\partial P}{\partial y} - \frac{\partial Q}{\partial x})dxdy$ (错误，符号反了)
    (B) $\iint_D (\frac{\partial Q}{\partial y} - \frac{\partial P}{\partial x})dxdy$ (错误，求导变量错了)
    (C) $\iint_D (\frac{\partial P}{\partial x} - \frac{\partial Q}{\partial y})dxdy$ (错误，求导变量和顺序都错了)
    (D) $\iint_D (\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y})dxdy$ (**正确**)

**4. 设 $\lim_{n \to \infty} nu_n = 0$，则 $\sum_{n=1}^\infty u_n$ ( C )**

*   **思路:** 这个条件比级数收敛的必要条件 $\lim_{n \to \infty} u_n=0$ 要强，但它能告诉我们什么？
*   **分析:**
    $\lim_{n \to \infty} nu_n = 0$ 意味着 $u_n$ 比 $1/n$更快地趋向于0。但是，它没有 $u_n$ 比 $1/n^p$ ($p>1$) 更快地趋向于0 那么强。
    *   **考虑 $u_n = \frac{1}{n\ln n}$**:
        $\lim_{n \to \infty} n u_n = \lim_{n \to \infty} \frac{1}{\ln n} = 0$，满足条件。
        但是，级数 $\sum \frac{1}{n\ln n}$ 根据积分判别法是**发散**的。
    *   **考虑 $u_n = \frac{1}{n^2}$**:
        $\lim_{n \to \infty} n u_n = \lim_{n \to \infty} \frac{1}{n} = 0$，满足条件。
        级数 $\sum \frac{1}{n^2}$ 是**收敛**的 (p-级数，p=2>1)。
    *   **考虑交错级数 $u_n = \frac{(-1)^n}{n\ln n}$**:
        $\lim_{n \to \infty} n u_n = 0$ 满足条件。
        级数 $\sum \frac{(-1)^n}{n\ln n}$ 根据莱布尼茨判别法是**收敛**的。
    
    因为我们既能举出收敛的例子，也能举出发散的例子，所以级数的敛散性是**不确定的**。
*   **选项分析:**
    (A) 条件收敛 (不一定)
    (B) 发散 (不一定)
    (C) 不一定 (正确)
    (D) 绝对收敛 (不一定)

**5. 设平面区域 $D: (x-2)^2+(y-1)^2 \le 1$，若 $I_1 = \iint_D (x+y)^2 d\sigma, I_2 = \iint_D (x+y)^3 d\sigma$ 则有 ( A )**

*   **思路:**
    分析被积函数在积分区域 $D$ 上的符号。
*   **步骤:**
    1.  **积分区域 D:** 这是一个以 $(2,1)$ 为圆心，半径为 1 的圆盘。
    2.  **分析 $x+y$ 的符号:**
        在区域 $D$ 内， $x$ 的范围是 $[2-1, 2+1] = [1,3]$，$y$ 的范围是 $[1-1, 1+1] = [0,2]$。
        因此，在整个区域 $D$ 上，$x>0$ 且 $y \ge 0$。
        这意味着 $x+y > 0$ 在整个区域 $D$ 上恒成立（除了可能的边界点(0,0)，但该点不在D内）。
    3.  **比较被积函数:**
        令 $f(x,y) = x+y$。我们知道在 $D$ 上 $f(x,y) > 0$。
        $I_1 = \iint_D f^2 d\sigma$
        $I_2 = \iint_D f^3 d\sigma$
        我们需要比较 $f^2$ 和 $f^3$ 在区域 $D$ 上的大小。
        *   当 $f > 1$ 时，$f^3 > f^2$。
        *   当 $0 < f < 1$ 时，$f^3 < f^2$。
        *   当 $f = 1$ 时，$f^3 = f^2$。
    4.  **分析 $f=x+y$ 的取值范围:**
        $x \in [1,3], y \in [0,2]$。
        $x+y$ 的最小值在点 $(1,0)$ 附近取得，为 $1$。
        $x+y$ 的最大值在点 $(3,2)$ 附近取得，为 $5$。
        更精确地，在圆盘 $(x-2)^2+(y-1)^2 \le 1$ 上，$x+y$ 的最小值在离直线 $x+y=k$ 最近的点，最大值在最远的点。
        圆心 $(2,1)$ 处，$x+y=3$。
        在整个区域 $D$ 上，$x+y \ge 1$ 恒成立 (当 $x=1, y=1$ 时 $x+y=2$; 当 $x=2, y=0$ 时 $x+y=2$)。
        总之，在区域 $D$ 内，$x+y > 1$。
    5.  **比较积分值:**
        因为在区域 $D$ 上 $(x+y) > 1$，所以 $(x+y)^3 > (x+y)^2$ 恒成立。
        根据二重积分的保号性，如果 $g(x,y) > h(x,y)$ 在区域 D 上恒成立，则 $\iint_D g(x,y)d\sigma > \iint_D h(x,y)d\sigma$。
        因此，$I_2 > I_1$。
*   **选项分析:**
    (A) $I_1 < I_2$ (正确)
    (B) $I_1 = I_2$
    (C) $I_1 > I_2$
    (D) 不能比较

---

### **三、计算题 (每题 8 分, 共 24 分)**

**1. 计算 $I = \int_{0}^{2} dx \int_{x}^{2} e^{-y^2} dy$**

*   **思路:**
    被积函数 $e^{-y^2}$ 没有简单的原函数，无法直接对 $y$ 积分。这强烈暗示我们需要**交换积分次序**。

*   **步骤:**
    1.  **画出积分区域 D:**
        $0 \le x \le 2$
        $x \le y \le 2$
        这是一个由 $y=x$, $x=0$ (y轴) 和 $y=2$ 围成的直角三角形。顶点为 $(0,0), (0,2), (2,2)$。
    2.  **交换积分次序 (描述为 Y-型):**
        *   $y$ 的范围是 $[0, 2]$。
        *   对于固定的 $y$， $x$ 的左边界是 $x=0$，右边界是 $x=y$。
        *   所以，新范围是 $0 \le y \le 2, 0 \le x \le y$。
    3.  **建立新积分:**
        $$
        I = \int_{0}^{2} dy \int_{0}^{y} e^{-y^2} dx
        $$
    4.  **计算内层积分 (关于 x):**
        $$
        \int_{0}^{y} e^{-y^2} dx = e^{-y^2} \int_{0}^{y} 1 \, dx = e^{-y^2} [x]_0^y = y e^{-y^2}
        $$
    5.  **计算外层积分 (关于 y):**
        $$
        I = \int_{0}^{2} y e^{-y^2} dy
        $$
        使用换元法，令 $u = -y^2$，则 $du = -2y\,dy \implies y\,dy = -\frac{1}{2}du$。
        $$
        I = \int_{0}^{-4} e^u \left(-\frac{1}{2}du\right) = \frac{1}{2} \int_{-4}^{0} e^u du = \frac{1}{2} [e^u]_{-4}^0 = \frac{1}{2}(e^0 - e^{-4}) = \frac{1}{2}(1 - e^{-4})
        $$

**最终答案:** $\frac{1-e^{-4}}{2}$