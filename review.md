# 局部敏感哈希的下界

王文昊231880275  7.18 

## 解决了什么问题

首先，我要介绍一下我们面临什么样的难题。在计算机科学中有一个很重要的问题被称为最近邻搜索：给定一个点，搜索出和它相聚最近的一个点。在低维度时，问题有很多种解决办法，但是这些办法在维度逐渐升高时，其时空复杂度都会急剧上升，所带来的开销是人难以忍受的。所以我们不去追求精确地解，而是求解他的近似最近邻：可以搜索返回一个距离它不超过cr的点。这样在一定程度上缓急了维度增加带来的骤然升高的开销。

所以我们构造了一个被称为局部敏感哈希（Locality Sensitive Hashing）的工具（后文用LSH简称代替），用于解决近似最近邻问题。该篇论文讨论了局部敏感哈希的下界，从理论上确认了这类函数的最优性能。这个函数是这样定义的，在l1空间中任选一些点，这些点距离一个确定点x的距离小于r时，将会通过哈希函数H以不小于p的概率映射到H(x)上，而距离大于cr的点，将会以不高于q的概率映射到H(x)上。这个哈希函数H被称为局部敏感哈希函数，那么我们再定义一个性能指标ρ=$\frac{log\frac{1}{p}}{log\frac{1}{q}}$,因为搜索的空间复杂度为O(dn+$n^{1+ρ}$)，时间复杂度为O($n^ρ$)+O($n^{ρ}logn$)，所以这个指标的定义是合理的。在本篇论文之前，最好的构造可以达到ρ<=1/c，会不会存在一个更好的构造，让这个指标变得更小的，从理论上进行分析，这个指标的最小值为1/2c，即：无论我们构造出多么优秀的哈希函数，它的性能指标都不会低于1/2c。

## 该问题为什么重要

从理论上讲，LSH是解决高维度下近似最近邻问题最常用，最有效的办法。所以研究它的最优性能自然是一个较为重要的核心问题。只有理解了它的性能极限，我们才能意识到它可能存在的局限性，进而提出更加优秀的办法。

从实践上讲，当我们知道他的最优性能后，我们便有了一个目标方向：当前的最优构造达到了1/c，我们可以尝试进行突破，使其达到1/2c，也可以不再花费时间在这上面，因为最优情况也就是1/2c，并没有在数量级上实现优化。这种确定的最优性能可以帮助我们避免在不可能实现的优化上浪费时间。

同时本文的证明还运用了傅里叶分析和随机游走问题，将此问题进行转化，转化为已知的可解决的问题，这本身也提供了一些证明思路，有助于其他论题的证明。(具体的证明细节在本文末尾)

## 相关研究

有**Indyk 和 Motwani（1998）** ：提出了 LSH 框架，并给出了 l1 空间的上界构造 ρ≤1/c。

**Andoni 和 Indyk（2006）**：对于 l2空间，给出了 ρ≤1/$c^2 $的构造

**Bonami（1970）和 Beckner（1975）**：Bonami-Beckner 不等式是本文证明的核心工具。这个不等式是傅里叶分析中的重要结果，本文将其创造性地应用到了 LSH 下界证明中。

**O‘Donnell、Wu 和 Zhou（2011/2014）** ：这是最重要的后续工作。他们将本文的下界从 1/(2c) 改进到了最优的 1/c（在 Hamming 空间下），这与Indyk 和 Motwani的构造相同，所以ρ的下界精确值被确定为1/c。

## 本文的解决思路

### 先证明一个定理1.3

对任意c,s>=1：
$$
\rho_s(c) \ge \frac{e^{1 / c^s} - 1}{e^{1 / c^s} + 1} \ge \frac{e - 1}{e + 1} \cdot \frac{1}{c^s} \ge \frac{0.462}{c^s}
$$
对于该定理后半部分的不等号很好推理，我们可以设$1/C^S$=x,那么对函数f(x)=$\frac{(e^x-1)}{(e^x+1)}$进行分析，进行二次求导，f‘’(x)=$\frac{2e^x}{(e^x+1)^2}*(1-(e^2)^x)$

由于x>0，所以该二阶导恒小于0，且$c^s$>=1,所以0<=x<=1，那么f(0)=0,f(1)=$\frac{e-1}{e+1}$,所以f(x)恒大于这两点所连接的弦，即f(x)=$\frac{e-1}{e+1}*x$，所以第二个不等号成立，对e取近似值可得到第三个不等号，那么我们只需要证明第一个不等号，后面的自然成立。

那么我们只需要证明第一个等号，就可以推出ρ>=1/2c。

### 为了证明第一个不等号，我们引入以下命题2.1

设 $\mathcal{H}$ 是汉明立方体 $(\{0,1\}^d, \|\cdot\|_1)$ 上的一个 $(r, R, p, q)$-敏感哈希族。假设 $r$ 是奇数且 $R < d/2$，那么
$$
p \le \left( q + e^{-\frac{1}{2}\left(\frac{d}{2} - R\right)^2} \right)^{\frac{e^{2r/d}-1}{e^2r/d+1}}.
$$
取 $R \approx \frac{d}{2} - \sqrt{d \log d}$ 以及 $r \approx R/c$，并令 $d \to \infty$，就得到 $s = 1$ 情形的定理 1.3。对于一般的 $s \ge 1$ 的情形，可由对任意 $x, y \in \{0,1\}^d$ 有 $\|x - y\|_s = \|x - y\|_1^{1/s}$ 这一事实得出。

下面给出命题2.1的证明：

这个命题证明仍旧比较困难，所以我们将他拆分为几个引理

#### **引理 2.2** 

设 $H$ 是汉明立方体 $(\{0,1\}^d, \|\cdot\|_1)$ 上的一个 $(r,R,p,q)$-敏感哈希族，并固定 $x \in \{0,1\}^d$。则
$$
\mathbb{E} \bigl|H^{-1}(H(x))\bigr| \le \sum_{k=0}^{\lfloor R \rfloor} \binom{d}{k} + q \cdot \sum_{k=\lfloor R \rfloor + 1}^{d} \binom{d}{k}.
$$
H(x)表示x映射到的位置，$H^-1$(y)表示映射到位置y上的点，那么$H^-1(H(x))$表示映射位置和x相同的点，在d维汉明立方体中，一个点和另一个点的距离表示为他们两个之间的不同位数，所以固定一个点，距离他为k的点的个数为$\binom{d}{k}$,相当于从d个数中挑出k位进行取反，设距离该点距离大于r小于R的点映射到这个点的概率为p'，那么不超过R的点的个数总和可以表示为$p\sum_{k=0}^{r}\binom{d}{k}+p'\sum_{k=r+1}^{R}\binom{d}{k}+q\sum_{k=R+1}^{d}\binom{d}{k}$,由于p<=1，且p'<=1，所以对p，p'取1，$\mathbb{E} \bigl|H^{-1}(H(x))\bigr|<=$$\sum_{k=0}^{R}\binom{d}{k}+q\sum_{k=R+1}^{d}\binom{d}{k}$

引理2.2证明完毕

根据引理2.2我们还可以进一步得到一个推论

#### **推论 2.3** 

假设 $R < \frac{d}{2}$。我们有
$$
\mathbb{E}\left|\mathcal{H}^{-1}(\mathcal{H}(x))\right| \le 2^d \left( q + e^{-\frac{1}{d}\left(\frac{d}{2} - R\right)^2} \right).
$$
这个推论的证明涉及到一个标准估计
$$
\sum_{k \le \frac{d}{2} - a} \binom{d}{k} \le 2^d \cdot e^{-\frac{a^2}{d}}
$$
此标准估计的推导可见[概率论与数理统计 (Spring 2023)/Entropy and volume of Hamming balls - TCS Wiki](https://tcs.nju.edu.cn/wiki/index.php?direction=next&oldid=11313)

不过此推导给出的是k<=pd的形式，我们令p=1/2-a/d，正好相对应，令c=a/d

并且为了方便求导，将原式中的对数换成以e为底

同时对h(p)在p=1/2处进行二阶泰勒展开可得，h(1/2-c)=ln2+0*c+1/2*(-4)*c^2+O(c^4)=ln2-2*c^2+O(C^4)

即$e^{d \cdot h(1/2-c)}=e^{d*(ln2-2*c^2)}=2^d*e^{-2dc^2}$

再将c=a/d代入可得右边为$2^d*e^{-2a^2/d}<2^d*e^{-a^2/d}$

标准估计推导完毕$$\sum_{k \le \frac{d}{2} - a} \binom{d}{k} \le 2^d \cdot e^{-\frac{a^2}{d}}$$

下面我将进行推论2.3的证明

我们有$$\mathbb{E} \bigl|H^{-1}(H(x))\bigr| \le \sum_{k=0}^{\lfloor R \rfloor} \binom{d}{k} + q \cdot \sum_{k=\lfloor R \rfloor + 1}^{d} \binom{d}{k}.$$由于所有组合数加起来等于2^d，那么右边可以写成$\sum_{k=0}^{\lfloor R \rfloor} \binom{d}{k} + q*(2^d-\sum_{k=0}^{\lfloor R \rfloor} \binom{d}{k})$,化简得，E<=$q*2^d+(1-q)*\sum_{k=0}^{\lfloor R \rfloor} \binom{d}{k}$<=$q*2^d+\sum_{k=0}^{\lfloor R \rfloor} \binom{d}{k}$.

令a=d/2-R代入标准估计，可得

组合数的求和<$2^d (e^{-\frac{1}{d}\left(\frac{d}{2} - R\right)^2} )$加上$q*2^d$可得，E<=$2^d \left( q + e^{-\frac{1}{d}\left(\frac{d}{2} - R\right)^2} \right)$

### 仅有推论2.3还不够，我们还需要引入一个引理2.4

设 $r$ 为奇数。给定非空集合 $B \subseteq \{0,1\}^d$，考虑如下定义的随机变量 $Q_B \in \{0,1\}^d$：从 $B$ 中均匀随机选取一点 $z$，然后从 $z$ 出发在汉明立方体上执行 $r$ 步标准随机游走。将这样得到的点记为 $Q_B$。则
$$
\Pr[Q_B \in B] \le \left( \frac{|B|}{2^d} \right)^{\frac{e^{2r/d} \,-\, 1}{e^{2r/d} \,+\, 1}}.
$$
下面给出它的证明：

**证明。** 我们首先回顾汉明立方体上傅里叶分析的一些背景与记号。
给定 $S \subseteq \{1,\dots , d\}$，沃尔什函数 $W_S : \{0, 1\}^d \to \{-1, 1\}$ 定义为
$$
W_S(u) = (-1)^{\sum_{j\in S} u_j}.
$$
这个函数的含义为，对u的第i位(i属于S)进行统计，统计这些位为1的个数，如果是偶数，则输出为1，如果是奇数则输出为-1



对 $f : \{0, 1\}^d \to \mathbb{R}$，记
$$
\widehat{f}(S) = \frac{1}{2^d} \sum_{u\in\{0,1\}^d} f(u) W_S(u),
$$


计算匹配度，如果f将u映射为正数，且w也将u映射为正数，那么f尖表现为很大的正数，如果相反，表现为很小的负数，如果二者没什么联系，那么求和时相互抵消，最终表现为0附近的数

在本次证明中f是B的指示函数，当u属于B是，f取1，否则取0，所以f尖计算了B中向量，属于S位置的奇偶性

从而 $f$ 可分解为
$$
f = \sum_{S \subseteq \{1,\dots,d\}} \widehat{f}(S) W_S.
$$
对任意 $f, g : \{0, 1\}^d \to \mathbb{R}$，定义内积
$$
\langle f, g\rangle = \frac{1}{2^d} \sum_{u\in\{0,1\}^d} f(u) g(u).
$$
由帕塞瓦尔恒等式，
$$
\langle f, g\rangle = \sum_{S \subseteq \{1,\dots,d\}} \widehat{f}(S) \widehat{g}(S).
$$
帕塞瓦尔恒等式的证明：因为Ws是标准正交基，那么fg可唯一展开为$f = \sum_{S \subseteq \{1,\dots,d\}} \widehat{f}(S) W_S.$

$g = \sum_{T \subseteq \{1,\dots,d\}} \widehat{g}(T) W_T.$

又由于正交性：
$$
\langle W_S, W_T\rangle = \begin{cases} 1, & S = T, \\ 0, & S \neq T. \end{cases}
$$
因此
$$
\sum_{S}\sum_{T} \widehat{f}(S)\,\widehat{g}(T)\,\langle W_S, W_T\rangle = \sum_{S} \widehat{f}(S)\,\widehat{g}(S) = \langle f, g\rangle .
$$


对 $\varepsilon \in [0, 1]$，博纳米–贝克纳算子 $T_\varepsilon$ 定义为
$$
T_\varepsilon f = \sum_{S \subseteq \{1,\dots,d\}} \varepsilon^{|S|} \widehat{f}(S) W_S.
$$
博纳米–贝克纳不等式 [3, 2] 指出，对任意 $f : \{0, 1\}^d \to \mathbb{R}$，
$$
\sum_{S \subseteq \{1,\dots,d\}} \varepsilon^{2|S|} \widehat{f}(S)^2 = \|T_\varepsilon f\|_2^2 = \frac{1}{2^d} \sum_{u\in\{0,1\}^d} (T_\varepsilon f(u))^2 \le \|f\|_{1+\varepsilon^2}^2 = \left( \frac{1}{2^d} \sum_{u\in\{0,1\}^d} f(u)^{1+\varepsilon^2} \right)^{\frac{2}{1+\varepsilon^2}}.
$$
下面给出上式的推理过程：由于$T_\varepsilon f = \sum_{S \subseteq \{1,\dots,d\}} \varepsilon^{|S|} \widehat{f}(S) W_S.$那么，$T_\varepsilon f$在基W下的傅里叶系数为$\varepsilon^{|S|} \widehat{f}(S)$,通过帕塞瓦尔恒等式可得，$\sum_{S \subseteq \{1,\dots,d\}} \varepsilon^{2|S|} \widehat{f}(S)^2 = \|T_\varepsilon f\|_2^2$

通过L2范数的定义$$
\|g\|_2^2 = \mathbb{E}_u [g(u)^2] = \frac{1}{2^d} \sum_u g(u)^2
$$，将g=$T_\varepsilon f$代入，可得第二个等号

第三个不等号是根据博纳米-贝克纳不等式推出（来自deepseek）

第四个等号为范数的展开求和。

将此不等式应用于 $B \subseteq \{0, 1\}^d$ 的示性函数，得到
$$
\sum_{S \subseteq \{1,\dots,d\}} \varepsilon^{2|S|} \widehat{1_B}(S)^2 \le \left( \frac{|B|}{2^d} \right)^{\frac{2}{1+\varepsilon^2}}. \tag{3}
$$
现在，令 $P$ 为 $\{0, 1\}^d$ 上标准随机游走的转移矩阵，即若 $u$ 和 $v$ 恰在一个坐标上不同，则 $P_{uv} = 1/d$，否则 $P_{uv} = 0$。直接计算可知，对每个 $S \subseteq \{1,\dots , d\}$，
$$
P W_S = \left(1 - \frac{2|S|}{d}\right) W_S,
$$
即 $W_S$ 是 $P$ 的特征向量，对应特征值 $1 - \frac{2|S|}{d}$。从 $B$ 中随机一点出发，经过 $r$ 步后仍落在 $B$ 中的概率等于
$$
\begin{aligned}
\Pr[Q_B \in B] &= \frac{1}{|B|} \sum_{a,b\in B} (P^r)_{ab} \\
&= \frac{2^d}{|B|} \langle P^r 1_B, 1_B\rangle \\
&= \frac{2^d}{|B|} \sum_{S \subseteq \{1,\dots,d\}} \widehat{1_B}(S)^2 \left(1 - \frac{2|S|}{d}\right)^r \\
&\le \frac{2^d}{|B|} \sum_{\substack{S \subseteq \{1,\dots,d\} \\ |S| \le d/2}} \widehat{1_B}(S)^2 \left(1 - \frac{2|S|}{d}\right)^r,
\end{aligned}
$$
其中利用了 $r$ 是奇数这一事实（即我们舍去了负项）。
于是，利用 (3) 式可得
$$
\Pr[Q_B \in B] \le \frac{2^d}{|B|} \sum_{S \subseteq \{1,\dots,d\}} \widehat{1_B}(S)^2 \cdot e^{-2r|S|/d} \le \frac{2^d}{|B|} \cdot \left( \frac{|B|}{2^d} \right)^{\frac{2}{1+e^{-2r/d}}} = \left( \frac{|B|}{2^d} \right)^{\frac{1-e^{-2r/d}}{1+e^{-2r/d}}}.
$$
引理2.4证明完毕

### 最后我们来证明命题2.1

对于命题2.1，由问题的定义，映射到同一点的概率至少为p，我们可以得出
$$
p \le \mathop{\mathbb{E}}\limits_{x \in \{0,1\}^d} \Pr\!\big[ H(W_r(x)) = H(x) \big]
= \mathop{\mathbb{E}}\limits_{x} \Pr\!\big[ W_r(x) \in H^{-1}(H(x)) \big].
$$
由全概率公式，上式
$$
= \mathbb{E}_x \sum_{k \in \mathbb{N}} \Pr\!\big[ W_r(x) \in H^{-1}(k) \land H(x) = k \big]
$$
=$Pr[H(x)=k]*Pr[Wr(x)\in H^-1(k)|H(x)=k]$

=$\frac{|H^-1(k)|}{2^d}*Pr[Q_{H^{-1}(k)} \in H^{-1}(k)]$

代入引理2.4
$$
\le \mathbb{E}_x \sum_{k \in \mathbb{N}} \frac{|H^{-1}(k)|}{2^d} \cdot \left( \frac{|H^{-1}(k)|}{2^d} \right)^{\frac{2r/d - 1}{2r/d + 1}}
$$

$$
=  \mathbb{E}_{H}\mathbb{E}_{x \in \{0,1\}^d} \left( \frac{|H^{-1}(H(x))|}{2^d} \right)^{\frac{e^{2r/d} - 1}{e^{2r/d} + 1}}
$$

根据推论2.3
$$
\le \mathbb{E}_{x \in \{0,1\}^d} \left( \frac{\mathbb{E}_H |H^{-1}(H(x))|}{2^d} \right)^{\frac{e^{2r/d} - 1}{e^{2r/d} + 1}}
$$

$$
\le (q + e^{-\frac{1}{d}\left(\frac{d}{2} - R\right)^2}) ^ \frac{e^{2r/d} - 1}{e^{2r/d} + 1}
$$

证明完毕

### 接着我们来完成利用命题2.1推导第一个不等号

当d趋近于无穷是，R/d=1/2，所以2r/d=2R/cd=1/c

而$e^{-\frac{1}{2}(\frac{d}{2}-R)^{2}}=e^{-\frac{dlogd}{2}=d^{-d/2}}=0$

所以命题2.1可以写成$p<=(q+0)^\frac{e^ \frac{1}{c}-1}{e^\frac{1}{c}+1}$

ρ=$\frac{log\frac{1}{p}}{log\frac{1}{q}}>=\frac{e^{1/c}-1}{e^{1/c}+1}*\frac{logq}{logq}=\frac{e^{1/c}-1}{e^{1/c}+1}$

这正是s=1的情况。那么第一个不等号也成立，整个式子成立，ρ>=1/2c。

证明完毕
