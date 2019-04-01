---
title: FFT简洁笔记
categories: 知识点
tags: 数论
abbrlink: 41766
date: 2019-03-24 18:56:27
---
<script type="text/javascript" src="/js/src/bai.js"></script>

# 前言

主要学自 [Menci 的博客](https://oi.men.ci/fft-notes/)

**不讲 $\text{FFT}$，这里只是讲一些其他的东西以及一些例题**

<!-- more -->

# 定义相关

## 多项式乘法本质

$$
C _i = \sum _{j = 0} ^i A _j B _{i - j}
$$

其中 $A, B, C$ 表示多项式的系数集合，$i$ 表示为第 $i$ 次项。

若没有第 $i$ 次项，系数视为 $0$。

## 复数乘法

$$
\begin{aligned}
& (a + bi)(c + di) \\\\
=& ac + adi + bci + bdi^2 \\\\
=& ac + (ad + bc)i - bd \\\\
=& (ac - bd) + (ad + bc)i
\end{aligned}
$$

## 单位根

$$
\omega _n ^k = \cos (\frac{2 k \pi} n) + \sin (\frac{2 k \pi} n)i \\\\
\text{(From: } e^{xi} = \cos x + \sin x i \text{)}
$$

### 两个性质

$$
\begin{aligned}
\omega ^{2k} _{2n} &= \cos (\frac{4 k \pi} {2n}) + \sin (\frac{4 k \pi} {2n})i \\\\
&= \cos (\frac{2 k \pi} n) + \sin(\frac{2 k \pi} n)i \\\\
&= \omega ^k _n \\\\
\omega ^{k + \frac n 2} _n &= \cos (\frac{2 (k + \frac n 2) \pi} n) + \sin (\frac{2 (k + \frac n 2) \pi} n)i \\\\
&= \cos (\frac{2 k \pi + n \pi} n) + \sin (\frac{2 k \pi + n \pi} n)i \\\\
&= -\cos (\frac{2 k \pi} n) + -\sin (\frac{2 k \pi} n)i \\\\
&= -\omega _n ^k
\end{aligned}
$$

# 过程图

![](https://ws1.sinaimg.cn/large/006Ftn1Tgy1g1mujanv8zj30ix07q75w.jpg)

（其中求值和插值朴素算法均为 $\mathcal O(n^2)$，用 $\text{FFT}$ 优化后为 $\mathcal O(n \log n)$）

# 代码实现

**注意**：$n$ 为 $2$ 的方幂。

递归比较慢，考虑如何迭代。

## 雷德算法

发现在分治到底划分完毕后，每个数刚好位于其二进制翻转后的位置：

```
000 001 010 011 100 101 110 111
 0   1   2   3   4   5   6   7
 0   2   4   6 - 1   3   5   7
 0   4 - 2   6 - 1   5 - 3   7
 0 - 4 - 2 - 6 - 1 - 5 - 3 - 7
000 100 010 110 001 101 011 111
```

于是有了如下代码（很抱歉，我并不知道如何证明正确性）：

```cpp
for (lim = 1, l = 0; lim <= n + m; lim <<= 1, ++l);
for (register int i = 0; i < n; ++i)
    rev[i] = rev[i >> 1] >> 1 | (i & 1) << (l - 1); // mind the priority of operation

// when calculate *a:
for (register int i = 0; i < n; ++i)
    if (i < rev[i])
        std::swap(a[i], a[rev[i]]);
```

这种倒序位实现迭代 $\text{FFT}$ 的算法称为**雷德算法**。

于是便可以迭代求解：

```cpp
inline void DFT(Comp* a, int n, int flag) // Comp: Complex
{
    for (register int i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);
    for (register int nn = 2, m; nn <= n; nn <<= 1) // nn: the length
    {
        m = nn >> 1;
        for (register Comp* curr = a; curr != a + n; curr += nn) // combining the A_1 and A_2
            for (register int i = 0; i < m; ++i)
            {
                static Comp tmp;
                tmp = Omega(nn, flag * i) * curr[i + m]; // attention!
                curr[i + m] = curr[i] - tmp, curr[i] += tmp;
            }
    }
}
```

注意到中间的这段代码：

```cpp
tmp = Omega(nn, flag * i) * curr[i + m];
curr[i + m] = curr[i] - tmp, curr[i] += tmp;
```

其中 `flag` 表示其为 $\text{DFT}$ ( `1` ) 还是 $\text{IDFT}$ ( `-1` )。

假设现在正在进行 $\text{DFT}$，那么代码抽出来看就是：

$$  
\begin{aligned}
t &\gets \omega _n ^k \times a _{k + \frac n 2} \\\\
a _{k + \frac n 2} &\gets a_k - t \\\\
a _k &\gets a_k + t
\end{aligned}
$$

省去了之前需要复制数组的麻烦。

这一过程被称为**蝴蝶操作**。

# 例题

## [【模板】多项式乘法](https://www.luogu.org/problemnew/show/P3803)

略。

我跑了大概 5000ms...慢的很

部分代码：

```cpp
const int maxN = 1 << 21 | 1;
const double pi = acos(-1);

int n, m, lim, l;
int rev[maxN];

struct Comp
{
    double real, imag;

    Comp() { }

    Comp(double real, double imag) : real(real), imag(imag) { }

    Comp operator-(const Comp& x) const
    { return Comp(real - x.real, imag - x.imag); }

    Comp operator+(const Comp& x) const
    { return Comp(real + x.real, imag + x.imag); }

    Comp& operator+=(const Comp& x) 
    {
        *this = *this + x;
        return *this;
    }

    Comp operator*(const Comp& x) const
    { return Comp(real * x.real - imag * x.imag, real * x.imag + imag * x.real); }

    Comp& operator*=(const Comp& x)
    {
        *this = *this * x;
        return *this;
    }
} a[maxN], b[maxN];

inline Comp Omega(int n, int k)
{ return Comp(cos((k << 1) * pi / n), sin((k << 1) * pi / n)); }

inline void DFT(Comp* a, int n, int flag)
{
    for (register int i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);
    for (register int nn = 2, m; nn <= n; nn <<= 1)
    {
        m = nn >> 1;
        for (register Comp* curr = a; curr != a + n; curr += nn)
            for (register int i = 0; i < m; ++i)
            {
                static Comp tmp;
                tmp = Omega(nn, flag * i) * curr[i + m];
                curr[i + m] = curr[i] - tmp, curr[i] += tmp;
            }
    }
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("LGP3803.in", "r", stdin);
    freopen("LGP3803.out", "w", stdout);
#endif
    n = read(), m = read();
    for (register int i = 0; i <= n; ++i)
        a[i] = Comp(read(), 0);
    for (register int i = 0; i <= m; ++i)
        b[i] = Comp(read(), 0);
    for (lim = 1, l = 0; lim <= n + m; lim <<= 1, ++l);
    for (register int i = 0; i < lim; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (l - 1);
    DFT(a, lim, 1), DFT(b, lim, 1);
    for (register int i = 0; i < lim; ++i)
        a[i] *= b[i];
    DFT(a, lim, -1);
    for (register int i = 0; i <= n + m; ++i)
        write(int(a[i].real / lim + 0.5)), putchar(' ');
    return 0;
}
```

## [[ZJOI2014]力](https://www.luogu.org/problemnew/show/P3338)

### 题意

~~自己看，特别清晰明了~~

### 题解

~~省选考裸题了解一下~~

设下标从 $0$ 开始。

颓柿子：

首先，$E_j$ 的表达式：

$$ 
E_j = \frac { \sum\limits _{i = 0} ^ {j - 1} \frac {q_i q_j} {(i - j)^2} - \sum\limits _{i = j + 1} ^ {n - 1} \frac {q_i q_j} {(i - j)^2} } {q_j}
$$

上下同时约去 $q_j$：

$$ 
\begin{aligned}
E_j
&= \frac { \sum\limits _{i = 0} ^ {j - 1} \frac {q_i q_j} {(i - j)^2} - \sum\limits _{i = j + 1} ^ {n - 1} \frac {q_i q_j} {(i - j)^2} } {q_j} \\\\
&= \sum _{i = 0} ^{j - 1} \frac {q_i} {(i - j)^2} - \sum _{i = j + 1} ^{n - 1} \frac {q_i} {(i - j)^2}
\end{aligned}
$$

令 $g_i = \frac 1 {i ^ 2}$，特别地，$g_0 = 0$。显然其为偶函数。

则上式可以变为：

$$ 
\begin{aligned}
E_j
&= \sum _{i = 0} ^{j - 1} \frac {q_i} {(i - j)^2} - \sum _{i = j + 1} ^{n - 1} \frac {q_i} {(i - j)^2} \\\\
&= \sum _{i = 0} ^{j} {q_i} g _{j - i} - \sum _{i = j} ^{n - 1} {q_i} g _{i - j}
\end{aligned}
$$

（$j - 1$ 变为 $j$ 是为了形式变得和多项式乘法一样，实际上 $j = i$ 时因 $g_0 = 0$ 不会对答案产生贡献~~才不是我打错了~~）

回忆下开头说的多项式乘法本质，可以发现减号前面的式子就是卷积形式了；至于后面，可以用一张图来感性理解：

![来自GoldenPotato dalao的博客](https://s2.ax1x.com/2019/04/01/As9wLQ.md.jpg)

（来自[GoldenPotato dalao的博客](https://www.cnblogs.com/GoldenPotato/p/10288217.html)）

如果要理性理解的话，请继续看：

考虑将 $q$ 反向变为 $q'$，即令 $q' _{n - i - 1} = q _i$。

将减号两边拆开算，即令 $A_j = \sum\limits _{i = 0} ^{j} {q_i} g _{j - i}, B_j = \sum\limits _{i = j} ^{n - 1} {q_i} g _{i - j}$。

考虑单独计算 $B_j$，此时用 $q'$ 代替 $q$：

$$
\begin{aligned}
B_j
&= \sum _{i = j} ^{n - 1} {q_i} g _{i - j} \\\\
&= \sum _{i = j} ^{n - 1} {q' _{n - i - 1}} g _{i - j}
\end{aligned}
$$

改变和式的枚举范围：

$$
\begin{aligned}
B_j
&= \sum _{i = j} ^{n - 1} q' _{n - i - 1} g _{i - j} \\\\
&= \sum _{i = 0} ^{n - j - 1} q' _{n - j - i - 1} g _i
\end{aligned}
$$

（如果不明白为什么 $g _{i - j}$ 变为 $g_i$，$q' _{n - i - 1}$ 却变为 $q' _{n - j - i - 1}$的话请自己认真好好想一想，两个式子 $i$ 的符号是不一样的）

再令 $t = n - j - 1$，则 $j = n - t - 1$

并且令 $B'$ 为 $B$ 的反向，即 $B' _{n - i - 1} = B _i$

上式可再次变为：

$$
\begin{aligned}
B _j
&= \sum _{i = 0} ^{n - j - 1} q' _{n - j - i - 1} g _i \\\\
B _{n - t - 1}
&= \sum _{i = 0} ^t q' _{t - i} g _i \\\\
B' _t
&= \sum _{i = 0} ^t q' _{t - i} g _i \\\\
\end{aligned}
$$

可以看到，再次变为多项式卷积形式，求得 $B'$ 后反向即可得到 $B$。

$t$ 跟 $n$ 一样就已经可以求出正确的 $B$ 了，再卷出 $A$ 相减即可。

部分代码：

```cpp
const int maxN = 1e5 + 2;
const int maxM = 1 << 18 | 1;
const double pi = acos(-1);

int n, lim, log_n;
double q[maxN], _q[maxN], g[maxN];
int rev[maxM];

struct Comp
{
    double real, imag;

    Comp() { }

    Comp(double real, double imag) : real(real), imag(imag) { }

    Comp operator-(const Comp& x) const
    { return Comp(real - x.real, imag - x.imag); }

    Comp operator+(const Comp& x) const
    { return Comp(real + x.real, imag + x.imag); }

    Comp operator*(const Comp& x) const
    { return Comp(real * x.real - imag * x.imag, imag * x.real + real * x.imag); }

    Comp& operator+=(const Comp& x)
    {
        *this = *this + x;
        return *this;
    }

    Comp& operator*=(const Comp& x) 
    {
        *this = *this * x;
        return *this;
    }
};

inline Comp Omega(int n, int k)
{ return Comp(cos((k << 1) * pi / n), sin((k << 1) * pi / n)); }

inline void Init()
{
    for (lim = 1, log_n = 0; lim < (n << 1); lim <<= 1, ++log_n);
    for (register int i = 1; i < lim; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (log_n - 1);
}

inline void DFT(Comp* a, int n, int flag)
{
    for (register int i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);
    for (register int nn = 2, m; nn <= n; nn <<= 1)
    {
        m = nn >> 1;
        for (register Comp* curr = a; curr != a + n; curr += nn)
            for (register int i = 0; i < m; ++i)
            {
                static Comp tmp;
                tmp = Omega(nn, flag * i) * curr[i + m];
                curr[i + m] = curr[i] - tmp, curr[i] += tmp;
            }
    }
}

inline void Mul(double* a, double* b, double* ans)
{
    static Comp _a[maxM], _b[maxM];
    for (register int i = 0; i < n; ++i)
        _a[i] = Comp(a[i], 0), _b[i] = Comp(b[i], 0);
    for (register int i = n; i < lim; ++i)
        _a[i] = _b[i] = Comp(0, 0);
    DFT(_a, lim, 1), DFT(_b, lim, 1);
    for (register int i = 0; i < lim; ++i)
        _a[i] *= _b[i];
    DFT(_a, lim, -1);
    for (register int i = 0; i < n; ++i)
        ans[i] = _a[i].real / lim;
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("LGP3338.in", "r", stdin);
    freopen("LGP3338.out", "w", stdout);
#endif
    scanf("%d", &n);
    for (register int i = 0; i < n; ++i)
    {
        if (i)
            g[i] = 1.0 / i / i;
        scanf("%lf", q + i), _q[i] = q[i];
    }
    std::reverse(_q, _q + n), Init();
    Mul(q, g, q), Mul(_q, g, _q);
    for (register int i = 0; i < n; ++i)
        printf("%lf\n", q[i] - _q[n - i - 1]);
    return 0;
}
```

~~有没有觉得代码很好看~~

## [【模板】A * B Problem](https://www.luogu.org/problemnew/show/P1919)

将大整数的每一位视作多项式系数，然后再卷积即可。

代码：

```py
input(); print(int(input()) * int(input()))
```

~~我偏不放 `C++` 代码~~

## [BZOJ 2194](https://www.lydsy.com/JudgeOnline/problem.php?id=2194)

### 题意

~~Please contact lydsy2012@163.com!~~

给定一个数 $n$ 以及长度为 $n$ 的序列 $a, b$，要求计算

$$
c_k = \sum _{i = k} ^{n - 1} a _i b _{i - k}
$$

数组下标从 0 开始。

$n \le 10^5$，$\forall a _i, b _i \in N$ 且 $\forall a _i, b _i \le 100$。

### 题解

看了上面推的式子，这边看起来就特别裸了...

直接颓柿子了：

$$
\begin{aligned}  
c _k
&= \sum _{i = k} ^{n - 1} a _i b _{i - k} \\\\
&= \sum _{i = 0} ^{n - k - 1} a _{i + k} b _i \\\\
c' _{n - k - 1}
&= \sum _{i = 0} ^{n - k - 1} a' _{n - i - k - 1} b _i \\\\
c' _{t}
&= \sum _{i = 0} ^t a' _{t - i} b _i
\end{aligned}
$$

然后就没有然后了。

部分代码：

```cpp
const int maxN = 1e5 + 2;
const int maxM = 1 << 18 | 1;
const double pi = acos(-1);

int n, lim, log_n;
int a[maxN], b[maxN], rev[maxM];

struct Comp
{
    double real, imag;

    Comp() { }

    Comp(double real, double imag) : real(real), imag(imag) { }

    Comp operator+(const Comp& x) const
    { return Comp(real + x.real, imag + x.imag); }

    Comp& operator+=(const Comp& x)
    {
        *this = *this + x;
        return *this;
    }

    Comp operator-(const Comp& x) const
    { return Comp(real - x.real, imag - x.imag); }

    Comp operator*(const Comp& x) const
    { return Comp(real * x.real - imag * x.imag, imag * x.real + x.imag * real); }

    Comp& operator*=(const Comp& x)
    {
        *this = *this * x;
        return *this;
    }
};

inline Comp Omega(int n, int k)
{ return Comp(cos((k << 1) * pi / n), sin((k << 1) * pi / n)); }

inline void Init()
{
    for (lim = 1, log_n = 0; lim < (n << 1); lim <<= 1, ++log_n);
    for (register int i = 0; i < lim; ++i)
        rev[i] = rev[i >> 1] >> 1 | (i & 1) << (log_n - 1);
}

inline void DFT(Comp* a, int n, int flag)
{
    for (register int i = 0; i < n; ++i)
        if (i < rev[i])
            std::swap(a[i], a[rev[i]]);
    for (register int nn = 2, m; nn <= n; nn <<= 1)
    {
        m = nn >> 1;
        for (register Comp* curr = a; curr != a + n; curr += nn)
            for (register int i = 0; i < m; ++i)
            {
                static Comp tmp;
                tmp = Omega(nn, i * flag) * curr[i + m];
                curr[i + m] = curr[i] - tmp, curr[i] += tmp;
            }
    }
}

inline void Mul(int* a, int* b, int* ans)
{
    static Comp _a[maxM], _b[maxM];
    for (register int i = 0; i < n; ++i)
        _a[i] = Comp(a[i], 0), _b[i] = Comp(b[i], 0);
    for (register int i = n; i < lim; ++i)
        _a[i] = _b[i] = Comp(0, 0);
    DFT(_a, lim, 1), DFT(_b, lim, 1);
    for (register int i = 0; i < lim; ++i)
        _a[i] *= _b[i];
    DFT(_a, lim, -1);
    for (register int i = 0; i < n; ++i)
        ans[i] = int(_a[i].real / lim + 0.5);
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("BZOJ2194.in", "r", stdin);
    freopen("BZOJ2194.out", "w", stdout);
#endif
    n = read();
    for (register int i = 0; i < n; ++i)
        a[n - i - 1] = read(), b[i] = read();
    Init();
    Mul(a, b, a);
    for (register int i = n - 1; ~i; --i)
        write(a[i]), putchar('\n');
    return 0;
}
```

以后可能还会补充例题的- -

# Thanks for your consideration!
