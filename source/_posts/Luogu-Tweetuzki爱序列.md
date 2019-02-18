---
title: '[Luogu]Tweetuzki爱序列'
tags:
  - 图论
  - DAG
  - DP
categories: 题解
abbrlink: 601
date: 2018-12-12 21:22:40
---

# 前言

题面传送门: [Luogu-5080](https://www.luogu.org/problemnew/show/P5080)

洛咕里面有人发了跟我类似的题解了- -

用这篇题解来纪念我想出来但是打挂成 $50$ 的 $100pts$

# 题目描述

给定序列 $a_1, a_2, \dots a_n$.

在原数列中任取 $k$ 个数, 排列组合成 $b_1,b_2, \dots b_k$, 使得

$$
\forall i \in [1,k), b_i=2b_{i-1}\ or\ b_i=\frac{b_{i-1}}{3}(3|b_{i-1})
$$

$(3|b_{i-1} \Longleftrightarrow b_{i-1}\ \text{mod}\ 3=0)$

求出最大的 $k$, 并输出这个序列.

$n \le 10^5, 3 \le a_i \le 10^{18}$.

`Special Judge`.

<!-- more -->

# 题目分析

很容易想到的是将这个数向可接上的所有数连上一条有向边.

也就是说:

```cpp
for i:=1...n
    if(a[i] * 2存在)
        Link(i, a[i] * 2所在位置)
    if(a[i] % 3 == 0 && a[i] / 3存在)
        Link(i, a[i] / 3所在位置)
```

这样子就得出了一张图.

对其求最长路即可.

但这样就出现了问题:

1. 时空复杂度怎么保证?

见下文中的`时空复杂度分析`.

2. 最长路不是`NP`问题吗?

实际上, 可以证明我们得出的图是一张`DAG`(有向无环图).

证明过程如下:

首先, 图是有向图. 证明其无环即可.

考虑反证法, 假设有向图中有环.

设 $i$ 为环上某一节点, $len$ 为环长度.

则必有:

$$
a_i \cdot (\frac{1}{3})^x \cdot 2^y=a_i \\\\
x,y \in N_+, x+y=len
$$

约去 $a_i$, 整理:

$$
2^y=3^x
$$

对等式两边取以 $2$ 为底的对数, 整理:
$$
\frac{y}{x}= \text{log}_2 3
$$

可知 $\text{log}_2 3$ 为无理数, 又 $x,y \in N_+$, 故 $x,y$ 无解, 与假设矛盾, 故不成立.

证毕.

这样就变成了一个`DAG`上的最长链问题了, 拓扑排序中`DP`一下即可.

# 时空复杂度分析

由刚才的证明, 我们可以进一步得出一个结论:

$a_i$ 重复对答案没有帮助, 可以直接去重.

这个也很好证，因为无环, 一个数无论如何变化都不会变化回自己了, 所以有多少个都没有用.

观察刚刚的伪代码, 因为 $a_i$ 可以去重, 那么每个节点至多连接两条有向边(实际连边还会更少), 所以可知空间复杂度为 $O(n)$;

而`DAG`的最长链~~显然~~是 $O(n)$ 的.

加上前面连边时 $O(n \text{log} n)$ 的复杂度(确认数是否存在需要 $O(\text{log} n)$~~哈希的话当我没说~~)

总时间复杂度为 $O(n \text{log} n)$.

# 实现细节

+ 可以用`std::map`判断数存不存在. 我用的是`std::tr1::unordered_map`.
+ 因边数不确定, 可以用`std::vector`存图.
+ 读入和输出数据都较大, 请避免使用未关同步的`std::cin`和`std::cout`.

# 代码

```cpp
/**********************************************************
 * Author        : EndSaH
 * Email         : hjxhb1@gmail.com
 * Created Time  : 2018-12-12 20:53
 * FileName      : new.cpp
 * Website       : https://endsah.cf
 * *******************************************************/

#include <sys/mman.h>
#include <unistd.h>
#include <cstdio>
#include <cctype>
#include <vector>
#include <queue>
#include <stack>
#include <tr1/unordered_map>

typedef long long LL;

class Istream
{
    char *ipos;

public :
    Istream()
    {
#ifndef ONLINE_JUDGE
        freopen("new.in", "r", stdin);
        freopen("new.out", "w", stdout);
#endif
        ipos = (char*)mmap(NULL, lseek(0, 0, SEEK_END), PROT_READ, MAP_PRIVATE, 0, 0);
    }

    Istream& operator>>(int& n)
    {
        n = 0;
        while(!isdigit(*ipos))
            ++ipos;
        while(n = (n << 3) + (n << 1) + (*ipos++ & 15), isdigit(*ipos));
        return *this;
    }

    Istream& operator>>(LL& n)
    {
        n = 0;
        while(!isdigit(*ipos))
            ++ipos;
        while(n = (n << 3) + (n << 1) + (*ipos++ & 15), isdigit(*ipos));
        return *this;
    }
} in;

char _obuf[1 << 20], _stk[20];
class Ostream
{
    char *opos, *oend, *stkpos;

public :
    Ostream()
    {
        oend = (opos = _obuf) + (1 << 20) - 1;
        stkpos = _stk;
    }

    ~Ostream()
    { fwrite(_obuf, 1, opos - _obuf, stdout); }

    void Putchar(char ch)
    {
        *opos = ch;
        if(opos == oend)
        {
            fwrite(_obuf, 1, 1 << 20, stdout);
            opos = _obuf;
        }
        ++opos;
    }

    Ostream& operator<<(int n)
    {
        do
        {
            *++stkpos = n % 10 ^ 48;
            n /= 10;
        } while(n);
        while(stkpos != _stk)
            Putchar(*stkpos--);
        return *this;
    }

    Ostream& operator<<(LL n)
    {
        do
        {
            *++stkpos = n % 10 ^ 48;
            n /= 10;
        } while(n);
        while(stkpos != _stk)
            Putchar(*stkpos--);
        return *this;
    }

    Ostream& operator<<(char c)
    {
        Putchar(c);
        return *this;
    }
} out;

template<typename _Tp>
inline bool Chkmin(_Tp& x, const _Tp& y)
{ return x > y ? x = y, true : false; }

template<typename _Tp>
inline bool Chkmax(_Tp& x, const _Tp& y)
{ return x < y ? x = y, true : false; }

const int maxN = 1e5 + 2;

int ind[maxN]/*indegree*/, F[maxN], pre[maxN];
LL a[maxN];
std::vector<int> G[maxN];
std::queue<int> q;
std::stack<int> stk;
std::tr1::unordered_map<LL, int> MAP;

int main()
{
    register int n, k = 0, pos/*position*/, u, v;

    in >> n;
    for(register int i = 1; i <= n; ++i)
    {
        in >> a[i];
        MAP[a[i]] = i;
    }
    for(register int i = 1; i <= n; ++i)
    {
        if(MAP.count(a[i] << 1))
        {
            G[i].push_back(MAP[a[i] << 1]);
            ++ind[MAP[a[i] << 1]];
        }
        if(a[i] % 3 == 0 and MAP.count(a[i] / 3))
        {
            G[i].push_back(MAP[a[i] / 3]);
            ++ind[MAP[a[i] / 3]];
        }
    }

    for(register int i = 1; i <= n; ++i)
        if(!ind[i])
        {
            q.push(i);
            F[i] = 1;
        }
    while(!q.empty())
    {
        u = q.front();
        q.pop();
        for(register unsigned int i = 0; i < G[u].size(); ++i)
        {
            v = G[u][i];
            if(Chkmax(F[v], F[u] + 1))
                pre[v] = u;
            --ind[v];
            if(!ind[v])
                q.push(v);
        }
    }
    for(register int i = 1; i <= n; ++i)
        if(Chkmax(k, F[i]))
            pos = i;

    out << k << '\n';
    for(register int i = pos; i; i = pre[i])
        stk.push(i);
    while(!stk.empty())
    {
        out << a[stk.top()] << ' ';
        stk.pop();
    }
    return 0;
}
```

$\huge\text{Thanks for your consideration!}$