---
title: '[CF1138E]Museums Tour'
categories: 题解
tags:
  - 图论
  - 分层图
  - DAG
  - DP
abbrlink: 32916
date: 2019-04-01 09:00:40
---
<script type="text/javascript" src="/js/src/bai.js"></script>

# 前言

[题面传送门](https://codeforces.com/contest/1138/problem/E)

# 题意

给定一个 $n$ 个点 $m$ 条边的**有向**图（无重边自环），定义一次周期为 $d$ 天。

现在你会在一个周期的第一天从 $1$ 号点开始行走，每经过一条边，就会花费一天时间；你**不能**在一个点逗留不走，但是所有点你都可以重复经过。

每个点在一个周期中，会有部分时间开放，而其他时间不开放；注意，**不开放不意味着你不能走到这个点**。

（举个例子，`10100` 表示一个周期的第一天和第三天开放，$2, 4, 5$ 天不开放）

求你能访问到的最多的位于开放时间的点（同样的点只算一次）。

$n, m \le 10^5, d \le 50$

<!-- more -->

# 题解

这题比较套路... 没有灵光一闪的话是根本想不到的...

首先可以发现强连通分量里面是可以乱走一通的，还是有向图，看数据范围也像是 $\cal O(nd)$，基本可以确定是缩点成 $\rm{DAG}$ 后 $\rm{DP}$ 了

但是那个天数不太好处理，缩成点之后以多少天进去可以转换成很多不同的出来的天数，对答案产生的贡献也不同，不加处理的话 $\rm{DP}$ 状态怎么设都是个很大的问题

你会发现 $\rm{DAG + DP}$ 这条路已经差不多走死了，只能往建模上想了

再仔细观察 $d$... $50$ 的数据范围大有用意啊

这种直接建图再处理特别复杂的，如果你对于这个 $d$ 的范围再敏感一点，你就会马上反应出一个名词：

**分层图**！

估计这一个词就可以点醒所有人了。

天数难以表示？那就换个建模方式，用图上的信息来表示！

具体而言，现在我们将一个点 $u$ 拆作 $d$ 个点 $u_{1, \dots , d}$，第 $i$ 个点表示第 $i$ 天位于 $u$ 点；对于原来的一条有向边 $(u, v)$，这样重新连边：$u_1 \longrightarrow v_2, u_2 \longrightarrow v_3, \cdots, u_d \longrightarrow v_1$

若点 $u$ 在第 $i$ 天开放，那么 $u_i$ 点权为 $1$；否则其点权为 $0$。

又有一个小问题：不会算重吗？

实际上，只有当被拆出来的点在同一个强连通分量中，才会算重，这时候就只需要加一次就可以了；

如果不在，就不会有重复的点权计入答案（想想为什么）

~~显然可以证明~~，若 $u _i$ 到 $u _{i'}$ 有一条路径，那么 $u _{i'}$ 到 $u _i$ 也必定有一条路径

也就是说，拆出来的点在一个强连通分量内，就需要去重；不在的话，他们之间必定没有路径可达，肯定不会再重复统计。

然后问题就变为了：求从 $1$ 号点出发的一条点权最大的路径。

缩点成 $\rm DAG\ DP$ 即可（这里我用的是记忆化搜索，建反图拓扑排序也可以）。

# 代码

```cpp
/**********************************************************
 * Author        : EndSaH
 * Email         : hjxhb1@gmail.com
 * Created Time  : 2019-03-31 11:27
 * FileName      : CF1138.cpp
 * Website       : https://endsah.cf
 * *******************************************************/

#include <cstdio>
#include <cstring>
#include <cctype>

#define debug(...) fprintf(stderr, __VA_ARGS__)
#define Debug(s) debug("The message in line %d, Function %s: %s\n", __LINE__, __FUNCTION__, s)
#define getchar() (ipos == iend and (iend = (ipos = _ibuf) + fread(_ibuf, 1, __bufsize, stdin), ipos == iend) ? EOF : *ipos++)
#define putchar(ch) (opos == oend ? fwrite(_obuf, 1, __bufsize, stdout), opos = _obuf : 0, *opos++ = (ch))
#define __bufsize (1 << 21)

char _ibuf[__bufsize], _obuf[__bufsize], _stk[20];
char *ipos = _ibuf, *iend = _ibuf, *opos = _obuf, *oend = _obuf + __bufsize, *stkpos = _stk;

struct END
{ ~END() { fwrite(_obuf, 1, opos - _obuf, stdout); } }
__;

inline int read()
{
    register int x = 0;
    register char ch;
    while (!isdigit(ch = getchar()));
    while (x = (x << 3) + (x << 1) + (ch & 15), isdigit(ch = getchar()));
    return x;
}

template <typename _INT>
inline void write(_INT x)
{
    while (*++stkpos = x % 10 ^ 48, x /= 10, x);
    while (stkpos != _stk)
        putchar(*stkpos--);
}

template <typename _Tp>
inline bool Chkmax(_Tp& x, const _Tp& y)
{ return x < y ? x = y, true : false; }

template <typename _Tp>
inline bool Chkmin(_Tp& x, const _Tp& y)
{ return x > y ? x = y, true : false; }

const int maxN = 5e6 + 2;
const int maxD = 52;

int n, m, d, dfst, colcnt, top;
char opt;
int stk[maxN], dfn[maxN], low[maxN], col[maxN], sum[maxN], F[maxN];
int head[2][maxN];

struct Chain
{ int v, next; }
chain[2][maxN];

inline void Link(int u, int v)
{
    static int ecnt = 0;
    chain[0][++ecnt] = (Chain){v, head[0][u]};
    head[0][u] = ecnt;
}

inline void Newlink(int u, int v)
{
    static int ecnt = 0;
    chain[1][++ecnt] = (Chain){v, head[1][u]};
    head[1][u] = ecnt;
}

void Tarjan(int u)
{
    dfn[u] = low[u] = ++dfst, stk[++top] = u;
    for (int i = head[0][u], v; i; i = chain[0][i].next)
    {
        v = chain[0][i].v;
        if (!dfn[v])
            Tarjan(v), Chkmin(low[u], low[v]);
        else if (!col[v])
            Chkmin(low[u], dfn[v]);
    }
    if (dfn[u] == low[u])
    {
        col[u] = ++colcnt;
        while (stk[top] != u)
            col[stk[top--]] = colcnt;
        --top;
    }
}

int DFS(int u)
{
    if (F[u])
        return F[u];
    for (int i = head[1][u]; i; i = chain[1][i].next)
        Chkmax(F[u], DFS(chain[1][i].v));
    return F[u] += sum[u];
}

int main()
{
#ifndef ONLINE_JUDGE
    freopen("CF1138E.in", "r", stdin);
    freopen("CF1138E.out", "w", stdout);
#endif
    n = read(), m = read(), d = read();
    while (m--)
    {
        int u = read(), v = read();
        for (register int i = 1; i < d; ++i)
            Link(u + (i - 1) * n, v + i * n);
        Link(u + (d - 1) * n, v);
    }
    for (register int i = 1; i <= n * d; ++i)
        if (!dfn[i])
            Tarjan(i);
    for (register int i = 1; i <= n * d; ++i)
        for (register int v, j = head[0][i]; j; j = chain[0][j].next)
            if (col[i] != col[v = chain[0][j].v])
                Newlink(col[i], col[v]);
    for (register int i = 1; i <= n; ++i)
    {
        for (register int j = 0; j < d; ++j)
        {
            while (isspace(opt = getchar()));
            if ((opt & 1) and F[col[i + j * n]] != i)
                F[col[i + j * n]] = i, ++sum[col[i + j * n]];
        }
    }
    memset(F + 1, 0, colcnt * sizeof(int));
    write(DFS(col[1]));
    return 0;
}
```

# Thanks for your consideration!
