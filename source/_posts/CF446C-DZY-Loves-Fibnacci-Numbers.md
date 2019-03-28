---
title: '[CF446C]DZY Loves Fibnacci Numbers'
categories: 题解
tags:
  - 线段树
  - 斐波那契数列
abbrlink: 42251
date: 2019-03-28 14:27:53
---
<script type="text/javascript" src="/js/src/bai.js"></script>

# $\text{Description}$

给定数列 $a_1,a_2,\cdots,a_n$.

需要支持两种操作:

+ `1 l r`: 对 $l$ 到 $r$ 的每个数 $a_i$ 加上 $Fib_{i - l + 1}$.
+ `2 l r`: 询问区间 $l$ 到 $r$ 的和.

$Fib$ 即斐波那契数列,  此处 $Fib_1 = Fib_2 = 1$.

$n \le 3 \ast 10^5$, 答案对 $10^9 + 9$ 取模.

<!-- more -->

# $\text{Solution}$

对于这种区间修改和区间求和, 首先想到线段树.

但是该如何打标记, 这才是最关键的问题.

单单是维护区间和挺简单的:

对于一个修改操作, 记录 $Fib$ 数列的前缀和, 每次对应的加上就可以了.

具体来说, 设当前需要修改的线段树节点对应区间为 $[l, r]$, 修改区间为 $[L, R]$ (~~显然~~, $[l, r] \subseteqq [L, R]$)

那么当前区间应该加上的值是: $sum_{r - L + 1} - sum_{l - L}$

$sum$ 表示 $Fib$ 数列的前缀和.

`Pushdown`操作同理, 只不过是 $L=l, R=r$.

但这么做有个问题: **如何快速合并标记呢**?

我可能是那个头最铁的人, 想到这个思路直接开打, 直接**不合并标记**, 而是用一个`std::stack`存放所有标记, `Pushdown`再一个个下放... 可想而知这样的程序效率是多低, $5 \ast 10^4$ 的数据都跑了 $\text{3s}$ ...

那么, 我们是时候稍微转变一下标记的维护思路了.

在此之前, 你需要先了解[**广义斐波那契数列**](https://baike.baidu.com/item/%E5%B9%BF%E4%B9%89Fibonacci%E6%95%B0%E5%88%97).

通俗的来讲, 就是把 $Fib$ 数列的第一项和第二项变为 $a$ 和 $b$, 而不再是固定的 $1,1$.

它有一些非常优美的性质:

$$\sum_{i = 1}^nS_i=S_{n+2}-S_2$$
$$S_i=bFib_{i - 1}+aFib_{i - 2}$$

$S$ 即为广义 $Fib$ 数列, $S_1 = a, S_2 = b$.

~~证明? 自己想去~~

于是考虑记录广义 $Fib$ 数列的前两项作为懒标记.

为什么? 因为这样就可以**快速合并**了!

> **要合并(对应数各自相加)两个广义 $Fib$ 数列, 只需要将 $a$ 和 $b$ 各自相加即可.**

上面这句话请着重注意!  ~~证明是不存在的~~

因为它的这些性质, 那么合并标记, 快速求和等等一系列问题, 我们已经全部解决了.

接下来, 就是恶心的各类细节操作了- -

详见[代码](https://www.luogu.org/paste/4kxm4tgf).

重点看几个函数: `Add`, `Pushdown`, `Change中return部分`.

# $\text{Thank you for your consideration!}$
