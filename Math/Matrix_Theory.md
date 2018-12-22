---
title: Matrix Theory
author: yubaoliu89@gmail.com
---

# The Rank of a Matrix

The rank of a matrix is defined as

(a) the maximum number of linearly independent column vectors in the matrix or

(b) the maximum number of linearly independent row vectors in the matrix. Both definitions are equivalent.

1. **Proof**:
若A是可逆矩阵，证明: $r(AB)=r(B)$

**Solution 1:**

$$
r(AB) \leq r(B) \\
B = A^{-1}*AB \\
r(B) = r(A^{-1}*AB) \leq r(AB)
$$

**Solution 2:**
A 可逆，则其可以寫成初等矩陣的乘積:
$$
AB = (A_1 A_2  ... A_s)B
$$
Fact: 左乘初等矩陣相當於進行一次行變換,初等行變換不改變矩陣的秩。

1. **Proof**: 若A是n階冪等矩陣，證明:$r(A)+r(I-A) = n$

**Solution:**
$$
r(A)+r(I-A) \geq r(A+I-A) = R(I) = n \\
r(A(I-A)) = r(A-A^2) = r(A-A) = 0 \\
$$

Fact: $r(AB) \geq r(A) + r(B) -n$

Therefor: $n  \geq r(A) + r(B)$

## 矩陣的等價標準形
1. 如果sxn矩陣A秩爲r, 則A必與
 $$
 M =
 \begin{bmatrix}
I & 0 \\
0 & 0 \\
 \end{bmatrix}
 $$等價。

1. sxn矩陣A,B等價
$$
  \Leftrightarrow r(A) = r(B) \\
  \Leftrightarrow 存在可逆矩陣P,Q使得 B=PAQ
$$

## 矩陣的滿秩分解
1. 假設存在$r(A_{sn}) = r$, 證明存在$B_{sr},$C_{rn}$, Let $A = BC$




















# Resources
- [【东南大学 周建华教授】矩阵论（分析）（全29集](https://www.bilibili.com/video/av14731888/?p=3&t=75)
