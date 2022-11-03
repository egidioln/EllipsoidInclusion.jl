# EllipsoidInclusion.jl

This is a Julia implementation of [this C++ library](https://github.com/egidioln/ellipsoidInclusion).

This module implements a function that checks the inclusion of one n-ellipsoid in another. For a positive definite matrix $P\succ0\in\mathbb{R}^{n\times n}$ and a vector $c\in\mathbb{R}^{n}$, an *ellipsoid shaped by* $P$ *and cetered at* $c$ is defined as $E(P,c) := \\{x\in\mathbb{R}^{n}:(x-c)^\top P(x-c)\leq 1\\}$.
