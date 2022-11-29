# Ellipsoids.jl

This is a Julia implementation of [this C++ library](https://github.com/egidioln/ellipsoidInclusion).

This module implements a function that checks the inclusion of one n-ellipsoid in another. For a positive definite matrix $P\succ0\in\mathbb{R}^{n\times n}$ and a vector $c\in\mathbb{R}^{n}$, an *ellipsoid shaped by* $P$ *and cetered at* $c$ is defined as $E(P,c) := \\{x\in\mathbb{R}^{n}:(x-c)^\top P(x-c)\leq 1\\}$.


The method and exmaples implemented in this library are available in [this paper](https://arxiv.org/abs/2211.06237). Please, cite it as:
```
@misc{calbert2022efficient,
  doi = {10.48550/ARXIV.2211.06237},
  url = {https://arxiv.org/abs/2211.06237},
  author = {Calbert, Julien and Egidio, Lucas N. and Jungers, Raphaël M.},
  title = {An Efficient Method to Verify the Inclusion of Ellipsoids},
  publisher = {preprint (arXiv)},
  year = {2022},
}

```


Some examples are still written in Python and some help rewriting them in Julia would be appreciated.