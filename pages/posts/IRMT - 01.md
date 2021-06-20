---
title: 随机矩阵入门
date: 2021/6/21
description: IRMT 第一章
tag: RMT
author: You
---

- # Chapter 1 入门

  对于一个$N$阶方形随机矩阵$H$，若其每个元素都独立抽样自标准正态分布，若将其对称化：
  $$
  H_S=\frac12(H+H^T)
  $$
  我们知道，得到的这个实对称阵$H_S$的特征值就一定是实数。这样的随机矩阵$H_S$我们称其抽样自高斯正交系综（Gaussian Orthogonal Ensemble, GOE）。

  可以用代码生成一个这样的矩阵，以6阶方阵为例：

  ```julia
  using LinearAlgebra
  
  n = 6
  H = rand(n, n)
  Hs = (H + H')/2
  ```

  > julia> H
  > 6×6 Matrix{Float64}:
  > 0.258289   0.918626   0.0594209  0.0798908  0.444571  0.089318
  > 0.736687   0.0993677  0.0850151  0.0453633  0.750808  0.540066
  > 0.562899   0.292135   0.315492   0.773703   0.2393    0.685847
  > 0.233285   0.158351   0.389546   0.0584985  0.861182  0.127726
  > 0.622874   0.507238   0.829647   0.614576   0.419945  0.525829
  > 0.0837589  0.363942   0.327212   0.0013731  0.121012  0.646936
  >
  > julia> Hs
  > 6×6 Matrix{Float64}:
  > 0.258289   0.827657   0.31116   0.156588   0.533722  0.0865384
  > 0.827657   0.0993677  0.188575  0.101857   0.629023  0.452004
  > 0.31116    0.188575   0.315492  0.581624   0.534474  0.50653
  > 0.156588   0.101857   0.581624  0.0584985  0.737879  0.0645497
  > 0.533722   0.629023   0.534474  0.737879   0.419945  0.32342
  > 0.0865384  0.452004   0.50653   0.0645497  0.32342   0.646936

  对这个矩阵做谱分解，得到六个特征值。

  ```julia
  eig_H = eigvals(Hs)
  ```

  > julia> eig_H
  > 6-element Vector{Float64}:
  > -0.4268506734034891
  > -0.3211450762357553
  > -0.029873622037362626
  > 0.4747160387218481
  > 0.9164903934803225
  > 3.497654804051602

  类似地，矩阵$H$的元素推广为复数或四元数时，我们也可以做出类似的定义，此时$H_S$是埃尔米特(hermitian)矩阵或自对偶(self-dual)矩阵，我们称这样的矩阵分别来自高斯酉系综(Gaussian Unitary Ensemble, GUE)和高斯辛系综(Gaussian Symplectic Ensemble, GSE).

  这些特征值的分布存在一定规律。当样本量比较大时，我们可以观察到特征值的规律。我们使用RandomMatrices库中的方法生成随机特征值，观察他们的分布。

  ```julia
  using LinearAlgebra, RandomMatrices, StatsPlots
  
  d_GOE = GaussianHermite(1)
  d_GUE = GaussianHermite(2)
  d_GSE = GaussianHermite(4)
  # N不大时 概率密度存在鞍点(Saddle-Point)
  eig_samp_1 = hcat([vcat([eigvalrand(ensem, 8) for i=1:50000]...) for ensem=[d_GOE,d_GUE,d_GSE]]...)
  density(eig_samp_1, label=["GOE" "GUE" "GSE"])
  
  # N很大时 概率密度呈现半圆律(Wigner’s semicircle law)的形状
  eig_samp_2 = hcat([eigvalrand(ensem, 10000) for ensem=[d_GOE,d_GUE,d_GSE]]...)
  density(eig_samp_2, label=["GOE" "GUE" "GSE"])
  ```

  > ![01-1](D:\读书笔记\随机矩阵\随机矩阵引论-理论与实践-code\01-1.png)
  >
  > ![01-1](D:\读书笔记\随机矩阵\随机矩阵引论-理论与实践-code\01-2.png)

  以GOE为例，很显然$H_S$的非对角线元素$(H_S)_{ij}$都来自正态分布$N(0,1)$，对角线元素$(H_S)_{ii}$都来自$N(0,\frac12)$，故而联合密度为
  $$
  p\left( H_S \right) =\prod_{i=1}^N{\left[ \frac{1}{\sqrt{2\pi}}\exp \left( -\frac{\left( H_S \right) _{ii}^{2}}{2} \right) \right]}\prod_{i<j}{\left[ \frac{1}{\sqrt{\pi}}\exp \left( -\left( H_S \right) _{ij}^{2} \right) \right]}
  $$

