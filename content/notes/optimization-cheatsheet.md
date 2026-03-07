---
title: "Optimization Cheat Sheet"
date: 2026-03-03T10:00:00Z
tags: ["optimization", "math"]
categories: ["notes"]
draft: true
math: true
---

For gradient descent on a smooth function, the update is

$$
\theta_{t+1} = \theta_t - \eta \nabla_\theta \mathcal{L}(\theta_t).
$$

A common L2-regularized objective is

$$
\mathcal{L}_{\lambda}(\theta) = \mathcal{L}(\theta) + \lambda \lVert\theta\rVert_2^2.
$$

This page enables KaTeX rendering via `math: true` in front matter.
