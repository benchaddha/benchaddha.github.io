---
title: "Grokking in Modular Addition: Phase-Transition Dynamics and Fourier Circuits in Toy Models"
date: 2026-03-07T11:25:00Z
tags: ["mechanistic-interpretability", "grokking", "experimental-design"]
categories: ["notes"]
draft: false
math: true
---
## Premise
In an environment increasingly shaped by generated outputs, the structure of model behavior is often difficult to inspect directly. This project examines that problem through a specific phenomenon: grokking.

When trained on small algorithmic datasets, neural networks often exhibit a delayed phase transition. They first memorize the training data, reaching high training accuracy while test accuracy remains low. Then, long after loss appears to have plateaued, they abruptly generalize, recovering the underlying generating rule and performing well on unseen examples.

## Testbed
I study this in a 1-layer transformer (TransformerLens) trained on modular addition: $a + b \pmod{113}$. The goal is not only to reproduce the grokking transition, but to intervene on it and identify the mechanisms that support generalization.

### Working hypotheses
**H1 (Stochastic dynamics):** If grokking behaves like an escape process from a “memorization basin,” then injecting gradient noise (e.g., Langevin-style perturbations) should change the time-to-grok distribution. I will evaluate this by comparing time-to-grok across paired seeds under matched training budgets.

**H2 (Circuit dependence):** After grokking, the model may implement modular addition via Fourier-like representations. I will test causal dependence by targeted ablations: selectively suppressing candidate frequency components and measuring whether generalization fails while training-set performance remains comparatively intact.

## Method Notes
This log will record the experimental setup, intervention details, and partial results. I plan to use paired-seed comparisons and time-to-event analysis (time-to-grok) to quantify effects rather than relying on single-run anecdotes.

> Note: The toy transformer was trained locally on an Apple silicon Mac.