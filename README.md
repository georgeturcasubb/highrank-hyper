# highrank-hyper

Magma scripts for producing infinite families of high-rank hyperelliptic curves.

## Overview

This repository contains a Magma script that generates genus 2 hyperelliptic curves over Q whose Jacobians achieve high Mordell–Weil rank after a small-degree base change. In particular, the script produces examples where, after base change to a number field K with [K:Q] ≤ 8, the Jacobian has rank at least 30.

- Curves: genus 2 hyperelliptic curves C/M, where M is an extension of degree up to 4 of Q. 
- Guarantee: exhibits independent points on Jac(C)(K) to give a lower bound rank(Jac(C)/K) ≥ 30 for some number field K with degree ≤ 8 over Q.
- Note: The rank over M may be smaller; the claim concerns the rank after base change to K.

## Requirements

- Magma (https://magma.maths.usyd.edu.au)
- A machine with sufficient RAM/CPU for high-rank searches (these computations can be intensive).

## Usage

1. Open Magma.
2. Load the script:
   ```magma
   load "highrankgen2.m";
   ```
3. Optional: adjust any parameters in the script.
4. Run the entry point/function as indicated in the script comments to start the search or to reproduce included examples.

## Output

The script reports, for each example it finds:

- A genus 2 curve C given by an equation y^2 = f(x) with f ∈ M[x].
- A number field K/Q with [K:M] <= 2 and [K:Q] ≤ 8 used for the base change.
- A call to the function RankBoudnds(J).
- Additional diagnostics that may include height pairings, independence checks, or saturation information (depending on how the script is configured).

## Reproducibility

- If the search involves randomness, the script may allow setting a seed to reproduce results. Check the top of the script for a parameter like `SetSeed(...)` or a similar option.
- Deterministic reconstructions of specific examples can typically be rerun by loading the script and invoking the function that prints/constructs stored examples.

## Caveats

- High-rank searches are computationally demanding and can take significant time and memory.
- The guaranteed lower bound concerns rank over K with [K:Q] ≤ 8, not necessarily over Q.
- Some verification steps (e.g., height pairing or saturation) can be expensive; you may disable or enable them depending on your needs and resources.

## Citation

If you use or build on these examples in academic work, please cite this repository. When possible, include the specific curve(s) and the degree of the base field used to attain rank ≥ 30.

## Acknowledgments

Thanks to the developers of Magma and the community working on high-rank Jacobians and hyperelliptic curves, whose techniques inspire and inform this search.
