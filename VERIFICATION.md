# Verification ledger

The authoritative statement is `manuscript.md` and its build `hyperbolic_slicing_manuscript.pdf` (4 November 2026). Numbering refers to that PDF. Every theorem is supported by a written proof; no proof is machine-verified. The scripts are floating-point sanity checks only (no interval certificates). They use explicit checks rather than `assert`, so they also check under `python -O`, and they exit with an error on failure. Both scripts were run normally and under `-O`; logs are in `results/normal/` and `results/O/`.

| Claim | Proof | Numerical check (script) |
|---|---|---|
| Lemma 2.1, projection formula | Sec. 2.1 | 36 instances against direct minimisation (`verify_and_plot.py`) |
| Lemma 2.2, (2.13)–(2.14), Table 2 | Sec. 2.2–2.3 | six quadratic constants against trigamma (`verify_and_plot.py`) |
| Thm 2.3, completion, diameter | Sec. 2.2 | — (proof only) |
| Prop 3.1, counterexample to BDC 2025 Prop. 19 | Sec. 3 | 12 backprojection evaluations (`verify_and_plot.py`) |
| (4.1), (4.2) | Sec. 4.1 | Figure 1B |
| Thm 4.1, radial rates, (4.7) | Sec. 4.2 | 24 cases over eight $(d,p)$ pairs (`verify_and_plot.py`); near-critical constant $A_{3,1.001}$ (`verify_extensions.py`) |
| Cor 4.2, depth insensitivity | Sec. 4.2 | follows from Thm 4.1 |
| Prop 4.3, dilogarithm formula, $\beta_2=\pi\alpha$ | Sec. 4.3 | seven comparisons (`verify_and_plot.py`) |
| Thm 4.4(ii), $d=3$ series and log law | Sec. 4.3 | Legendre coefficients; series vs. direct quadrature; trend of $\beta_3/(\alpha^2\log(1/\alpha))$ (`verify_extensions.py`) |
| Thm 4.4(iii), $d\geq4$ constant $\kappa_d$ | Sec. 4.3 | $\beta_d(\alpha)/\alpha^2\to\kappa_d$ for $d=4,5,6$ (`verify_extensions.py`) |
| Thm 5.1, Thm 5.2, Cor 5.3 | Sec. 5 | exact binomial identity (5.6), ten sums (`verify_and_plot.py`) |
| Prop 5.4 | Sec. 5.3 | — (proof only) |
| Thm 6.1, Prop 6.2 | Sec. 6 | — (proof only) |

Independent checks during review (not shipped): five AI-assisted reviewers, an adversarial reviewer and a proof audit re-derived every proof and recomputed the constants in higher precision. See `REVIEW_SUMMARY.md`.
