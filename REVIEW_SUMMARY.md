# Summary of the internal review and the revision

The 1 November 2026 version of this paper went through one internal review round on 4 November 2026. Five reviewer roles took part: an editor, a methodology (proof-checking) reviewer, a domain reviewer, a cross-disciplinary reviewer and a devil's advocate. A separate adversarial check and an independent replay of the code were also run. A verification round then re-reviewed the revision and audited every new proof. All reviews were AI-assisted. This was **not** external peer review.

This file lists the objections and how this version answers them. It contains no scores or grades.

## Overall verdict on correctness

No reviewer or audit found an error in a stated result. Every proof was independently re-derived, and the load-bearing constants were recomputed (the exact quadratic constants, the radial rates (4.7), the counterexample (3.1), the minimax constant $C_{d,p}/8$). Four reviewers and the adversarial check read Proposition 19 and Appendix D.3 of Bonet, Drumetz and Courty (JMLR 2025) and confirmed that the paper quotes and locates them correctly.

## Objections and resolutions

| # | Objection | Resolution in this version |
|---|---|---|
| 1 | The contribution is modest relative to its framing; once the Klein formula is known the completion is close to standard. | Plain-language introduction with the practical consequences up front, results-at-a-glance table, conclusion section. New results added (rows 2–5). |
| 2 | The paper says it gives no rate for the unpowered distance, but its own two-point argument does. | New Corollary 5.3: $n^{-1/(2p)}$ is the minimax rate for $D_p$ itself. |
| 3 | "Sharp" overstates the minimax result: the constants differ by a factor of about 20–80. | Now called rate-optimal; the constant ratios are stated (about 21–19 for $p=1$, 82–62 for $p=2$). |
| 4 | Only expected errors are bounded. | New Theorem 5.1(ii): a distribution-free high-probability bound, by bounded differences applied to the distribution-function majorant. |
| 5 | The boundary geometry is known only for $d=2$. | New Theorem 4.4: the boundary-sphere metric in every dimension: $\pi\alpha$ ($d=2$), an exact Legendre series with $\tfrac12\alpha^2\log(1/\alpha)$ behaviour ($d=3$), and bi-Lipschitz to the angle with exact constant $(d-2)/((d-1)(d-3))$ ($d\ge4$). |
| 6 | The practical meaning of saturation for training is not stated. | New Corollary 4.2 (depth insensitivity) and a practical remark; Section 6 is motivated by it. |
| 7 | Classical hyperbolic Radon literature and the source of the failed Euclidean step are not cited; scope of the correction unclear. | Helgason and Boman–Lindskog cited; Proposition 19 quoted; stated that the counterexample also refutes the general Cartan–Hadamard statement; stated that other uses of the property were not audited. |
| 8 | Measurability of the integrand via clipping holds only almost everywhere; the random-projection bound is uniform only for fixed pairs. | Both statements corrected, with the counterexample to simultaneous uniformity. |
| 9 | Missing definitions and steps ($\SW_p$, $\varepsilon_M$ as a supremum, unbiasedness, fixed origin, Poincaré-to-Klein conversion). | All added. |
| 10 | Code: checks used `assert` (silent under `python -O`); the near-critical radial constant helper truncated its integral. | Explicit checks that exit with an error; whole-line integration; a new script checks the extensions and a near-critical constant. Both scripts pass normally and under `-O`. |
| 11 | Verification round: notation clash between projections and distribution functions; measurability details for the bounded-differences step; degenerate rotation plane; minor wording. | All fixed. No correctness issue remained. |

## What remains open

Sharp minimax constants; the horospherical construction and general Cartan–Hadamard manifolds; sampling rates for the radially augmented metric; whether a Gaussian kernel of $D_2$ gives a useful two-sample test on hyperbolic data.
