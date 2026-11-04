#!/usr/bin/env python3
"""Reproduce deterministic checks and figures for the research candidate.

The mathematical proofs are in manuscript.md. Numerical agreement is not a proof.
Run from any working directory: python code/verify_and_plot.py
Requires Python >=3.10, NumPy, SciPy, and Matplotlib.
"""
from __future__ import annotations

import csv
import json
import math
import os
import platform
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/hyperbolic_slicing_matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import beta, gammaln, polygamma, spence
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "figures"
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)


def check(ok, what: str) -> None:
    """Explicit check that also runs under python -O (no assert)."""
    if not ok:
        raise SystemExit(f"CHECK FAILED: {what}")


def c_dim(d: int) -> float:
    return math.exp(gammaln(d / 2) - 0.5 * math.log(math.pi) - gammaln((d - 1) / 2))


def log_sech(s: float) -> float:
    s = abs(s)
    return math.log(2) - s - math.log1p(math.exp(-2 * s))


def angular_moment(d: int, p: float) -> float:
    return 2 * c_dim(d) * quad(
        lambda s: s**p * math.exp((d - 1) * log_sech(s)),
        0, np.inf, epsabs=2e-12, epsrel=2e-12, limit=250,
    )[0]


def phi(z: complex | float) -> complex:
    # scipy.special.spence(w) = Li_2(1-w). The series fixes the disk branch.
    z = complex(z)
    return complex(spence(1 - z) - 0.25 * spence(1 - z * z))


def dirac_h2_formula(r: float, t: float, alpha: float) -> float:
    q = 1.0 if math.isinf(r) else math.tanh(r / 2)
    s = 1.0 if math.isinf(t) else math.tanh(t / 2)
    return max(0.0, (2 * phi(q*q) + 2 * phi(s*s) - 4 * phi(q*s*np.exp(1j*alpha))).real)


def slice_angle(r: float, theta: float) -> float:
    if math.isinf(r):
        # Values at the two singular angles are not sampled by quadrature.
        return math.log(abs(math.cos(theta / 2))) - math.log(abs(math.sin(theta / 2)))
    a = math.tanh(r) * math.cos(theta)
    return math.atanh(a)


def dirac_h2_quadrature(r: float, t: float, alpha: float) -> float:
    singular = [0.0, math.pi, 2 * math.pi]
    singular += [alpha % (2 * math.pi), (alpha + math.pi) % (2 * math.pi)]
    breaks = sorted(set(singular))
    f = lambda theta: (slice_angle(r, theta) - slice_angle(t, theta - alpha))**2
    return sum(quad(f, a, b, epsabs=1e-10, epsrel=1e-10, limit=250)[0]
               for a, b in zip(breaks[:-1], breaks[1:])) / (2 * math.pi)


def delta(r: float, s: float) -> float:
    # This softplus identity avoids tanh(r) rounding to one for large r.
    return 0.5 * float(np.logaddexp(0.0, 2*(s-r)) - np.logaddexp(0.0, -2*(s+r)))


def rate_constant(d: int, p: float) -> tuple[str, float]:
    m = d - 1
    if m < 2 * p:
        # y=s-r form of the limiting integral on the whole line, stable logarithms.
        def f(y: float) -> float:
            h = 0.5 * float(np.logaddexp(0.0, 2*y))
            if h == 0:
                return 0.0
            return math.exp(p * math.log(h) - m*y)
        val = (quad(f, -np.inf, 0, epsabs=1e-12, epsrel=1e-11, limit=400)[0]
               + quad(f, 0, np.inf, epsabs=1e-12, epsrel=1e-11, limit=400)[0])
        return "boundary", 2**(m+1) * c_dim(d) * val
    if m == 2*p:
        return "critical", 2**(p+1) * c_dim(d)
    return "interior", 2**p * c_dim(d) * beta((p+1)/2, m/2-p)


def scaled_radial_error(d: int, p: float, r: float) -> float:
    m = d - 1
    regime, _ = rate_constant(d, p)
    scale = min(m, 2*p)
    def f(s: float) -> float:
        dv = delta(r, s)
        if dv <= 0:
            return 0.0
        return math.exp(p*math.log(dv) + m*log_sech(s) + scale*r)
    val = 2*c_dim(d)*sum(quad(f, a, b, epsabs=1e-10, epsrel=2e-10, limit=250)[0]
                         for a, b in [(0, r), (r, r+60)])
    return val / r if regime == "critical" else val


def finite_radial_moment(d: int, p: float, r: float) -> float:
    # atanh(tanh(r)tanh(s)) = s-delta(r,s), integrated against angular law.
    return 2*c_dim(d)*quad(lambda s: abs(s-delta(r,s))**p
                          * math.exp((d-1)*log_sech(s)),
                          0, max(40, r+20), epsabs=2e-11, limit=200)[0]


def main() -> None:
    checks: dict = {
        "status": "floating-point diagnostic checks; not formal proofs",
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "normalisation": "curvature -1; normalised angular probability measure",
    }

    # Independent hyperboloid minimisation versus the coordinate formula.
    rng = np.random.default_rng(20261101)
    projection = []
    for d in [2, 3, 7]:
        for _ in range(12):
            u = rng.normal(size=d); u /= np.linalg.norm(u)
            v = rng.normal(size=d); v /= np.linalg.norm(v)
            r = float(rng.uniform(0.05, 4.0))
            x0, a = math.cosh(r), math.sinh(r)*float(u@v)
            expected = math.atanh(a/x0)
            opt = minimize_scalar(lambda t: x0*math.cosh(t)-a*math.sinh(t),
                                  bounds=(-r-1,r+1), method="bounded",
                                  options={"xatol":1e-12})
            projection.append(abs(opt.x-expected))
    checks["projection_formula"] = {"instances": len(projection),
                                      "max_coordinate_error": max(projection)}
    check(max(projection) < 2e-6, 'max(projection) < 2e-6')

    constants = []
    for d in [2,3,4,5,6,10]:
        numeric = angular_moment(d, 2)
        exact = 0.5*float(polygamma(1, (d-1)/2))
        constants.append({"d":d,"quadrature_C_d2":numeric,"trigamma_C_d2":exact,
                          "absolute_error":abs(numeric-exact),"diameter_D2":2*math.sqrt(exact)})
        check(abs(numeric-exact) < 2e-10, 'abs(numeric-exact) < 2e-10')
    checks["quadratic_constants"] = constants

    # Test the published C0 claim using independent integration in d=2 and d=3.
    c0 = []
    for d in [2,3]:
        for r in [0,0.5,1,2,5,10]:
            if d == 2:
                numeric = quad(lambda theta: 1-(math.tanh(r)*math.cos(theta))**2,
                               0,2*math.pi,epsabs=1e-12)[0]/(2*math.pi)
            else:
                numeric = quad(lambda u: 0.5*(1-(math.tanh(r)*u)**2),-1,1,epsabs=1e-12)[0]
            exact = 1-math.tanh(r)**2/d
            c0.append({"d":d,"r":r,"quadrature":numeric,"formula":exact,
                       "absolute_error":abs(numeric-exact)})
            check(abs(numeric-exact) < 1e-11, 'abs(numeric-exact) < 1e-11')
    checks["backprojection_counterexample"] = c0

    fourier = []
    pairs = [(0.4,1.1,0.2),(1.3,2.1,1.2),(3.0,0.6,2.8),
             (np.inf,np.inf,math.pi/6),(np.inf,np.inf,math.pi/2),
             (np.inf,np.inf,math.pi),(1.2,np.inf,0.7)]
    for r,t,alpha in pairs:
        numeric = dirac_h2_quadrature(r,t,alpha)
        formula = dirac_h2_formula(r,t,alpha)
        fourier.append({"r":"infinity" if math.isinf(r) else r,
                        "t":"infinity" if math.isinf(t) else t,
                        "alpha":alpha,"quadrature":numeric,"formula":formula,
                        "absolute_error":abs(numeric-formula)})
        check(abs(numeric-formula) < 3e-8, 'abs(numeric-formula) < 3e-8')
        if math.isinf(r) and math.isinf(t):
            check(abs(formula-math.pi*alpha) < 1e-11, 'abs(formula-math.pi*alpha) < 1e-11')
    checks["dilogarithm_and_boundary_formula"] = fourier

    radial = []
    for d,p in [(2,2),(3,2),(4,2),(5,2),(6,2),(3,1.3),(4,1.5),(8,2.2)]:
        regime, const = rate_constant(d,p)
        for r in [4.0,8.0,12.0]:
            scaled = scaled_radial_error(d,p,r)
            radial.append({"d":d,"p":p,"r":r,"regime":regime,
                           "predicted_coefficient":const,"scaled_error":scaled,
                           "ratio":scaled/const})
    checks["radial_asymptotics"] = radial
    for d,expected in [(2,4*math.log(2)),(3,math.pi**2/6),(5,6.0)]:
        check(abs(rate_constant(d,2)[1]-expected) < 2e-9, 'abs(rate_constant(d,2)[1]-expected) < 2e-9')
    # Check convergence trends, not the invalid demand that a finite-r limit be exact.
    for d,p in [(2,2),(3,2),(4,2),(5,2),(6,2),(3,1.3),(4,1.5),(8,2.2)]:
        rs = [row for row in radial if row["d"]==d and row["p"]==p]
        check(abs(rs[-1]["ratio"]-1) < abs(rs[0]["ratio"]-1) + 1e-7, 'abs(rs[-1]["ratio"]-1) < abs(rs[0]["ratio"]-1) + 1e-7')

    binomial_rows = []
    for n in [1,2,3,4,10,25,100,400,1600,6400]:
        k = np.arange(n+1)
        exact_sum = float(np.sum(np.abs(k/n-0.5)*binom.pmf(k,n,0.5)))
        closed = 0.5*float(binom.pmf((n-1)//2,n-1,0.5))
        eta = 1/(4*math.sqrt(n))
        tv = 0.5*float(np.sum(np.abs(binom.pmf(k,n,0.5+eta)-binom.pmf(k,n,0.5-eta))))
        binomial_rows.append({"n":n,"expected_imbalance":exact_sum,
                              "closed_formula":closed,"sqrt_n_scaled":math.sqrt(n)*exact_sum,
                              "asymptotic_constant":1/math.sqrt(2*math.pi),
                              "testing_total_variation":tv})
        check(abs(exact_sum-closed) < 5e-13, 'abs(exact_sum-closed) < 5e-13')
        check(tv <= 0.5 + 2e-12, 'tv <= 0.5 + 2e-12')
    checks["binomial_and_minimax_checks"] = binomial_rows

    (DATA/"verification_results.json").write_text(json.dumps(checks,indent=2)+"\n")
    with (DATA/"radial_asymptotics.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(radial[0])); w.writeheader(); w.writerows(radial)
    with (DATA/"binomial_errors.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(binomial_rows[0])); w.writeheader(); w.writerows(binomial_rows)

    # Non-standard, high-contrast colours chosen for reproducible vector figures.
    aubergine="#663C63"; copper="#BA673E"; teal="#187F82"; mustard="#A78C28"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
                         "axes.spines.right":False,"axes.titleweight":"bold",
                         "figure.facecolor":"#FCFAF6","axes.facecolor":"#FCFAF6",
                         "savefig.facecolor":"#FCFAF6","axes.labelcolor":"#25222A",
                         "xtick.color":"#38343A","ytick.color":"#38343A"})

    fig,axs=plt.subplots(3,1,figsize=(7.2,10.0),layout="constrained")
    radii=np.linspace(0,8,180)
    d2=np.array([math.sqrt(dirac_h2_formula(float(r),0,0)) for r in radii])
    axs[0].plot(radii,radii,color=copper,lw=2.3,label="Hyperbolic Wasserstein")
    axs[0].plot(radii,d2,color=aubergine,lw=2.5,label="Geodesic sliced")
    axs[0].axhline(math.pi/2,color=aubergine,lw=1,ls="--",alpha=0.65)
    axs[0].annotate(r"limit $\pi/2$",(5.2,math.pi/2),(5.2,2.45),arrowprops={"arrowstyle":"-","color":aubergine},color=aubergine)
    axs[0].set(xlabel="Hyperbolic radius r",ylabel="Distance from the origin",
               title="A. A point escapes to infinity",ylim=(0,8.5))
    axs[0].legend(loc="upper left",frameon=False,fontsize=8)
    n=np.geomspace(1,100,140)
    # r=n²; for large r the exact moment is at its analytic saturation value.
    sat=np.array([dirac_h2_formula(float(x*x),0,0) if x*x<30 else math.pi**2/4 for x in n])
    axs[1].loglog(n,n**3,color=copper,lw=2.3,label=r"$W_2^2=n^3$")
    axs[1].loglog(n,sat/n,color=aubergine,lw=2.4,label=r"$D_2^2\leq\pi^2/(4n)$")
    axs[1].set(xlabel="Sequence index n",ylabel="Squared distance to the origin",
               title="B. Rare distant mass")
    axs[1].legend(frameon=False,fontsize=8)
    ns=np.array([x["n"] for x in binomial_rows])
    vals=np.array([x["sqrt_n_scaled"] for x in binomial_rows])
    axs[2].semilogx(ns,vals,"o-",color=teal,ms=4,lw=1.8,label="Exact binomial sum")
    axs[2].axhline(1/math.sqrt(2*math.pi),color=mustard,ls="--",lw=1.6,label=r"$1/\sqrt{2\pi}$")
    axs[2].set(xlabel="Sample size n",ylabel=r"$\sqrt{n}\,\mathbb{E}|B_n/n-1/2|$",
               title="C. Sharp sampling exponent")
    axs[2].legend(frameon=False,fontsize=8)
    for ax in axs: ax.grid(alpha=0.15)
    fig.savefig(FIG/"geometry_and_sampling.pdf",bbox_inches="tight")
    fig.savefig(FIG/"geometry_and_sampling.png",dpi=190,bbox_inches="tight")
    plt.close(fig)

    fig,axs=plt.subplots(2,1,figsize=(7.2,7.3),layout="constrained")
    for d,color in [(2,aubergine),(3,teal),(5,copper),(6,mustard)]:
        rr=[x for x in radial if x["d"]==d and x["p"]==2]
        axs[0].plot([x["r"] for x in rr],[x["ratio"] for x in rr],"o-",lw=1.8,
                    color=color,label=f"d = {d}"+(" (critical)" if d==5 else ""))
    axs[0].axhline(1,color="#6B6770",ls="--",lw=1)
    axs[0].set(xlabel="Hyperbolic radius r",ylabel="Numerical integral / leading asymptotic",
               title="A. Radial boundary rates",ylim=(0.65,1.03))
    axs[0].legend(frameon=False,fontsize=8)
    aa=np.linspace(0,math.pi,250)
    axs[1].plot(aa,np.sqrt(math.pi*aa),color=aubergine,lw=2.5,label=r"$D_2=\sqrt{\pi\alpha}$")
    axs[1].plot(aa,aa,color=teal,lw=1.4,ls="--",label="Angular distance")
    axs[1].set(xlabel="Smaller angular separation α (radians)",ylabel="Boundary distance",
               title="B. Exact boundary circle geometry")
    axs[1].set_xticks([0,math.pi/2,math.pi],["0","π/2","π"])
    axs[1].legend(frameon=False,fontsize=9)
    for ax in axs: ax.grid(alpha=0.15)
    fig.savefig(FIG/"boundary_rates.pdf",bbox_inches="tight")
    fig.savefig(FIG/"boundary_rates.png",dpi=190,bbox_inches="tight")
    plt.close(fig)

    print(json.dumps({"checks_passed":True,"projection_instances":len(projection),
                      "max_projection_error":max(projection),
                      "max_constant_error":max(x["absolute_error"] for x in constants),
                      "max_dilogarithm_error":max(x["absolute_error"] for x in fourier),
                      "radial_cases":len(radial),"figures":2,
                      "data_file":str(DATA/"verification_results.json")},indent=2))


if __name__ == "__main__":
    main()
