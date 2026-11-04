"""Numerical sanity checks for the revision's extensions (not proofs).

Checks: (a) Legendre coefficients of atanh used in Theorem 4.3(ii);
(b) the d=3 boundary series against direct spherical quadrature;
(c) the d>=4 small-angle constant (d-2)/((d-1)(d-3));
(d) the d=3 logarithmic small-angle law (slow convergence, trend only);
(e) the near-critical radial constant A_{3,1.001} by formula (4.4).
Writes data/extension_checks.json (or --out DIR). Fails with SystemExit, also under -O.
"""
import json, sys, platform
from pathlib import Path
import mpmath as mp

mp.mp.dps = 15
ROOT = Path(__file__).resolve().parents[1]
OUT = Path(sys.argv[sys.argv.index("--out")+1]) if "--out" in sys.argv else ROOT/"data"

def check(ok, what):
    if not ok:
        raise SystemExit("CHECK FAILED: " + what)

def cd(d): return mp.gamma(mp.mpf(d)/2)/(mp.sqrt(mp.pi)*mp.gamma(mp.mpf(d-1)/2))

def boundary_sq(d, a):
    """D_2(delta_u, delta_w)^2 for boundary points at angle a, by quadrature."""
    ca, sa = mp.cos(a), mp.sin(a)
    def inner(t):
        r = mp.sqrt(1-t*t)
        if d == 3:
            return mp.quad(lambda ph: (mp.atanh(t)-mp.atanh(t*ca+r*sa*mp.cos(ph)))**2, [0, mp.pi])/mp.pi
        return mp.quad(lambda x: (mp.atanh(t)-mp.atanh(t*ca+r*sa*x))**2*cd(d-1)*(1-x*x)**(mp.mpf(d-4)/2), [-1, 0, 1])
    return mp.re(mp.quad(lambda t: inner(t)*cd(d)*(1-t*t)**(mp.mpf(d-3)/2), [-1, ca, 1]))

def series3(a, L=4001):
    """2*sum_{l odd<L} (2l+1)/(l^2(l+1)^2)*(1-P_l(cos a)), Legendre by recurrence."""
    c = mp.cos(a); p0, p1 = mp.mpf(1), c; tot = mp.mpf(0)
    for l in range(1, L):
        if l % 2 == 1:
            tot += mp.mpf(2*l+1)/(l*l*(l+1)**2)*(1-p1)
        p0, p1 = p1, ((2*l+1)*c*p1 - l*p0)/(l+1)
    return 2*tot

res = {"legendre": [], "d3_series": [], "d_ge4_small_angle": [], "d3_log_ratio": []}
for l in (1, 3, 5, 7, 9):
    v = mp.quad(lambda t: mp.atanh(t)*mp.legendre(l, t), [-1, 0, 1])
    res["legendre"].append([l, float(v), 2.0/(l*(l+1))]); check(abs(v-mp.mpf(2)/(l*(l+1))) < 1e-12, f"legendre {l}")
for a in (0.5, 1.5, float(mp.pi)):
    q, s = boundary_sq(3, a), series3(a)
    res["d3_series"].append([a, float(q), float(s)]); check(abs(q-s) < 5e-4, f"d3 series a={a}")
check(abs(series3(mp.pi, 40001)-mp.pi**2/3) < 1e-8, "d3 diameter pi^2/3")
for d in (4, 5, 6):
    k = mp.mpf(d-2)/((d-1)*(d-3)); r1, r2 = boundary_sq(d, 0.1)/0.01, boundary_sq(d, 0.02)/0.0004
    res["d_ge4_small_angle"].append([d, float(r1), float(r2), float(k)])
    check(abs(r2-k) < abs(r1-k) and abs(r2-k) < 0.01, f"d={d} small-angle constant")
for a in (1e-2, 1e-3, 1e-4):
    ratio = series3(a, 200001)/(a*a*mp.log(1/a))
    res["d3_log_ratio"].append([a, float(ratio)])
check(res["d3_log_ratio"][0][1] > res["d3_log_ratio"][1][1] > res["d3_log_ratio"][2][1] > 0.5, "d3 log ratio decreasing to 1/2")
d, p = 3, mp.mpf("1.001")
k = mp.mpf(d+1)/2  # integrand ~ x^(p-k) at 0: subtract it analytically (near-critical singularity)
I0 = 1/(p-k+1) + mp.quad(lambda x: (mp.log1p(x)**p - x**p)*x**(-k), [0, 1])
A = 2**(d-p-1)*cd(d)*(I0 + mp.quad(lambda x: mp.log1p(x)**p*x**(-k), [1, mp.inf]))
res["A_3_1.001"] = float(A); check(abs(A-mp.mpf("1000.30712861759")) < 1e-6, "A_{3,1.001}")
res["versions"] = {"python": platform.python_version(), "mpmath": mp.__version__}
OUT.mkdir(parents=True, exist_ok=True)
(OUT/"extension_checks.json").write_text(json.dumps(res, indent=2)+"\n")
print(json.dumps({"checks_passed": True, **{k: v for k, v in res.items() if k != "versions"}}, indent=1))
