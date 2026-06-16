#!/usr/bin/env python3
"""Render the primality-test cost figures from data/*.csv to figures/*.png + *.svg.
Pure matplotlib (no pandas). Run:  python plot.py
"""
import csv, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "data"
OUT  = "figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 160, "savefig.dpi": 160, "font.size": 11,
    "axes.grid": True, "grid.alpha": 0.3, "axes.titlesize": 12,
    "legend.framealpha": 0.9, "figure.autolayout": True,
})

def load(fname):
    with open(os.path.join(DATA, fname), newline="") as f:
        return list(csv.DictReader(f))

def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.join(OUT, name + ".png"))

METHOD_STYLE = {
    "Fermat":            ("#d62728", "o"),
    "Solovay-Strassen":  ("#2ca02c", "s"),
    "Miller-Rabin":      ("#1f77b4", "^"),
    "AKS":               ("#9467bd", "D"),
    "Lucas-Lehmer":      ("#ff7f0e", "v"),
    "Trial division":    ("#7f7f7f", "x"),
}

def avg_by_bits(rows, method, col):
    acc = collections.defaultdict(list)
    for r in rows:
        if r["method"] == method and r[col] not in ("", "-"):
            acc[int(r["bits"])].append(float(r[col]))
    xs = sorted(acc)
    ys = [sum(acc[b]) / len(acc[b]) for b in xs]
    return xs, ys

# ---------------- Figure 1: scaling ----------------
def fig1():
    rows = load("fig1_scaling.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for m in ("Fermat", "Solovay-Strassen", "Miller-Rabin"):
        c, mk = METHOD_STYLE[m]
        x, y = avg_by_bits(rows, m, "word_muls")
        ax1.plot(x, y, marker=mk, color=c, label=m, lw=1.6, ms=5)
        xm, ym = avg_by_bits(rows, m, "mod_muls")
        xs, ys = avg_by_bits(rows, m, "mod_sqrs")
        tot = [a + b for a, b in zip(ym, ys)]
        ax2.plot(xm, tot, marker=mk, color=c, label=m, lw=1.6, ms=5)
    # cubic reference line on left panel
    x, y = avg_by_bits(rows, "Miller-Rabin", "word_muls")
    if x:
        x0, y0 = x[-1], y[-1]
        ref = [y0 * (b / x0) ** 3 for b in x]
        ax1.plot(x, ref, "k--", lw=1, alpha=0.6, label=r"$\propto \mathrm{bits}^3$")
    ax1.set(xscale="log", yscale="log", xlabel="bit length of n", ylabel="word multiplications")
    ax1.set_title("(a) word multiplications")
    ax2.set(xscale="log", yscale="log", xlabel="bit length of n", ylabel="modular mults + squarings")
    ax2.set_title("(b) modular mults + squarings")
    ax1.legend(fontsize=9); ax2.legend(fontsize=9)
    save(fig, "fig1_scaling")

# ---------------- Figure 2: Fermat fooled ----------------
def fig2():
    rows = load("fig2b_fermat_seeds.csv")
    fooled = collections.Counter(); total = collections.Counter()
    for r in rows:
        total[r["method"]] += 1
        if r["verdict"] == "PROBABLY_PRIME":
            fooled[r["method"]] += 1
    methods = ["Fermat", "Miller-Rabin"]
    vals = [fooled[m] for m in methods]
    tot  = [total[m] for m in methods]
    fig, ax = plt.subplots(figsize=(6, 4.2))
    bars = ax.bar(methods, vals, color=[METHOD_STYLE[m][0] for m in methods], width=0.55)
    for b, v, t in zip(bars, vals, tot):
        ax.text(b.get_x() + b.get_width()/2, v + 0.4, f"{v}/{t}", ha="center", fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(tot) + 3)
    ax.set_ylabel("times declared PROBABLY PRIME (fooled)")
    save(fig, "fig2_fermat_fooled")

# ---------------- Figure 3: cost vs k / confidence ----------------
def fig3():
    rows = load("fig3_confidence.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    for m in ("Fermat", "Solovay-Strassen", "Miller-Rabin"):
        c, mk = METHOD_STYLE[m]
        pts = sorted((int(r["k"]), float(r["word_muls"])) for r in rows if r["method"] == m)
        ax1.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, color=c, label=m, lw=1.6)
    ax1.set(xlabel="rounds k", ylabel="word multiplications")
    ax1.set_title("(a) cost vs rounds k")
    ax1.legend(fontsize=9)
    # equal-confidence bars: target 2^-40  -> Miller k=20, Solovay k=40
    def wm(method, k):
        for r in rows:
            if r["method"] == method and int(r["k"]) == k:
                return float(r["word_muls"])
        return 0
    labels = ["Miller-Rabin\nk=20", "Solovay-Strassen\nk=40"]
    vals = [wm("Miller-Rabin", 20), wm("Solovay-Strassen", 40)]
    cols = [METHOD_STYLE["Miller-Rabin"][0], METHOD_STYLE["Solovay-Strassen"][0]]
    bars = ax2.bar(labels, vals, color=cols, width=0.55)
    for b, v in zip(bars, vals):
        ax2.text(b.get_x()+b.get_width()/2, v, f"{v/1e6:.1f}M", ha="center", va="bottom", fontsize=11)
    ax2.set_ylabel("word multiplications")
    ax2.set_title("(b) cost at equal confidence ($2^{-40}$)")
    save(fig, "fig3_confidence")

# ---------------- Figure 4: AKS blowup ----------------
def fig4():
    rows = load("fig4_aks.csv")
    xa, ya = avg_by_bits(rows, "AKS", "word_muls")
    xm, ym = avg_by_bits(rows, "Miller-Rabin", "word_muls")
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(xa, ya, marker="D", color=METHOD_STYLE["AKS"][0], label="AKS", lw=1.8, ms=6)
    ax.plot(xm, ym, marker="^", color=METHOD_STYLE["Miller-Rabin"][0], label="Miller-Rabin", lw=1.8, ms=6)
    # annotate ratio at last common point
    if xa and xm:
        b = xa[-1]; ratio = ya[-1] / ym[-1]
        ax.annotate(f"~{ratio:,.0f}x", xy=(b, ya[-1]), xytext=(b*0.6, ya[-1]),
                    fontsize=11, fontweight="bold", va="center")
    ax.set(yscale="log", xlabel="bit length of n", ylabel="word multiplications (log)")
    ax.legend()
    save(fig, "fig4_aks_blowup")

# ---------------- Figure 5: Lucas-Lehmer ----------------
def fig5():
    rows = load("fig5_lucas.csv")
    # restrict to prime-exponent Mersennes so both tests do full work (fair scaling)
    rows = [r for r in rows if "primeExp" in r["label"]]
    xl, yl = avg_by_bits(rows, "Lucas-Lehmer", "word_muls")
    xm, ym = avg_by_bits(rows, "Miller-Rabin", "word_muls")
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(xl, yl, marker="v", color=METHOD_STYLE["Lucas-Lehmer"][0], label="Lucas-Lehmer (deterministic)", lw=1.8, ms=6)
    ax.plot(xm, ym, marker="^", color=METHOD_STYLE["Miller-Rabin"][0], label="Miller-Rabin (k=10)", lw=1.8, ms=6)
    if xl and xm:
        ratio = ym[-1] / yl[-1]
        ax.annotate(f"~{ratio:,.0f}x cheaper", xy=(xl[-1], yl[-1]), xytext=(xl[-1]*0.35, yl[-1]),
                    fontsize=11, fontweight="bold", va="center")
    ax.set(xscale="log", yscale="log", xlabel="bit length of Mersenne number $M_p$",
           ylabel="word multiplications (log)")
    ax.legend()
    save(fig, "fig5_lucas_vs_miller")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5()
    print("all figures written to", OUT + "/")
