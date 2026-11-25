"""Render timings from the scaling example. Run from any directory."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator

root = Path(__file__).parent
with (root / "scaling.csv").open() as source:
    rows = list(csv.DictReader(source))

plt.rcParams.update({"font.size": 10, "svg.fonttype": "none", "svg.hashsalt": "chazelle"})
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), layout="constrained")
for panel, key, title in zip(axes, ["submap_ns", "sweep_ns"],
                             ["Chazelle submaps", "Comparison sweep"]):
    for family, color, symbol in [("star", "#9c4933", "o"), ("comb", "#76639a", "^"), ("convex", "#236988", "s")]:
        data = [r for r in rows if r["family"] == family]
        vertices = [int(r["vertices"]) for r in data]
        per_vertex = [int(r[key]) / (1000 * n) for r, n in zip(data, vertices)]
        panel.plot(vertices, per_vertex, marker=symbol, markersize=4, color=color,
                   label=family.capitalize(), linewidth=1.8)
    panel.set_xscale("log", base=2)
    panel.set_xticks([128, 1024, 8192, 65536, 1048576])
    panel.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))
    panel.yaxis.set_major_locator(MaxNLocator(nbins=5))
    panel.set_ylim(bottom=0)
    panel.grid(axis="y", alpha=0.25)
    panel.spines[["right", "top"]].set_visible(False)
    panel.set_title(title, loc="left", fontweight="bold")
    panel.set_xlabel("Polygon vertices")
    panel.set_ylabel("Time per vertex (µs)")
axes[0].legend(frameon=False)
fig.savefig(root / "scaling.svg", metadata={"Date": None})
