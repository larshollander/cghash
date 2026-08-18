import re
import numpy as np
import matplotlib.pyplot as plt

from collections import defaultdict

def parse_results_file(filepath):
    """
    Parse hash function benchmark results from the specified file.
    
    Returns:
        times_by_algo: dict[(algorithm, len_D, len_a)] -> list of time values
        steps_by_params: dict[(len_D, len_a)] -> list of step values
    """
    
    times_by_algo = defaultdict(list)
    steps_by_params = defaultdict(list)
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Split into blocks (each block separated by empty lines)
    blocks = re.split(r'\n\s*\n', content.strip())
    
    for block in blocks:
        lines = [l.strip() for l in block.split('\n') if l.strip()]
        
        # Skip header line
        for line in lines[1:]:
            parts = line.split('|')
            if len(parts) >= 5:
                algorithm = parts[0].strip()
                len_D = int(parts[1].strip())
                len_a = int(parts[2].strip())
                avg_steps = int(parts[3].strip())
                avg_time_str = parts[4].strip()
                
                # Convert scientific notation to float
                avg_time = float(avg_time_str) * 10**6
                
                # Store raw values for averaging across runs
                times_by_algo[(algorithm, len_D, len_a)].append(avg_time)
                steps_by_params[(len_D, len_a)].append(avg_steps)
    
    return times_by_algo, steps_by_params


def compute_statistics(times_by_algo, steps_by_params):
    """
    Compute mean and variability statistics for all parameter combinations.
    
    Returns:
        times_stats: dict[(algorithm, len_D, len_a)] -> {'mean': float, 'std': float, 'count': int}
        steps_stats: dict[(len_D, len_a)] -> {'mean': int, 'std': float, 'count': int}
    """
    
    times_stats = {}
    for key, values in times_by_algo.items():
        times_stats[key] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'count': len(values)
        }
    
    steps_stats = {}
    for key, values in steps_by_params.items():
        steps_stats[key] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'count': len(values)
        }
    
    return times_stats, steps_stats


def print_summary(times_stats, steps_stats):
    """Print a formatted summary of the results."""
    
    print("=" * 70)
    print("HASH FUNCTION BENCHMARK RESULTS SUMMARY")
    print("=" * 70)
    
    print("\n--- TIME STATISTICS BY ALGORITHM AND PARAMETERS ---")
    print(f"{'Algorithm':<15} {'len(D)':>8} {'len(a)':>8} {'Mean Time(μs)':>14} {'Std Dev(μs)':>14} {'Runs':>6}")
    print("-" * 70)
    
    sorted_times = sorted(times_stats.keys(), key=lambda k: (k[1], k[2], k[0]))
    for algo, len_d, len_a in sorted_times:
        stats = times_stats[(algo, len_d, len_a)]
        print(f"{algo:<15} {len_d:>8} {len_a:>8} {stats['mean']:.6e} {stats['std']:.6e} {stats['count']:>6}")
    
    print("\n--- STEP COUNT STATISTICS BY PARAMETERS ---")
    print(f"{'len(D)':>8} {'len(a)':>8} {'Mean Steps':>12} {'Std Dev':>12} {'Runs':>6}")
    print("-" * 70)
    
    sorted_steps = sorted(steps_stats.keys(), key=lambda k: (k[0], k[1]))
    for len_d, len_a in sorted_steps:
        stats = steps_stats[(len_d, len_a)]
        print(f"{len_d:>8} {len_a:>8} {stats['mean']:>10.1f} {stats['std']:>10.1f} {stats['count']:>6}")


# Alternative: Return as numpy arrays for further processing
def convert_to_numpy_arrays(times_stats, steps_stats):
    """
    Convert stats dictionaries to numpy-friendly structured arrays.
    
    Useful for plotting with matplotlib/seaborn or saving to CSV/Parquet.
    """
    
    # Time array: [algo_idx, len_D, len_a] -> mean_time
    algorithms = sorted(set(k[0] for k in times_stats.keys()))
    len_D_values = sorted(set(k[1] for k in times_stats.keys()))
    len_a_values = sorted(set(k[2] for k in times_stats.keys()))
    
    algo_to_idx = {a: i for i, a in enumerate(algorithms)}
    len_D_to_idx = {d: i for i, d in enumerate(len_D_values)}
    len_a_to_idx = {a: i for i, a in enumerate(len_a_values)}
    
    time_means = np.zeros((len(algorithms), len(len_D_values), len(len_a_values)))
    time_stds = np.zeros((len(algorithms), len(len_D_values), len(len_a_values)))
    
    for key, stats in times_stats.items():
        algo, len_d, len_a = key
        time_means[algo_to_idx[algo], len_D_to_idx[len_d], len_a_to_idx[len_a]] = stats['mean']
        time_stds[algo_to_idx[algo], len_D_to_idx[len_d], len_a_to_idx[len_a]] = stats['std']
    
    # Step count array: [len_D, len_a] -> mean_steps
    steps_means = np.zeros((len(len_D_values), len(len_a_values)))
    steps_stds = np.zeros((len(len_D_values), len(len_a_values)))
    
    for key, stats in steps_stats.items():
        len_d, len_a = key
        steps_means[len_D_to_idx[len_d], len_a_to_idx[len_a]] = stats['mean']
        steps_stds[len_D_to_idx[len_d], len_a_to_idx[len_a]] = stats['std']
    
    return {
        'algorithms': algorithms,
        'len_D': len_D_values,
        'len_a': len_a_values,
        'time_means': time_means,
        'time_stds': time_stds,
        'step_means': steps_means,
        'step_stds': steps_stds
    }

ALGO_COLORS = {
    "pr":          "#e74c3c",
    "pr-dma":      "#e67e22",
    "pp":          "#3498db",
    "pp-dma":      "#2ecc71",
    "pp-mr2":      "#9b59b6",
    "pp-mr2-dma":  "#1abc9c",
}

# ─────────────────────────────────────────────────────────
# 1. Grouped bar chart: avg. time vs len(a), one panel per len(D)
# ─────────────────────────────────────────────────────────
def plot_time_grouped_bars(numpy_data):
    algorithms = numpy_data["algorithms"]
    len_D_vals = numpy_data["len_D"]
    len_a_vals = numpy_data["len_a"]
    means      = numpy_data["time_means"]   # (n_algos, n_D, n_a)
    stds       = numpy_data["time_stds"]

    n_algos = len(algorithms)
    bar_w = 0.13
    x = np.arange(len(len_a_vals))

    fig, axes = plt.subplots(1, len(len_D_vals), figsize=(5 * len(len_D_vals), 5),
                             sharey=False)

    for ax_idx, len_d in enumerate(len_D_vals):
        ax = axes[ax_idx]
        for i, algo in enumerate(algorithms):
            offsets = (i - n_algos / 2) * bar_w + bar_w / 2
            ax.bar(x + offsets, means[i, ax_idx],
                   width=bar_w, yerr=stds[i, ax_idx],
                   capsize=2, label=algo,
                   color=ALGO_COLORS.get(algo, None))
        ax.set_title(f"len(D) = {len_d}")
        ax.set_xlabel("len(a)")
        ax.set_xticks(x)
        ax.set_xticklabels(len_a_vals)
        if ax_idx == 0:
            ax.set_ylabel("Avg. time (μs)")

    # Single legend at top
    fig.legend(*axes[0].get_legend_handles_labels(),
               loc="upper center", ncol=n_algos,
               bbox_to_anchor=(0.5, 1.06), frameon=False)

    fig.suptitle("Hash function time by algorithm and len(a)", y=1.02)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────
# 2. Line chart: avg. step count vs len(a), one line per len(D)
# ─────────────────────────────────────────────────────────
def plot_steps_line(numpy_data):
    len_D_vals = numpy_data["len_D"]
    len_a_vals = numpy_data["len_a"]
    step_means = numpy_data["step_means"]  # (n_D, n_a)
    step_stds  = numpy_data["step_stds"]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    for j, len_d in enumerate(len_D_vals):
        ax.errorbar(len_a_vals, step_means[j], yerr=step_stds[j],
                    marker="o", capsize=4, linewidth=2,
                    label=f"len(D) = {len_d}")

    ax.set_xlabel("len(a)")
    ax.set_ylabel("Avg. number of steps")
    ax.set_title("Step count vs len(a)")
    ax.legend(frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────
# 3. Heatmap: speedup (%) of each algorithm vs baseline
#    Rows = algorithms, columns = (len(D), len(a)) combos
# ─────────────────────────────────────────────────────────
def plot_speedup_heatmap(numpy_data, baseline_algo):
    algorithms = numpy_data["algorithms"]
    len_D_vals = numpy_data["len_D"]
    len_a_vals = numpy_data["len_a"]
    means      = numpy_data["time_means"]

    bl_idx = algorithms.index(baseline_algo)
    baseline = means[bl_idx]                      # (n_D, n_a)
    speedup_pct = (baseline[None, :, :] - means) / baseline[None, :, :] * 100

    # Flatten (len_D, len_a) into column labels
    col_labels = [f"D={d}, a={a}" for d in len_D_vals for a in len_a_vals]
    row_labels = [a for a in algorithms if a != baseline_algo]
    non_bl_mask = np.array([a != baseline_algo for a in algorithms])
    data = speedup_pct[non_bl_mask].reshape(non_bl_mask.sum(), -1)

    fig, ax = plt.subplots(figsize=(10, 4))
    limit = max(data.max(), -data.min())
    im = ax.imshow(data, cmap="RdYlGn", aspect="auto",
                   vmin=-limit, vmax=limit)

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)

    # Annotate cells
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{data[i,j]:.1f}%",
                    ha="center", va="center", fontsize=7,
                    color="white" if abs(data[i, j]) > data.max()*0.6 else "black")

    ax.set_title("Speedup vs baseline (%)")
    fig.colorbar(im, ax=ax, label="Speedup (%)", shrink=0.8)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────
# 4. Box plot: variability across runs (needs multiple runs)
#    Pass the raw times_by_algo dict, pick a specific len(D)
# ─────────────────────────────────────────────────────────
def plot_variability_box(times_by_algo, len_d):
    """
    Box plot of raw time values per algorithm at a fixed len(D).
    Only meaningful when you have multiple runs per config.
    """
    # Find all len_a values for this len_d
    len_a_set = sorted({k[2] for k in times_by_algo if k[1] == len_d})
    algorithms = sorted({k[0] for k in times_by_algo})

    fig, ax = plt.subplots(figsize=(9, 5))

    positions = []
    box_data = []
    labels = []
    pos = 0

    for len_a in len_a_set:
        for algo in algorithms:
            vals = times_by_algo.get((algo, len_d, len_a), [])
            if len(vals) >= 2:  # need ≥2 runs for a box
                box_data.append(vals)
                positions.append(pos)
                labels.append(f"{algo}\na={len_a}")
                pos += 1
        pos += 0.5  # gap between len(a) groups

    if not box_data:
        ax.text(0.5, 0.5, "Not enough runs for box plots (need ≥2 per config)",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
    else:
        bp = ax.boxplot(box_data, positions=positions, widths=0.7,
                        patch_artist=True, showmeans=True,
                        meanprops=dict(marker="D", markerfacecolor="white",
                                       markeredgecolor="black", markersize=5))
        for patch, algo in zip(
                [p for p in bp["boxes"]],
                [lbl.split("\n")[0] for lbl in labels]):
            patch.set_facecolor(ALGO_COLORS.get(algo, "grey"))
            patch.set_alpha(0.7)

        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_ylabel("Avg. time (μs)")
        ax.set_title(f"Run-to-run variability at len(D) = {len_d}")

    fig.tight_layout()
    return fig


if __name__ == "__main__":
    filename = "results.txt"

    # Parse raw data
    times_by_algo, steps_by_params = parse_results_file(filename)
    
    # Compute statistics
    times_stats, steps_stats = compute_statistics(times_by_algo, steps_by_params)
    
    # Optionally convert to numpy arrays for further analysis
    numpy_data = convert_to_numpy_arrays(times_stats, steps_stats)

    fig1 = plot_time_grouped_bars(numpy_data)
    fig1.savefig("time_grouped_bars.png", dpi=150, bbox_inches="tight")

    fig2 = plot_steps_line(numpy_data)
    fig2.savefig("steps_vs_len_a.png", dpi=150, bbox_inches="tight")

    fig3 = plot_speedup_heatmap(numpy_data, "pr")
    fig3.savefig("speedup_heatmap.png", dpi=150, bbox_inches="tight")

    fig4 = plot_variability_box(times_by_algo, len_d=2048)
    fig4.savefig("variability_box_D2048.png", dpi=150, bbox_inches="tight")

    plt.show()
