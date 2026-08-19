import re
import sys

from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ─────────────────────────────────────────────────────────
# PARSING - Handle both formats automatically
# ─────────────────────────────────────────────────────────
def parse_results_file(filepath):
    """
    Parse hash function benchmark results from either format:
    - Old format: hash algorithm | len(D) | len(a) | avg. #steps | avg. time (s)
    - New CLR format: hash algorithm | k | len(D) | len(a) | avg. #steps | avg. time (s)
    
    Returns:
        raw_data: list of dicts with all fields
    """
    raw_data = []
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Split into blocks
    blocks = re.split(r'\n\s*\n', content.strip())
    
    for block in blocks:
        lines = [l.strip() for l in block.split('\n') if l.strip()]
        
        # Detect header to know field positions
        header_line = lines[0]
        fields = header_line.split('|')
        
        # Determine if k column exists
        has_k = 'k' in header_line.lower() or len(fields) == 6
        
        for line in lines[1:]:
            parts = line.split('|')
            if len(parts) >= 5:
                entry = {}
                
                if has_k:
                    entry['algorithm'] = parts[0].strip()
                    entry['k'] = int(parts[1].strip())
                    entry['len_D'] = int(parts[2].strip())
                    entry['len_a'] = int(parts[3].strip())
                    entry['avg_steps'] = int(parts[4].strip())
                    entry['avg_time'] = float(parts[5].strip()) * 10**6
                else:
                    entry['algorithm'] = parts[0].strip()
                    entry['k'] = 0  # Default for single-factor
                    entry['len_D'] = int(parts[1].strip())
                    entry['len_a'] = int(parts[2].strip())
                    entry['avg_steps'] = int(parts[3].strip())
                    entry['avg_time'] = float(parts[4].strip()) * 10**6
                
                raw_data.append(entry)
    
    return raw_data

# Group data structures
def group_raw_data(raw_data):
    """Organize raw data for easy access."""
    
    # By algorithm type: single-factor vs CLR
    single_factor = [d for d in raw_data if not d['algorithm'].startswith('clr-')]
    clr = [d for d in raw_data if d['algorithm'].startswith('clr-')]
    
    # Time storage: (algorithm, k, len_D, len_a) -> list of times
    times_by_config = defaultdict(list)
    steps_by_params = defaultdict(list)  # (k, len_D, len_a) -> list of steps
    
    for d in raw_data:
        key = (d['algorithm'], d['k'], d['len_D'], d['len_a'])
        times_by_config[key].append(d['avg_time'])
        step_key = (d['k'], d['len_D'], d['len_a'])
        steps_by_params[step_key].append(d['avg_steps'])
    
    return {
        'single_factor': single_factor,
        'clr': clr,
        'times_by_config': times_by_config,
        'steps_by_params': steps_by_params
    }

def compute_statistics(data_groups):
    """Compute mean/std/min/max/count for all configurations."""
    
    stats = {}
    for key, values in data_groups['times_by_config'].items():
        stats[key] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'count': len(values)
        }
    
    return stats

# ─────────────────────────────────────────────────────────
# SINGLE-FACTOR VISUALIZATIONS (original functionality preserved)
# ─────────────────────────────────────────────────────────
def plot_single_factor_grouped_bars(stats, output_path="single_factor_bars.png"):
    """Time vs len(a), grouped by algorithm, one panel per len(D)."""
    
    # Filter single-factor only (k=0)
    sf_stats = {k: v for k, v in stats.items() if k[1] == 0}
    if not sf_stats:
        print("No single-factor data found")
        return None
    
    algorithms = sorted(set(k[0] for k in sf_stats.keys()))
    len_D_vals = sorted(set(k[2] for k in sf_stats.keys()))
    len_a_vals = sorted(set(k[3] for k in sf_stats.keys()))
    
    n_algos = len(algorithms)
    bar_w = 0.13
    x = np.arange(len(len_a_vals))
    
    fig, axes = plt.subplots(1, len(len_D_vals), sharey=False)
    if len(len_D_vals) == 1:
        axes = [axes]
    
    ALGO_COLORS = {
        "pr": "#e74c3c", "pr-dma": "#e67e22",
        "pp": "#3498db", "pp-dma": "#2ecc71",
        "pp-mr2": "#9b59b6", "pp-mr2-dma": "#1abc9c"
    }
    
    for ax_idx, len_d in enumerate(len_D_vals):
        ax = axes[ax_idx]
        for i, algo in enumerate(algorithms):
            offsets = (i - n_algos / 2) * bar_w + bar_w / 2
            means = np.array([sf_stats[(algo, 0, len_d, la)]['mean'] 
                            for la in len_a_vals])
            stds = np.array([sf_stats[(algo, 0, len_d, la)]['std'] 
                           for la in len_a_vals])
            ax.bar(x + offsets, means, width=bar_w, yerr=stds, capsize=2,
                   label=algo, color=ALGO_COLORS.get(algo, None))
        ax.set_title(f"len(D) = {len_d}")
        ax.set_xlabel("len(a)")
        ax.set_xticks(x)
        ax.set_xticklabels(len_a_vals)
        if ax_idx == 0:
            ax.set_ylabel("Avg. time (μs)")
    
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center",
               ncol=n_algos, bbox_to_anchor=(0.5, 1.06), frameon=False)
    fig.suptitle("Single-Factor Hash Time by Algorithm and len(a)", y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    return fig

# ─────────────────────────────────────────────────────────
# CLR-SPECIFIC VISUALIZATIONS
# ─────────────────────────────────────────────────────────
def plot_clr_time_vs_k(stats, output_path="clr_time_vs_k.png"):
    """
    CLR time vs k, with separate lines for each algorithm variant,
    grouped by len(D). Shows which k is optimal and compares implementations.
    """
    
    # Filter CLR only (k > 0)
    clr_stats = {k[:3]: v for k, v in stats.items() if k[1] > 0 and k[0].startswith('clr-')}
    if not clr_stats:
        print("No CLR data found")
        return None
    
    algorithms = sorted(set(k[0] for k in clr_stats.keys()))
    len_D_vals = sorted(set(k[2] for k in clr_stats.keys()))
    k_vals = sorted(set(k[1] for k in clr_stats.keys()))
    
    CLR_COLORS = {"clr-pr": "#e74c3c", "clr-pp-dma": "#2ecc71", "clr-pp-mr2-dma": "#1abc9c"}
    
    fig, axes = plt.subplots(1, len(len_D_vals), sharey=True, figsize=(12,5), layout='constrained')
    if len(len_D_vals) == 1:
        axes = [axes]
    
    for ax_idx, len_d in enumerate(len_D_vals):
        ax = axes[ax_idx]
        for algo in algorithms:
            times = np.array([clr_stats[(algo, k, len_d)]['mean'] 
                            if isinstance(clr_stats.get((algo, k, len_d)), dict)
                            else clr_stats[(algo, k, len_d)]['mean']
                            for k in k_vals])
            stds = np.array([clr_stats[(algo, k, len_d)]['std'] for k in k_vals])
            
            marker = 'o' if algo == "clr-pr" else 's'
            ax.plot(k_vals, times, marker=marker, linewidth=2, markersize=8,
                   label=algo.replace('clr-', ''), color=CLR_COLORS.get(algo, None))
            
            # Highlight min point
            min_k = k_vals[np.argmin(times)]
            min_val = np.min(times)
            ax.scatter(min_k, min_val, s=150, zorder=5, marker='*',
                      color='gold', edgecolor='black', linewidth=1.5)
        
        ax.set_title(f"len(D) = {len_d}")
        ax.set_xlabel("k (number of small factors)")
        ax.set_ylabel("Avg. time (μs)")
        ax.set_xticks(k_vals)
        ax.grid(True, alpha=0.3, linestyle='--')
    
    # Legend using first axis
    handles = []
    labels = []
    for algo in algorithms:
        h = plt.Line2D([], [], marker='o', color=CLR_COLORS.get(algo, None),
                      linewidth=2, markersize=8, label=algo.replace('clr-', ''))
        handles.append(h)
        labels.append(algo.replace('clr-', ''))
    
    fig.legend(handles, labels, loc="lower left", ncol=len(algorithms))
    #,
    #          bbox_to_anchor=(0.5, 1.05), frameon=True, fancybox=True)
    
    fig.suptitle("CLR Hash Time vs k (star marks optimum per curve)", y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    return fig

def plot_clr_speedup_heatmap(stats, output_path="clr_speedup_heatmap.png"):
    """
    Heatmap showing speedup (%) of clr-pp-dma and clr-pp-mr2-dma
    vs clr-pr baseline across k and len(D).
    """
    
    clr_stats = {k[:3]: v for k, v in stats.items() if k[1] > 0 and k[0].startswith('clr-')}
    if not clr_stats:
        print("No CLR data found")
        return None
    
    len_D_vals = sorted(set(k[2] for k in clr_stats.keys()))
    k_vals = sorted(set(k[1] for k in clr_stats.keys()))

    pr_opt_k = {}

    for len_d in len_D_vals:
        pr_times = [(k, clr_stats.get(('clr-pr', k, len_d), {'mean': None})['mean']) 
                   for k in set(k[1] for k in clr_stats.keys() if k[2] == len_d)]
        pr_opt_k[len_d], _ = min(pr_times, key=lambda x: x[1] or float('inf'))

    # Extract baselines (clr-pr-3) and optimized versions
    baseline = np.zeros((len(len_D_vals), len(k_vals)))
    pp_dma = np.zeros((len(len_D_vals), len(k_vals)))
    mr2_dma = np.zeros((len(len_D_vals), len(k_vals)))
    
    for di, len_d in enumerate(len_D_vals):
        for ki, k in enumerate(k_vals):
            baseline[di, ki] = clr_stats[('clr-pr', pr_opt_k[len_d], len_d)]['mean']
            pp_dma[di, ki] = clr_stats[('clr-pp-dma', k, len_d)]['mean']
            mr2_dma[di, ki] = clr_stats[('clr-pp-mr2-dma', k, len_d)]['mean']
    
    # Compute speedup (%)
    pp_speedup = (baseline - pp_dma) / baseline * 100
    mr2_speedup = (baseline - mr2_dma) / baseline * 100
    
    fig, axes = plt.subplots(1, 2, figsize=(15,5))
    
    titles = ["clr-pp-dma vs clr-pr", "clr-pp-mr2-dma vs clr-pr"]
    arrays = [pp_speedup, mr2_speedup]
    
    for ax_idx, (ax, arr, title) in enumerate(zip(axes, arrays, titles)):
        limit = max(np.max(arr), -np.min(arr))*1.1
        im = ax.imshow(arr, cmap="RdYlGn", aspect="auto",
                       vmin=-limit, vmax=limit)
        
        # Y-axis: len(D)
        ax.set_yticks(np.arange(len(len_D_vals)))
        ax.set_yticklabels([f"{d}" for d in len_D_vals])
        
        # X-axis: k
        ax.set_xticks(np.arange(len(k_vals)))
        ax.set_xticklabels(k_vals)
        ax.set_xlabel("k")
        
        # Annotate cells
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                val = arr[i, j]
                ax.text(j, i, f"{val:.1f}%", ha="center", va="center",
                       fontsize=9, color="black")
        
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.8)
    
    fig.suptitle("CLR Optimization Speedup vs Baseline (%)", y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    return fig

# ─────────────────────────────────────────────────────────
# SUMMARY REPORT
# ─────────────────────────────────────────────────────────
def print_summary_report(stats, data_groups):
    """Print a human-readable summary of findings."""
    
    print("=" * 70)
    print("HASH FUNCTION BENCHMARK SUMMARY REPORT")
    print("=" * 70)
    
    # Single-factor section
    sf_keys = [k for k in stats if k[1] == 0]
    if sf_keys:
        print("\n--- SINGLE-FACTOR HASHES (k=0) ---")
        print(f"{'Algorithm':<18} {'len(D)':>8} {'len(a)':>8} {'Mean Time (μs)':>16} {'Std Dev':>9}")
        print("-" * 70)
        for key in sorted(sf_keys, key=lambda x: (x[2], x[3], x[0]))[:24]:  # Limit output
            algo, k, len_d, len_a = key
            s = stats[key]
            print(f"{algo:<18} {len_d:>8} {len_a:>8} {s['mean']:16.3f} {s['std']:9.3f}")
    
    # CLR section
    clr_stats = {k[:3]: v for k, v in stats.items() if k[1] > 0 and k[0].startswith('clr-')}
    clr_keys = clr_stats.keys()
    if clr_stats:
        print("\n--- CLR HASHES (k>0) ---")

        len_D_vals = sorted(set(k[2] for k in clr_keys))
        
        # Best performing combinations
        best = [min((k for k in clr_keys if k[2] == len_d), key=lambda k: clr_stats[k]['mean']) for len_d in len_D_vals]
        best = [[k for k in stats.keys() if k[:3] == target][0] for target in best]

        print(f"\nOptimal configurations:\n")

        for k in best:
            print(f"  {k[0]} at k={k[1]}, len(D)={k[2]}")
            print(f"  Time: {int(stats[k]['mean'])} μs, Steps: {data_groups['steps_by_params'].get((k[1], k[2], k[3]), [])[0] if data_groups['steps_by_params'].get((k[1], k[2], k[3])) else 'N/A'}")
        
        # Optimal k per len(D)
        print("\nOptimal k by len(D):")
        print(f"{'len(D)':>8} {'clr-pr':<8} {'clr-pp-dma':<12} {'clr-pp-mr2-dma':<16} {'Speedup':<9}")
        print("-" * 55)
            
        for len_d in len_D_vals:
            keys_at_d = [(k, clr_stats[k]) for k in clr_keys if k[2] == len_d]
            
            pr_times = [(k, clr_stats.get(('clr-pr', k, len_d), {'mean': None})['mean']) 
                       for k in set(k[1] for k in clr_keys if k[2] == len_d)]
            pr_opt_k, pr_opt_t = min(pr_times, key=lambda x: x[1] or float('inf'))

            pp_times = [(k, clr_stats.get(('clr-pp-dma', k, len_d), {'mean': None})['mean']) 
                       for k in set(k[1] for k in clr_keys if k[2] == len_d)]
            pp_opt_k, pp_opt_t = min(pp_times, key=lambda x: x[1] or float('inf'))
            
            mr2_times = [(k, clr_stats.get(('clr-pp-mr2-dma', k, len_d), {'mean': None})['mean'])
                        for k in set(k[1] for k in clr_keys if k[2] == len_d)]
            mr2_opt_k, mr2_opt_t = min(mr2_times, key=lambda x: x[1] or float('inf'))

            baseline_t = clr_stats[('clr-pr', pr_opt_k, len_d)]['mean']
            speedup = (baseline_t - pp_opt_t) / baseline_t * 100
            
            print(f"{len_d:>8} k={pr_opt_k:<6} k={pp_opt_k:<10} k={mr2_opt_k:<14} {f'{speedup:<.1f}%':<9}")
    
    print("\n" + "=" * 70)

# ─────────────────────────────────────────────────────────
# MAIN DRIVER
# ─────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Parse file (handles both old and new formats automatically)
    filepath = sys.argv[1]
    
    print(f"Parsing {filepath}...")
    raw_data = parse_results_file(filepath)
    print(f"Found {len(raw_data)} entries")
    
    # Organize and compute statistics
    data_groups = group_raw_data(raw_data)
    stats = compute_statistics(data_groups)
    
    # Print summary
    print_summary_report(stats, data_groups)
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    
    # Single-factor (if present)
    if data_groups['single_factor']:
        plot_single_factor_grouped_bars(stats, "single_factor_bars.png")
        print("  ✓ single_factor_bars.png")
    
    # CLR-specific
    if data_groups['clr']:
        plot_clr_time_vs_k(stats, "clr_time_vs_k.png")
        print("  ✓ clr_time_vs_k.png")
        
        plot_clr_speedup_heatmap(stats, "clr_speedup_heatmap.png")
        print("  ✓ clr_speedup_heatmap.png")
        
    plt.show()
    print("\nDone!")
