import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def plot_scenario(csv_path, scenario_name):
    df = pd.read_csv(csv_path)
    seq_df = df[df["seq_id"] != "stats"]

    # Key features
    features = ["peak_mag", "snr", "phase_coherence"]

    # Group by label
    groups = seq_df.groupby("label")

    # Generate histograms for each feature
    fig, axes = plt.subplots(1, len(features), figsize=(15, 5))
    if len(features) == 1:
        axes = [axes]

    for i, feat in enumerate(features):
        ax = axes[i]
        for label, group in groups:
            ax.hist(group[feat], alpha=0.7, label=label, bins=20)
        ax.set_title(f"{feat.title()} by Group - {scenario_name}")
        ax.set_xlabel(feat)
        ax.set_ylabel("Frequency")
        ax.legend()

    plt.tight_layout()
    png_path = f"reports/analysis/plots/{scenario_name}_{feat}_hist.png"
    plt.savefig(png_path)
    plt.close()

    # FDR for multiple tests (p-values per feature between groups)
    pvals = []
    for feat in features:
        g1 = groups.get_group("high_gc")[feat].dropna()
        g2 = groups.get_group("low_gc")[feat].dropna()
        _, p = stats.ttest_ind(g1, g2)
        pvals.append(p)

    # Benjamini-Hochberg FDR
    reject, p_adjusted = stats.fdrcorrection(pvals, alpha=0.05)

    # Stats row
    stats_df = pd.DataFrame(
        {
            "feature": features,
            "p_value": pvals,
            "p_adjusted_fdr": p_adjusted,
            "significant": reject,
        }
    )
    stats_df.to_csv(f"results/scenarios/{scenario_name}/fdr_stats.csv", index=False)

    print(f"Plots saved to {png_path.replace(feat, '')[:-4]}*.png")
    print("FDR stats:", stats_df)


# Usage: plot_scenario('path/to/csv', 'scenario1_gc_bias')
