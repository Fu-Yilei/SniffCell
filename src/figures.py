import matplotlib.pyplot as plt
import re
import seaborn as sns
import numpy as np
import pandas as pd


def assign_variant_with_cell_type_names(df, cell_type_names, epsilon=1e-6):
    # Initialize dictionary with cell type names as keys and 0 as default values
    cell_type_dicts = {cell_type: [] for cell_type in cell_type_names}
    for _, row in df.iterrows():
        cell_type_probs = row["cell_type_prob_em"]
        for cell_type, prob in zip(cell_type_names, cell_type_probs):
            cell_type_dicts[cell_type].append(prob)
    return cell_type_dicts

    
def plot_cell_type_box_distributions(cell_type_dicts, output_file=None):
    # build a wide DataFrame and then melt it
    cell_type_df = pd.DataFrame({k: pd.Series(v) for k, v in cell_type_dicts.items()})
    df_long = cell_type_df.melt(var_name="Cell Type", value_name="Probability")

    plt.figure(figsize=(8, 9))
    # boxplot with hue equal to the same variable as x
    sns.boxplot(
        data=df_long,
        x="Cell Type", y="Probability",
        hue="Cell Type",           # map palette onto the categories
        palette="Set2",
        showfliers=False,
        legend=False               # turn off the redundant legend
    )
    # overlay grey dots—either for all cell types or just a specific one
    overlay = df_long


    sns.stripplot(
        data=overlay,
        x="Cell Type", y="Probability",
        color="grey",
        jitter=True,
        size=4,
        alpha=0.6
    )

    plt.xlabel("Cell Type")
    plt.ylabel("Probability")
    plt.title("Cell Type Probability Distributions")
    # plt.xticks(rotation=45, ha="right")
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    if output_file:
        plt.savefig(output_file)