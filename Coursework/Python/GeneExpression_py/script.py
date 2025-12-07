# Load packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

%matplotlib inline

# My two DESeq2 result files, I've used shortened comparison names
file_ij = "IJ.tsv"
file_ie = "IE.tsv"

# Read the TSV files
res_ij = pd.read_csv(file_ij, sep="\t")
res_ie = pd.read_csv(file_ie, sep="\t")

# look at
res_ij.head()

# Check basic info
res_ij.info()
res_ie.info()

# Check missing values
print("Missing values in IJ file:")
print(res_ij.isna().sum())

print("Missing values in IE file:")
print(res_ie.isna().sum())

# Basic statistics
res_ij.describe()

# count significant genes (up / down)

# my threshold
p_cut = 0.05
fc_cut = 1

# make short name
fc_ij = res_ij["log2FoldChange"]
padj_ij = res_ij["padj"]

# IJ: up and down
up_ij = res_ij[(padj_ij < p_cut) & (fc_ij > fc_cut)]
down_ij = res_ij[(padj_ij < p_cut) & (fc_ij < -fc_cut)]

print("IJ result:")
print("up genes:", len(up_ij))
print("down genes:", len(down_ij))

# same thing for IE

p_cut = 0.05
fc_cut = 1

# make short name for IE
fc_ie = res_ie["log2FoldChange"]
padj_ie = res_ie["padj"]

up_ie = res_ie[(padj_ie < p_cut) & (fc_ie > fc_cut)]
down_ie = res_ie[(padj_ie < p_cut) & (fc_ie < -fc_cut)]

# volcano plot for IJ 
df = res_ij.copy()

# make a new column for the group
df["sig"] = "not_sig"

# make short names
fc = df["log2FoldChange"]
p = df["padj"]

# up genes
df["sig"][(p < 0.05) & (fc > 1)] = "up"

# down genes
df["sig"][(p < 0.05) & (fc < -1)] = "down"

# volcano for IJ
plt.figure(figsize=(6,5))

sns.scatterplot(
    data=df,
    x="log2FoldChange",
    y=-np.log10(df["padj"]),
    hue="sig",
    palette={"not_sig":"grey", "up":"red", "down":"green"},
    alpha=0.6
)

plt.title("Volcano plot (IJ)")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(padj)")

plt.tight_layout()
plt.show()

# Summary statistics for IJ
print("Summary for IJ")

print("Significant up genes:", up_ij)
print("Significant down genes:", down_ij)

print("log2FoldChange summary:")
print(res_ij["log2FoldChange"].describe())

print("padj summary:")
print(res_ij["padj"].describe())


# Summary statistics for IE
print("Summary for I")

print("Significant up genes:", up_ie)
print("Significant down genes:", down_ie)

print("log2FoldChange summary:")
print(res_ie["log2FoldChange"].describe())

print("padj summary:")
print(res_ie["padj"].describe())

import matplotlib.pyplot as plt
import numpy as np

# make a volcano plot for IJ comparison
df = res_ij.copy()

# make a new column
df["sig"] = "no"
df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] > 1), "sig"] = "up"
df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] < -1), "sig"] = "down"

# set colors
colors = {"no": "grey", "up": "red", "down": "green"}

# volcano plot
plt.figure(figsize=(6,5))

# volcano plot for IE
df_ie = res_ie.copy()

df_ie["sig"] = "no"
fc_ie = df_ie["log2FoldChange"]
p_ie = df_ie["padj"]

df_ie.loc[(p_ie < 0.05) & (fc_ie > 1), "sig"] = "up"
df_ie.loc[(p_ie < 0.05) & (fc_ie < -1), "sig"] = "down"

colors_ie = {"no": "grey", "up": "red", "down": "blue"}

plt.figure(figsize=(6, 5))
plt.scatter(
    df_ie["log2FoldChange"],
    -np.log10(df_ie["padj"]),
    c=df_ie["sig"].map(colors_ie),
    s=10,
    alpha=0.6,
)

plt.xlabel("log2 fold change")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano plot (IE)")
plt.tight_layout()

plt.savefig("volcano_IE.png", dpi=300)
plt.show()

# MA plot for IJ
plt.figure(figsize=(6, 5))

x_ij = res_ij["baseMean"]
y_ij = res_ij["log2FoldChange"]

plt.scatter(x_ij, y_ij, s=5, alpha=0.5)
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)

plt.xlabel("baseMean")
plt.ylabel("log2 fold change")
plt.title("MA plot (IJ)")
plt.tight_layout()

plt.savefig("MA_IJ.png",)
plt.show()

# MA plot for IE
plt.figure(figsize=(6, 5))

x_ie = res_ie["baseMean"]
y_ie = res_ie["log2FoldChange"]

plt.scatter(x_ie, y_ie, s=5, alpha=0.5)
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)

plt.xlabel("baseMean")
plt.ylabel("log2 fold change")
plt.title("MA plot (IE)")
plt.tight_layout()
plt.savefig("MA_IE.png", dpi=300)
plt.show()

# volcano plot for IE
df_ie = res_ie.copy()

df_ie["sig"] = "no"
fc_ie = df_ie["log2FoldChange"]
p_ie = df_ie["padj"]

df_ie.loc[(p_ie < 0.05) & (fc_ie > 1), "sig"] = "up"
df_ie.loc[(p_ie < 0.05) & (fc_ie < -1), "sig"] = "down"

colors_ie = {"no": "grey", "up": "red", "down": "blue"}

plt.figure(figsize=(6, 5))
plt.scatter(
    df_ie["log2FoldChange"],
    -np.log10(df_ie["padj"]),
    c=df_ie["sig"].map(colors_ie),
    s=10,
    alpha=0.6,
)

plt.xlabel("log2 fold change")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano plot (IE)")
plt.tight_layout()

plt.savefig("volcano_IE.png", dpi=300)
plt.show()


up_ie = res_ie[(p_ie < p_cut) & (fc_ie > fc_cut)]
down_ie = res_ie[(p_ie < p_cut) & (fc_ie < -fc_cut)]

print("IE:")
print("  up =", len(up_ie))
print("  down =", len(down_ie))

# summary table for both comparisons

summary = pd.DataFrame({
    "comparison": ["IJ", "IE"],
    "up_genes": [len(up_ij), len(up_ie)],
    "down_genes": [len(down_ij), len(down_ie)]
})

summary

# volcano plot for IJ 

df = res_ij.copy()

# make a new column for the group
df["sig"] = "not_sig"

# make short names
fc = df["log2FoldChange"]
p = df["padj"]

# up genes
df["sig"][(p < 0.05) & (fc > 1)] = "up"

# down genes
df["sig"][(p < 0.05) & (fc < -1)] = "down"

# volcano for IJ
plt.figure(figsize=(6,5))

sns.scatterplot(
    data=df,
    x="log2FoldChange",
    y=-np.log10(df["padj"]),
    hue="sig",
    palette={"not_sig":"grey", "up":"red", "down":"green"},
    alpha=0.6
)

plt.title("Volcano plot (IJ)")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(padj)")

plt.tight_layout()
plt.show()

# Summary statistics for IJ
print("Summary for IJ")

print("Significant up genes:", up_ij)
print("Significant down genes:", down_ij)

print("log2FoldChange summary:")
print(res_ij["log2FoldChange"].describe())

print("padj summary:")
print(res_ij["padj"].describe())


# Summary statistics for IE
print("Summary for I")

print("Significant up genes:", up_ie)
print("Significant down genes:", down_ie)

print("log2FoldChange summary:")
print(res_ie["log2FoldChange"].describe())

print("padj summary:")
print(res_ie["padj"].describe())

import matplotlib.pyplot as plt
import numpy as np

# make a volcano plot for IJ comparison
df = res_ij.copy()

# make a new column
df["sig"] = "no"
df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] > 1), "sig"] = "up"
df.loc[(df["padj"] < 0.05) & (df["log2FoldChange"] < -1), "sig"] = "down"

# set colors
colors = {"no": "grey", "up": "red", "down": "green"}

# volcano plot
plt.figure(figsize=(6,5))
plt.scatter(
    df["log2FoldChange"],
    -np.log10(df["padj"]),
    c=df["sig"].map(colors),
    alpha=0.6,
    s=10
)

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano Plot for IJ (simple version)")

plt.tight_layout()
plt.savefig("volcano_IJ.png",)
plt.show()

# volcano plot for IE
df_ie = res_ie.copy()

df_ie["sig"] = "no"
fc_ie = df_ie["log2FoldChange"]
p_ie = df_ie["padj"]

df_ie.loc[(p_ie < 0.05) & (fc_ie > 1), "sig"] = "up"
df_ie.loc[(p_ie < 0.05) & (fc_ie < -1), "sig"] = "down"

colors_ie = {"no": "grey", "up": "red", "down": "blue"}

plt.figure(figsize=(6, 5))
plt.scatter(
    df_ie["log2FoldChange"],
    -np.log10(df_ie["padj"]),
    c=df_ie["sig"].map(colors_ie),
    s=10,
    alpha=0.6,
)

plt.xlabel("log2 fold change")
plt.ylabel("-log10 adjusted p-value")
plt.title("Volcano plot (IE)")
plt.tight_layout()

plt.savefig("volcano_IE.png", dpi=300)
plt.show()

# MA plot for IJ
plt.figure(figsize=(6, 5))

x_ij = res_ij["baseMean"]
y_ij = res_ij["log2FoldChange"]

plt.scatter(x_ij, y_ij, s=5, alpha=0.5)
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)

plt.xlabel("baseMean")
plt.ylabel("log2 fold change")
plt.title("MA plot (IJ)# MA plot for IE
plt.figure(figsize=(6, 5))

x_ie = res_ie["baseMean"]
y_ie = res_ie["log2FoldChange"]

plt.scatter(x_ie, y_ie, s=5, alpha=0.5)
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)

plt.xlabel("baseMean")
plt.ylabel("log2 fold change")
plt.title("MA plot (IE)")
plt.tight_layout()

plt.savefig("MA_IE.png", dpi=300)
plt.show()
")
plt.tight_layout()

plt.savefig("MA_IJ.png",)
plt.show()

# MA plot for IE
plt.figure(figsize=(6, 5))

x_ie = res_ie["baseMean"]
y_ie = res_ie["log2FoldChange"]

plt.scatter(x_ie, y_ie, s=5, alpha=0.5)
plt.axhline(0, color="black", linestyle="--", linewidth=0.5)

plt.xlabel("baseMean")
plt.ylabel("log2 fold change")
plt.title("MA plot (IE)")
plt.tight_layout()

plt.savefig("MA_IE.png", dpi=300)
plt.show()

# p-value histogram for IJ
plt.figure(figsize=(6, 4))

plt.hist(res_ij["pvalue"].dropna(), bins=40, edgecolor="black")
plt.xlabel("p-value")
plt.ylabel("count")
plt.title("P-value histogram (IJ)")
plt.tight_layout()

plt.savefig("p_hist_IJ.png",)
plt.show()

# p-value histogram for IE
plt.figure(figsize=(6, 4))

plt.hist(res_ie["pvalue"].dropna(), bins=40, edgecolor="black")
plt.xlabel("p-value")
plt.ylabel("count")
plt.title("P-value histogram (IE)")
plt.tight_layout()

plt.savefig("p_hist_IE.png",)
plt.show()

# heatmap for top DE genes (IJ and IE)

top_n = 50

# take top genes 
top_ij = res_ij.sort_values("padj").head(top_n)
top_ie = res_ie.sort_values("padj").head(top_n)

# find common genes
common_genes = set(top_ij["gene_id"]) & set(top_ie["gene_id"])
common_genes = list(common_genes)

print("number of common top genes:", len(common_genes))

if len(common_genes) > 0:
    # make small matrix: rows = genes, columns = IJ/IE log2FC
    ij_fc = res_ij.set_index("gene_id").loc[common_genes, "log2FoldChange"]
    ie_fc = res_ie.set_index("gene_id").loc[common_genes, "log2FoldChange"]

    heat_df = pd.DataFrame({
        "IJ_log2FC": ij_fc,
        "IE_log2FC": ie_fc
    })

    heat_df = heat_df.sort_index()

    plt.figure(figsize=(5, 8))
    sns.heatmap(
        heat_df,
        cmap="RdBu_r",
        center=0,
        cbar_kws={"label": "log2 fold change"}
    )
    plt.title("Top DE genes (IJ vs IE)")
    plt.ylabel("gene_id")
    plt.xlabel("comparison")
    plt.tight_layout()

    plt.savefig("heatmap_top_genes.png",)
    plt.show()
else:
    print("no common genes in top lists, cannot make heatmap")

# save significant gene lists for IJ and IE

up_ij.to_csv("sig_up_IJ.csv", index=False)
down_ij.to_csv("sig_down_IJ.csv", index=False)
up_ie.to_csv("sig_up_IE.csv", index=False)
down_ie.to_csv("sig_down_IE.csv", index=False)    