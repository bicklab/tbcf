#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# In[2]:


bcf_cox = pd.read_csv("gs://bicklab-main-storage/Users/Hannah_Poisner/tcf_bcf_phewas/BCF_PheWAS_UKB_cox_all.csv")
tcf_cox = pd.read_csv("gs://bicklab-main-storage/Users/Hannah_Poisner/tcf_bcf_phewas/TCF_PheWAS_UKB_cox_all.csv")
ratio_cox = pd.read_csv('gs://bicklab-main-storage/Users/Hannah_Poisner/tcf_bcf_phewas/count_fraction_PheWAS_UKB_cox_all.csv')


# In[3]:


len(ratio_cox[ratio_cox['total_samples']>100])


# In[14]:


ratio_cox.head()


# In[4]:


ratio_sig = 0.05/len(ratio_cox[ratio_cox['total_samples']>100])
bcf_sig = 0.05/len(bcf_cox[bcf_cox['total_samples']>100])
tcf_sig = 0.05/len(tcf_cox[tcf_cox['total_samples']>100])


# # Sig in Both

# In[5]:


ratio_cox_sig = ratio_cox[(ratio_cox['p.value']<=ratio_sig) & (ratio_cox['total_samples']>100)]
bcf_cox_sig = bcf_cox[(bcf_cox['p.value']<=bcf_sig) & (bcf_cox['total_samples']>100)]
tcf_cox_sig = tcf_cox[(tcf_cox['p.value']<=tcf_sig) & (tcf_cox['total_samples']>100)]
ratio_cox_not_sig = ratio_cox[(ratio_cox['p.value']>ratio_sig) & (ratio_cox['total_samples']>100)]


# In[ ]:


len(ratio_cox_sig)


# In[7]:


bcf_cox[bcf_cox['disease']=='BI_170']


# In[ ]:





# In[6]:


ratio_tcf = ratio_cox_sig.merge(tcf_cox_sig, how='inner', on='disease')
len(ratio_tcf)


# In[7]:


ratio_tcf[(ratio_tcf['estimate_x']>0) & (ratio_tcf['estimate_y']<0)]


# In[8]:


custom_pal = sns.blend_palette(['lightgrey', "darkblue"], as_cmap=True, n_colors=5)
custom_pal


# In[11]:


ratio_tcf.head()


# In[ ]:


fig, ax = plt.subplots(figsize=(8, 7)) # width=10 inches, height=6 inches

ax =sns.scatterplot(data=ratio_tcf, x='estimate_x', y='estimate_y', size='log_p_value_y', hue='log_p_value_x', palette=custom_pal)
for index, row in ratio_tcf.iterrows():
    if (row['estimate_x'] >0) & (row['estimate_y'] <0):
        plt.text(row['estimate_x']-0.1, row['estimate_y']+.50, row['disease'], ha='center', va='bottom')
for index, row in ratio_tcf.iterrows():
    if(row['log_p_value_x']>47):
        plt.text(row['estimate_x']-.01, row['estimate_y']+.50, row['disease'], ha='center', va='bottom')
        
ax.axhline(y=0, color='black')
ax.axvline(x=0, color='black')
#plt.axis('equal')
plt.xlim(-17, 17)
plt.ylim(-17,17)
plt.xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=12, fontweight='bold')
plt.ylabel('TCF Effect Estimate', fontsize=12, fontweight='bold')
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_x':'Ratio -Log10(P)', 'log_p_value_y':'TCF -Log10(P)'}
replacer=my_map.get
labels2 = [replacer(n, n) for n in labels]
#ax.legend(handles, labels2)
sns.despine(bottom=True, left=True)
plt.tick_params(axis='x', which='both', bottom=False)

plt.tick_params(axis='y', which='both', left=False)
plt.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1,1), frameon=False)


# In[ ]:


from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(style="white", context="talk", font_scale=1.2)
fig, ax = plt.subplots(figsize=(8, 7))

# Scatterplot
ax = sns.scatterplot(
    data=ratio_tcf,
    x='estimate_x',
    y='estimate_y',
    size='log_p_value_y',
    hue='log_p_value_x',
    palette=custom_pal,
    sizes=(50, 400),
    alpha=0.8,
    edgecolor="black",
    linewidth=0.4
)

# Collect text labels for significant or interesting points
texts = []
for _, row in ratio_tcf.iterrows():
    # Example condition: opposite-direction effects or highly significant points
    if ((row['estimate_x'] > 0) & (row['estimate_y'] < 0)) or (row['log_p_value_x'] > 47):
        texts.append(
            ax.text(
                row['estimate_x'],
                row['estimate_y'],
                row['disease'],
                fontsize=10,
                fontweight='bold',
                color='black',
                ha='center',
                va='center',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=1)
            )
        )

# Automatically adjust label positions and add arrows
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.8),
    expand_points=(1.2, 1.2),
    force_points=0.5,
    force_text=0.5
)

# Reference lines
ax.axhline(y=0, color='gray', linestyle='--', linewidth=1)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)

# Axes limits and labels
ax.set_xlim(-17, 17)
ax.set_ylim(-17, 17)
ax.set_xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=16, fontweight='bold', labelpad=10)
ax.set_ylabel('T-cell Fraction Effect Estimate', fontsize=16, fontweight='bold', labelpad=10)

# Legend cleanup
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_x': 'Ratio -log₁₀(P)', 'log_p_value_y': 'TCF -log₁₀(P)'}
labels2 = [my_map.get(n, n) for n in labels]
ax.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1, 1), frameon=False, title="Significance")

# Style tweaks
sns.despine()
ax.tick_params(axis='both', which='major', length=0)
#plt.tight_layout()
plt.show()


# In[ ]:


ratio_tcf.head()


# In[12]:


from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(style="white", context="talk", font_scale=1.2)

fig, ax = plt.subplots(figsize=(8, 7))

# -------------------------------
# 1. Scatterplot (all black points)
# -------------------------------
ax.scatter(
    ratio_tcf['estimate_x'],
    ratio_tcf['estimate_y'],
    s=80,                    # fixed point size
    color='black',           # all black
    alpha=0.75,
    edgecolor='black',
    linewidth=0.4
)

# -------------------------------
# 2. Text labels for selected points
# -------------------------------
texts = []
for _, row in ratio_tcf.iterrows():
    if ((row['estimate_x'] > 0) & (row['estimate_y'] < 0)) or (row['log_p_value_x'] > 47):

        texts.append(
            ax.text(
                row['estimate_x'],
                row['estimate_y'],
                row['disease'],
                fontsize=10,
                fontweight='bold',
                color='black',
                ha='center',
                va='center',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=1)
            )
        )

adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.7),
    expand_points=(1.2, 1.2),
    force_points=0.5,
    force_text=0.5
)

# -------------------------------
# 3. Reference lines
# -------------------------------
ax.axhline(0, color='gray', linestyle='--', linewidth=1)
ax.axvline(0, color='gray', linestyle='--', linewidth=1)

# -------------------------------
# 4. Style and axes
# -------------------------------
ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)

ax.set_xlabel(
    r'$\beta_{\mathrm{Lymphocyte/WBC}}$',
    fontsize=18,
    fontweight='bold',
    labelpad=10
)
ax.set_ylabel(
    r'$\beta_{\mathrm{TCF}}$',
    fontsize=18,
    fontweight='bold',
    labelpad=10
)

sns.despine(top=True, right=True, left=False, bottom=False)  # keep left & bottom
ax.spines['left'].set_color('black')
ax.spines['bottom'].set_color('black')

ax.tick_params(axis='both', which='major', length=0)

plt.tight_layout()
plt.show()


# In[19]:


ratio_tcf.head()


# In[20]:


fig, ax = plt.subplots(figsize=(10, 8))

# -----------------------------
# POINT SIZE mapped to samples
# -----------------------------
sizes = (ratio_tcf['total_samples_x'] / ratio_tcf['total_samples_x'].max()) * 80

# -----------------------------
# SCATTER: all points black
# -----------------------------
ax.scatter(
    ratio_tcf['estimate_x'],
    ratio_tcf['estimate_y'],
    s=sizes,
    color='black',
    edgecolor='black',
    linewidth=0.6,
    alpha=0.8
)

# -----------------------------
# ERROR BARS: horizontal + vertical
# -----------------------------
ax.errorbar(
    ratio_tcf['estimate_x'],
    ratio_tcf['estimate_y'],
    xerr=ratio_tcf['std.error_x'],
    yerr=ratio_tcf['std.error_y'],
    fmt='none',
    ecolor='gray',
    elinewidth=1,
    capsize=3,
    zorder=0
)

# -----------------------------
# DIAGONAL REFERENCE LINE
# -----------------------------
lims = [-17, 17]
ax.plot(lims, lims, linestyle='--', color='gray', linewidth=1)

# -----------------------------
# GRID (light, thin) like example figure
# -----------------------------
ax.grid(True, color='lightgray', alpha=0.4, linewidth=0.7)

# -----------------------------
# AXIS LABELS (bold)
# -----------------------------
ax.set_xlabel(r'Lymphocyte/WBC $\beta$', fontsize=18, fontweight='bold')
ax.set_ylabel(r'TCF $\beta$', fontsize=18, fontweight='bold')

# ---------------------------------
# CLEAN SPINES: keep left/bottom
# ---------------------------------
sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)

# -----------------------------
# OPTIONAL: size legend (black)
# -----------------------------
for size_val in [20, 40, 60, 80]:
    ax.scatter([], [], s=size_val/80 * 80, color='black', label=f'N > {size_val}')

ax.legend(
    title="Sample Size",
    loc='upper left',
    bbox_to_anchor=(1.05, 1),
    frameon=False
)

plt.tight_layout()
plt.show()


# In[22]:


# Example cutpoints — adjust if needed
bins = [0, 300, 600, np.inf]
labels = ["Small (≤300)", "Medium (300–600)", "Large (>600)"]

ratio_tcf["sample_bin"] = pd.cut(
    ratio_tcf["total_samples_x"],  # or total_samples_y — whichever drives marker size
    bins=bins,
    labels=labels,
    include_lowest=True
)
size_map = {
    "Small (≤300)": 80,
    "Medium (300–600)": 160,
    "Large (>600)": 260
}

ratio_tcf["point_size"] = ratio_tcf["sample_bin"].map(size_map)


fig, ax = plt.subplots(figsize=(8, 7))

# Scatter with binned sizes
ax.scatter(
    ratio_tcf["estimate_x"], 
    ratio_tcf["estimate_y"],
    s=ratio_tcf["point_size"],
    color="black",
    alpha=0.7,
    edgecolor="black",
    linewidth=0.7
)

# Error bars (horizontal + vertical)
ax.errorbar(
    ratio_tcf["estimate_x"],
    ratio_tcf["estimate_y"],
    xerr=ratio_tcf["std.error_x"],
    yerr=ratio_tcf["std.error_y"],
    fmt="none",
    ecolor="black",
    elinewidth=0.8,
    capsize=3,
    alpha=0.8
)

# Zero-reference lines
ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)

# Labels
ax.set_xlabel(r"Lyphocyte/WBC $\beta$", fontsize=16, fontweight="bold")
ax.set_ylabel(r"TCF $\beta$", fontsize=16, fontweight="bold")

# Legend for sample-size bins (manual)
for label, size in size_map.items():
    ax.scatter([], [], s=size, color="black", label=label)

ax.legend(title="Sample Size (binned)", frameon=False, loc="upper left")

sns.despine(top=True, right=True)
plt.show()


# In[35]:


ratio_tcf


# In[8]:


from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# ---- 1. Define sample-size bins ----
bins = [0, 300, 600, np.inf]
labels = ["≤300", "300–600", ">600"]

ratio_tcf["sample_bin"] = pd.cut(
    ratio_tcf["total_samples_x"],   # or total_samples_y
    bins=bins,
    labels=labels,
    include_lowest=True
)

size_map = {
    "≤300": 80,
    "300–600": 160,
    ">600": 260
}

ratio_tcf["point_size"] = ratio_tcf["sample_bin"].map(size_map)

# ---- 2. Plot ----
plt.rcParams['svg.fonttype'] = 'none'
sns.set(style="white", context="paper", font_scale=1.2)
fig, ax = plt.subplots(figsize=(10, 7))

# Main scatter (black only)
ax.scatter(
    ratio_tcf["estimate_x"],
    ratio_tcf["estimate_y"],
    s=ratio_tcf["point_size"],
    color="black",
    alpha=0.50,
    edgecolor="black",
    linewidth=0.7
)

# ---- 3. Error bars ----
ax.errorbar(
    ratio_tcf["estimate_x"],
    ratio_tcf["estimate_y"],
    xerr=ratio_tcf["std.error_x"],
    yerr=ratio_tcf["std.error_y"],
    fmt="none",
    ecolor="black",
    elinewidth=0.8,
    capsize=3,
    alpha=0.8
)

# ---- 4. Add text labels (keep your conditions) ----

for _, row in ratio_tcf.iterrows():
    if ((row['estimate_x'] > 0) & (row['estimate_y'] < 0)) \
       or (row['log_p_value_x'] > 47) \
       or ((row['estimate_x'] < -5) & (row['estimate_y'] < -5)):

        # Determine offset direction based on quadrant
        # x_offset: move right for positive x, left for negative x
        x_offset = 0 if row["estimate_x"] >= 0 else 5
        # y_offset: move up for positive y, down for negative y
        y_offset = 5 if row["estimate_y"] >= 0 else 5

        ax.annotate(
            row["disease"],
            xy=(row["estimate_x"], row["estimate_y"]),
            xytext=(x_offset, y_offset),        # offset in points
            textcoords='offset points',
            fontsize=14,
            fontweight='bold',
            ha='center' if x_offset==0 else ('left' if x_offset>0 else 'right'),
            va='center' if y_offset==0 else ('bottom' if y_offset>0 else 'top'),
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.1, pad=1),
            arrowprops=dict(arrowstyle='->', color='gray', lw=0.8)
        )


# ---- 5. Reference lines ----
ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)

# ---- 6. Axes labels ----
ax.set_xlabel('\u03B2 (Lymphocyte/WBC)', fontsize=20, fontweight="bold")
ax.set_ylabel('\u03B2 (T-Cell Fraction)', fontsize=20, fontweight="bold")

ax.set_xlim(-14, 14)
ax.set_ylim(-14, 14)

# ---- 7. Right-Side Legend (binned sample size) ----
for label, size in size_map.items():
    ax.scatter([], [], s=size, color="black", label=label)

legend = ax.legend(
    title="Case Number",
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
    labelspacing=1.0,    # default 0.5 → increase to spread out entries
    handlelength=2.5,    # length of the marker
    handleheight=1.5     # height of the marker box
)
ax.tick_params(axis='both', labelsize=18)

sns.despine(top=True, right=True, bottom=True, left=True)
plt.tight_layout()
#plt.show()
fig.savefig('F5A_tcf_lymphpercentage.svg', format='svg', dpi=300)


# In[9]:


ratio_bcf = ratio_cox_sig.merge(bcf_cox_sig, how='inner', on='disease')
len(ratio_bcf)


# In[19]:


from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# ---- 1. Define sample-size bins ----
bins = [0, 300, 600, np.inf]
labels = ["≤300", "300–600", ">600"]

ratio_bcf["sample_bin"] = pd.cut(
    ratio_bcf["total_samples_x"],  # or total_samples_y
    bins=bins,
    labels=labels,
    include_lowest=True
)

size_map = {
    "≤300": 80,
    "300–600": 160,
    ">600": 260
}

ratio_bcf["point_size"] = ratio_bcf["sample_bin"].map(size_map)

# ---- 2. Setup figure ----
plt.rcParams['svg.fonttype'] = 'none'  # keep text editable in SVG
sns.set(style="white", context="talk", font_scale=1.2)

fig, ax = plt.subplots(figsize=(10, 7))

# ---- 3. Scatter points ----
ax.scatter(
    ratio_bcf["estimate_x"],
    ratio_bcf["estimate_y"],
    s=ratio_bcf["point_size"],
    color="black",
    alpha=0.5,
    edgecolor="black",
    linewidth=0.7
)

# ---- 4. Error bars (optional) ----
ax.errorbar(
    ratio_bcf["estimate_x"],
    ratio_bcf["estimate_y"],
    xerr=ratio_bcf["std.error_x"],
    yerr=ratio_bcf["std.error_y"],
    fmt="none",
    ecolor="black",
    elinewidth=0.8,
    capsize=3,
    alpha=0.8
)

# ---- 5. Add labels with arrows ----
texts = []
x = []
y = []

for _, row in ratio_bcf.iterrows():
    if row['estimate_y'] > 8:  # condition for labeling
        txt = ax.text(
            row["estimate_x"],
            row["estimate_y"],
            row["disease"],
            fontsize=14,
            fontweight='bold',
            ha='left',
            va='bottom',
            color='black',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.2, pad=1)
        )
        texts.append(txt)
        x.append(row["estimate_x"])
        y.append(row["estimate_y"])

# Adjust text to avoid overlap, with arrows
adjust_text(
    texts,
    x=x,
    y=y,
    arrowprops=dict(arrowstyle='->', color='k', lw=0.8),
    expand_points=(1.2, 1.2),
    expand_text=(1.2, 1.2),
    force_points=0.5,
    force_text=0.5,
    lim=2000
)

# ---- 6. Reference lines ----
ax.axhline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)
ax.axvline(0, color="black", linestyle="--", linewidth=1, alpha=0.5)

# ---- 7. Axes labels ----
ax.set_xlabel('\u03B2 (Lymphocyte/WBC)', fontsize=20, fontweight="bold")
ax.set_ylabel('\u03B2 (B-Cell Fraction)', fontsize=20, fontweight="bold")

ax.set_xlim(-14, 14)
ax.set_ylim(-14, 14)

# ---- 8. Legend for sample sizes ----
for label, size in size_map.items():
    ax.scatter([], [], s=size, color="black", label=label)

legend = ax.legend(
    title="Case Number",
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
    labelspacing=0.2,    # default 0.5 → increase to spread out entries
    handlelength=2.5,    # length of the marker
    handleheight=1.5     # height of the marker box
)

# ---- 9. Final tweaks ----
ax.tick_params(axis='both', labelsize=18)

sns.despine(top=True, right=True, bottom=True, left=True)
plt.tight_layout()
plt.show()

# ---- 10. Save figure if desired ----
fig.savefig('F5B_bcf_lymphpercentage.svg', format='svg', dpi=300)


# In[ ]:


ratio_bcf[(ratio_bcf['estimate_x']>0) & (ratio_bcf['estimate_y']<0)]


# # Not Sig Analysis

# In[24]:


not_ratio_tcf = tcf_cox_sig.merge(ratio_cox_not_sig, how='inner', on='disease')
not_ratio_tcf['disease']


# In[ ]:





# In[25]:


not_ratio_tcf.head(1)


# In[26]:


not_ratio_bcf = bcf_cox_sig.merge(ratio_cox_not_sig, how='left', on='disease').dropna(subset=['term_y'])
not_ratio_bcf


# In[ ]:


ratio_sig


# In[ ]:


tcf_sig


# In[10]:


import matplotlib_inline
matplotlib_inline.__version__


# In[ ]:


format_pval(5.208333333333334e-05)


# In[27]:


import pandas as pd
import matplotlib.pyplot as plt
import forestplot as fp
import numpy as np

# Sample data
df = pd.DataFrame({
    "Variable": ["TCF Sarcoidosis", "Lymph/WBC Sarcoidosis"],
    "Estimate": [-6.6, -4.43],
    "SE": [1.3, 1.295],  # Standard Error
    "P-value": [3.946027e-07, 0.000623]  # P-values
})

# Compute 95% Confidence Intervals
df["Lower CI"] = df["Estimate"] - 1.96 * df["SE"]
df["Upper CI"] = df["Estimate"] + 1.96 * df["SE"]
# Create forest plot with p-values

# Plot settings
# Plot settings
fig, ax = plt.subplots(figsize=(8, 1))
y_pos = np.arange(len(df))
ax.text(df["Upper CI"].max() + 0.8, max(y_pos) + 0.5, "Beta [95% CI]", fontsize=20, fontweight="bold", ha="left")
ax.text(df["Upper CI"].max() + 5.6, max(y_pos) + 0.5, "P-value", fontsize=20, fontweight="bold", ha="left")

# Plot estimates as points
ax.scatter(df["Estimate"], y_pos, color="black")

# Add confidence intervals as horizontal lines
for i in range(len(df)):
    ax.plot([df["Lower CI"][i], df["Upper CI"][i]], [y_pos[i], y_pos[i]], color="black")

# Move variable labels to the left
for i in range(len(df)):
    ax.text(df["Lower CI"].min() - 4.9, y_pos[i], f"{df['Variable'][i]}", fontsize=15, fontweight='bold', verticalalignment="center", ha="left")

# Move Beta & CI values to the right
for i in range(len(df)):
    ax.text(df["Upper CI"].max() + 0.8, y_pos[i], f"{df['Estimate'][i]:.2f} [{df['Lower CI'][i]:.2f}, {df['Upper CI'][i]:.2f}]", fontsize=15, verticalalignment="center", ha="left")

# Move p-values further outside the plot
for i, pval in enumerate(df["P-value"]):
    significance = "*" if pval < 5.208333333333334e-05 else ""  # Add a star for p < 0.05
    ax.text(df["Upper CI"].max() + 5.6, y_pos[i], f"{pval:.2e}{significance}", fontsize=15, verticalalignment="center", ha="left")

# Remove y-ticks and legend
ax.set_yticks([])
ax.set_xlabel("Effect Size", fontsize=20, fontweight='bold')
ax.tick_params(labelsize=15)
# Remove all spines for a clean look
for spine in ax.spines.values():
    spine.set_visible(False)


plt.tight_layout()
plt.show()


# In[ ]:


df


# # All

# In[11]:


ratio_cox_100 = ratio_cox[ratio_cox['total_samples']>100]
bcf_cox_100 = bcf_cox[bcf_cox['total_samples']>100]
tcf_cox_100 = tcf_cox[tcf_cox['total_samples']>100]


# In[12]:


ratio_bcf_all = ratio_cox_100.merge(bcf_cox_100, how='inner', on='disease').dropna()

ratio_tcf_all = ratio_cox_100.merge(tcf_cox_100, how='inner', on='disease').dropna()


# In[13]:


custom_pal2 = sns.blend_palette(['lightgrey', "darkblue"], as_cmap=True, n_colors=6)
custom_pal2


# In[ ]:


fig, ax = plt.subplots(figsize=(9, 8)) # width=10 inches, height=6 inches

ax =sns.scatterplot(data=ratio_bcf_all, x='estimate_x', y='estimate_y', size='log_p_value_y', hue='log_p_value_x', palette=custom_pal2)

for index, row in ratio_bcf_all.iterrows():
    if(row['log_p_value_y']>250):
        plt.text(row['estimate_x'], row['estimate_y']+.70, row['disease'], ha='center', va='bottom')



ax.axhline(y=0, color='black')
ax.axvline(x=0, color='black')
plt.xlim(-15, 15)
plt.ylim(-15,15)

plt.xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=12, fontweight='bold')
plt.ylabel('BCF Effect Estimate', fontsize=12, fontweight='bold')
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_x':'Ratio -log10(P)', 'log_p_value_y':'BCF -log10(P)'}
replacer=my_map.get
labels2 = [replacer(n, n) for n in labels]
#ax.legend(handles, labels2)
sns.despine(bottom=True, left=True)
plt.tick_params(axis='x', which='both', bottom=False)

plt.tick_params(axis='y', which='both', left=False)
plt.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1,1), frameon=False)


# In[ ]:


ratio_tcf_all.sort_values(by='log_p_value_x')


# In[ ]:


fig, ax = plt.subplots(figsize=(9, 8)) # width=10 inches, height=6 inches

ax =sns.scatterplot(data=ratio_tcf_all, x='estimate_x', y='estimate_y', size='log_p_value_y', hue='log_p_value_x', palette=custom_pal2)

for index, row in ratio_tcf_all.iterrows():
    if(row['log_p_value_x']>47):
        plt.text(row['estimate_x']-.01, row['estimate_y']+.50, row['disease'], ha='center', va='bottom')


ax.axhline(y=0, color='black')
ax.axvline(x=0, color='black')
plt.xlim(-15, 15)
plt.ylim(-15,15)

plt.xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=12, fontweight='bold')
plt.ylabel('TCF Effect Estimate', fontsize=12, fontweight='bold')
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_x':'Ratio -log10(P)', 'log_p_value_y':'TCF -log10(P)'}
replacer=my_map.get
labels2 = [replacer(n, n) for n in labels]
#ax.legend(handles, labels2)
sns.despine(bottom=True, left=True)
plt.tick_params(axis='x', which='both', bottom=False)

plt.tick_params(axis='y', which='both', left=False)
plt.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1,1), frameon=False)


# In[ ]:


not_ratio_bcf


# In[ ]:


fig, ax = plt.subplots(figsize=(9, 8)) # width=10 inches, height=6 inches

ax =sns.scatterplot(data=not_ratio_bcf, x='estimate_x', y='estimate_y', size='log_p_value_x', hue='log_p_value_y', palette=custom_pal)

for index, row in not_ratio_bcf.iterrows():
    if(row['log_p_value_x']>7):
        plt.text(row['estimate_x']-.01, row['estimate_y']+.30, row['disease'], ha='center', va='bottom')



ax.axhline(y=0, color='black')
ax.axvline(x=0, color='black')
plt.xlim(-15, 15)
plt.ylim(-15,15)

plt.xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=12, fontweight='bold')
plt.ylabel('BCF Effect Estimate', fontsize=12, fontweight='bold')
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_y':'Ratio -log10(P)', 'log_p_value_x':'BCF -log10(P)'}
replacer=my_map.get
labels2 = [replacer(n, n) for n in labels]
#ax.legend(handles, labels2)
sns.despine(bottom=True, left=True)
plt.tick_params(axis='x', which='both', bottom=False)

plt.tick_params(axis='y', which='both', left=False)
plt.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1,1), frameon=False)


# In[ ]:


fig, ax = plt.subplots(figsize=(9, 8)) # width=10 inches, height=6 inches

ax =sns.scatterplot(data=not_ratio_tcf, x='estimate_x', y='estimate_y', size='log_p_value_x', hue='log_p_value_y', palette=custom_pal)

for index, row in not_ratio_tcf.iterrows():
    if(row['log_p_value_x']>5):
        plt.text(row['estimate_x']-.01, row['estimate_y']+.30, row['disease'], ha='center', va='bottom')



ax.axhline(y=0, color='black')
ax.axvline(x=0, color='black')
plt.xlim(-15, 15)
plt.ylim(-15,15)

plt.xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=12, fontweight='bold')
plt.ylabel('TCF Effect Estimate', fontsize=12, fontweight='bold')
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_y':'Ratio -log10(P)', 'log_p_value_x':'TCF -log10(P)'}
replacer=my_map.get
labels2 = [replacer(n, n) for n in labels]
#ax.legend(handles, labels2)
sns.despine(bottom=True, left=True)
plt.tick_params(axis='x', which='both', bottom=False)

plt.tick_params(axis='y', which='both', left=False)
plt.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1,1), frameon=False)


# In[ ]:


from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(style="white", context="talk", font_scale=1.2)
fig, ax = plt.subplots(figsize=(9, 8))

# Round the columns used for size and hue
df = not_ratio_tcf.copy()
df['log_p_value_x_round'] = df['log_p_value_x'].round(1)
df['log_p_value_y_round'] = df['log_p_value_y'].round(1)

# Scatterplot using rounded values for legend
ax = sns.scatterplot(
    data=df,
    x='estimate_x',
    y='estimate_y',
    size='log_p_value_x_round',
    hue='log_p_value_y_round',
    palette=custom_pal,
    sizes=(50, 400),
    alpha=0.8,
    edgecolor="black",
    linewidth=0.4
)

# Labels for significant points
texts = []
for _, row in df.iterrows():
    if row['log_p_value_x'] > 5:
        texts.append(
            ax.text(
                row['estimate_x'],
                row['estimate_y'],
                row['disease'],
                fontsize=10,
                fontweight='bold',
                color='black',
                ha='center',
                va='center',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=1)
            )
        )

# Automatically adjust labels and add arrows
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.8),
    expand_points=(1.2, 1.2),
    force_points=0.5,
    force_text=0.5
)

# Reference lines
ax.axhline(y=0, color='gray', linestyle='--', linewidth=1)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)

# Axes limits and labels
ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)
ax.set_xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=16, fontweight='bold')
ax.set_ylabel('T-cell Fraction Effect Estimate', fontsize=16, fontweight='bold')

# Rounded legend labels
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_y_round': 'Ratio -log₁₀(P)', 'log_p_value_x_round': 'T-cell Fraction -log₁₀(P)'}
labels2 = [my_map.get(n, n) for n in labels]
ax.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1, 1), frameon=False, title="Significance")

# Style cleanup
sns.despine()
ax.tick_params(axis='both', which='major', length=0)
#plt.tight_layout()
plt.show()


# In[ ]:


from adjustText import adjust_text
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set(style="white", context="talk", font_scale=1.2)
fig, ax = plt.subplots(figsize=(9, 8))

# Create rounded copies for legend
df = not_ratio_bcf.copy()
df['log_p_value_x_round'] = df['log_p_value_x'].round(1)
df['log_p_value_y_round'] = df['log_p_value_y'].round(1)

# Scatterplot using rounded values for size and hue
ax = sns.scatterplot(
    data=df,
    x='estimate_x',
    y='estimate_y',
    size='log_p_value_x_round',   # rounded for size legend
    hue='log_p_value_y_round',    # rounded for color legend
    palette=custom_pal,
    sizes=(50, 400),
    alpha=0.8,
    edgecolor="black",
    linewidth=0.4
)

# Labels for significant points
texts = []
for _, row in df.iterrows():
    if row['log_p_value_x'] > 7:
        texts.append(
            ax.text(
                row['estimate_x'],
                row['estimate_y'],
                row['disease'],
                fontsize=10,
                fontweight='bold',
                color='black',
                ha='center',
                va='center',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.6, pad=1)
            )
        )

# Automatically adjust labels and add arrows
adjust_text(
    texts,
    arrowprops=dict(arrowstyle='-', color='gray', lw=0.8),
    expand_points=(1.2, 1.2),
    force_points=0.5,
    force_text=0.5
)

# Reference lines
ax.axhline(y=0, color='gray', linestyle='--', linewidth=1)
ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)

# Axes and labels
ax.set_xlim(-15, 15)
ax.set_ylim(-15, 15)
ax.set_xlabel('Lymphocyte/WBC Ratio Effect Estimate', fontsize=16, fontweight='bold')
ax.set_ylabel('B-cell Fraction Effect Estimate', fontsize=16, fontweight='bold')

# Rounded legend labels
handles, labels = ax.get_legend_handles_labels()
my_map = {'log_p_value_y_round': 'Ratio -log₁₀(P)', 'log_p_value_x_round': 'B-Cell Fraction -log₁₀(P)'}
labels2 = [my_map.get(n, n) for n in labels]
ax.legend(handles, labels2, loc="upper left", bbox_to_anchor=(1, 1), frameon=False, title="Significance")

# Style tweaks
sns.despine()
ax.tick_params(axis='both', which='major', length=0)
#plt.tight_layout()
plt.show()


# In[ ]:


tcf_cox.sort_values(by=['log_p_value'],ascending=True).tail(35)


# # Ratio Effect Comparison

# In[14]:


ratio_cox_sig_tcf_merge = ratio_cox_sig.merge(tcf_cox, how='left', on='disease')
ratio_cox_sig_cfs_merge = ratio_cox_sig_tcf_merge.merge(bcf_cox, how='left', on='disease').dropna()
ratio_cox_sig_cfs_merge.head()


# In[15]:


ratio_cox_sig_cfs_merge['effect_bias'] = ratio_cox_sig_cfs_merge['estimate_y'] - ratio_cox_sig_cfs_merge['estimate']
ratio_cox_sig_cfs_merge.head()


# In[16]:


strongest_tcf_bias = ratio_cox_sig_cfs_merge.sort_values(by='effect_bias', ascending=False).head(10)
strongest_bcf_bias = ratio_cox_sig_cfs_merge.sort_values(by='effect_bias', ascending=False).tail(10)
strongest_bias = pd.concat([strongest_bcf_bias, strongest_tcf_bias])
strongest_bias.head()


# In[17]:


strongest_bias


# In[18]:


strongest_bias.replace("DE_686", 'Chronic ulcer of skin', inplace=True)


# In[30]:


ratio_sig


# In[33]:


strongest_bias[(strongest_bias['disease'].isin(df_long_keep)) & (strongest_bias['p.value_x']<ratio_sig)]


# In[19]:


strongest_bias.rename(columns={'estimate_y':"beta_tcf", 'estimate':'beta_bcf', 'std.error_y':'se_tcf', 'std.error':'se_bcf',}, inplace=True)
df_long_beta = strongest_bias.melt(id_vars='disease', value_vars=['beta_tcf', 'beta_bcf'], var_name='source', value_name='beta')
df_long_se = strongest_bias.melt(id_vars='disease', value_vars=['se_tcf', 'se_bcf'], var_name='source', value_name='se')
df_long=pd.concat([df_long_beta, df_long_se[['se']]], axis=1)


# In[27]:


df_long.head()


# In[21]:


df_long_keep = ['Essential hypertension','Rheumatoid arthritis','Chronic lymphoid leukemia','Gout','Polycythemias','Essential thrombocythemia',]


# In[22]:


strongest_bias_fil = strongest_bias[strongest_bias['disease'].isin(df_long_keep)]


# In[25]:


df_long_fil = df_long[df_long['disease'].isin(df_long_keep)]
df_long_fil.to_csv('data_for_het_test_for_kun.txt',index=False,sep='\t')


# In[ ]:


import matplotlib.pyplot as plt

# Sorting for visual clarity
df_long_fil['disease'] = pd.Categorical(df_long_fil['disease'], categories=strongest_bias_fil['disease'], ordered=True)

fig, ax = plt.subplots(figsize=(8, 6))
sources = df_long_fil['source'].unique()
bar_width = 0.4
y_pos = range(0,6)
print(y_pos)

for i, source in enumerate(sources):
    offset = -bar_width/2 if source == 'beta_tcf' else bar_width/2
    subset = df_long_fil[df_long_fil['source'] == source]
    ax.barh(
        [y + offset for y in y_pos],
        subset['beta'],
        height=bar_width,
        label=source.upper().replace('_', ' '),
        color=("#1C4786" if "tcf" in source else "#E69F00")
    )

ax.set_yticks(y_pos)
ax.set_yticklabels(strongest_bias_fil['disease'])
ax.axvline(0, color='gray', linestyle='--')
ax.tick_params(axis='both', which='major', labelsize=16)
ax.spines[['top', 'right','left']].set_visible(False)  # cleaner frame
ax.set_xlabel("Beta Coefficient", fontweight='bold', fontsize=16)
ax.legend(title="Source")
plt.tight_layout()
plt.show()


# In[ ]:


df_long[df_long['disease']=='Hypertension']


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt

# Ensure proper categorical ordering based on strongest_beta values
df_long_fil['disease'] = pd.Categorical(
    df_long_fil['disease'],
    categories=strongest_bias_fil['disease'],
    ordered=True
)

# Plot setup
fig, ax = plt.subplots(figsize=(8, 6))
sources = df_long_fil['source'].unique()
y_pos = range(len(strongest_bias_fil['disease']))
point_offset = 0.15  # side-by-side spacing of points

# Loop through each source and plot points with 95% CI error bars
for i, source in enumerate(sources):
    offset = -point_offset if source == 'beta_tcf' else point_offset
    subset = df_long_fil[df_long_fil['source'] == source].copy()
    
    # Compute y positions with offset
    y_vals = [y + offset for y in y_pos]
    
    # Plot the points
    ax.scatter(
        subset['beta'],
        y_vals,
        label=source.replace({'beta_bcf':'B-cell Fraction', 'beta_tcf':'T-cell Fraction'}.items()),
        color="#1C4786" if "tcf" in source else "#E69F00",
        s=100,
        edgecolor='black',
        zorder=3
    )
    
    # Plot 95% confidence interval error bars
    ci = 1.96 * subset['se']
    ax.errorbar(
        x=subset['beta'],
        y=y_vals,
        xerr=ci,
        fmt='none',
        ecolor='black',
        elinewidth=1,
        capsize=3,
        zorder=2
    )

# Y-axis formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(strongest_bias_fil['disease'], fontweight='bold')
ax.axvline(0, color='gray', linestyle='--', linewidth=1)
ax.tick_params(axis='both', labelsize=14)
ax.spines[['top', 'right', 'left']].set_visible(False)

# Labels and legend
ax.set_xlabel("Beta Coefficient", fontweight='bold', fontsize=16)
ax.legend(title="Phenotype", fontsize=14, title_fontsize=16, loc='upper right', frameon=False)

plt.tight_layout()
plt.show()


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Ensure proper categorical ordering based on strongest_beta values
df_long_fil['disease'] = pd.Categorical(
    df_long_fil['disease'],
    categories=strongest_bias_fil['disease'],
    ordered=True
)

# Define color palette
COLORS = {
    'tcf': '#1C4786',  # Blue for T-cell
    'bcf': '#E69F00'   # Orange for B-cell
}

LABELS = {
    'beta_tcf': 'T-cell Fraction',
    'beta_bcf': 'B-cell Fraction'
}

# Plot setup with larger figure size for poster
fig, ax = plt.subplots(figsize=(12, 10))

sources = df_long_fil['source'].unique()
y_pos = np.arange(len(strongest_bias_fil['disease']))
point_offset = 0.2  # Adjusted spacing for clarity

# Loop through each source and plot points with 95% CI error bars
for i, source in enumerate(sources):
    offset = -point_offset if source == 'beta_tcf' else point_offset
    subset = df_long_fil[df_long_fil['source'] == source].copy()
    
    # Compute y positions with offset
    y_vals = y_pos + offset
    
    # Determine color based on source
    color = COLORS['tcf'] if 'tcf' in source else COLORS['bcf']
    label = LABELS.get(source, source)
    
    # Plot the points
    ax.scatter(
        subset['beta'],
        y_vals,
        label=label,
        color=color,
        s=200,  # Larger points for poster visibility
        edgecolor='black',
        linewidth=1.5,
        zorder=3,
        alpha=0.9
    )
    
    # Plot 95% confidence interval error bars
    ci = 1.96 * subset['se']
    ax.errorbar(
        x=subset['beta'],
        y=y_vals,
        xerr=ci,
        fmt='none',
        ecolor='black',
        elinewidth=2,  # Thicker error bars for visibility
        capsize=5,
        capthick=2,
        zorder=2,
        alpha=0.8
    )

# Add vertical line at zero
ax.axvline(0, color='#666666', linestyle='--', linewidth=2, zorder=1)

# Y-axis formatting
ax.set_yticks(y_pos)
ax.set_yticklabels(strongest_bias_fil['disease'], fontsize=20, fontweight='bold')
ax.set_ylim(-0.5, len(y_pos) - 0.5)

# X-axis formatting
ax.set_xlabel("Effect Size (β)", fontsize=22, fontweight='bold', labelpad=15)
ax.tick_params(axis='x', labelsize=20, width=2, length=6)
ax.tick_params(axis='y', width=0)  # Remove y-axis ticks

# Remove spines for cleaner look
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_linewidth(2)


# Legend formatting
legend = ax.legend(
    title="Cell Type Association",
    fontsize=20,
    title_fontsize=23,
    bbox_to_anchor=(0.70, 1),  # (x, y) - adjust x value to move right
    frameon=False,
    edgecolor='black',
    markerscale=1.2
)
legend.get_frame().set_linewidth(2)
legend.get_title().set_fontweight('bold')


# Adjust layout
plt.tight_layout()
plt.rcParams['svg.fonttype'] = 'none'  # Ensures text stays as text, not paths

# Save as high-resolution files for poster
plt.savefig('phewas_bias_forest_plot_poster.svg', format='svg', bbox_inches='tight', facecolor='white')

plt.show()


# In[ ]:


df_long_fil['ci'] = 1.96 * df_long_fil['se']
df_long_fil['low'] = df_long_fil['beta'] - df_long_fil['ci']
df_long_fil['high'] = df_long_fil['beta'] + df_long_fil['ci']
hper_only = df_long_fil[df_long_fil['disease']=='Essential hypertension']


# In[ ]:


df_long_fil


# In[ ]:




