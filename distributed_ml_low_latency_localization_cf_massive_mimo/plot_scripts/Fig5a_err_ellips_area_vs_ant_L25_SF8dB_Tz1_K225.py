import matplotlib.pyplot as plt
color_map = {
    'Centralized RSS': 'tab:blue',
    'Centralized Hybrid': 'tab:orange',
    'Centralized AOA': 'tab:green',
    'Distributed Median': 'tab:red',
    'Distributed Mean': 'tab:purple',
    'Distributed Bayesian': 'tab:brown',
    'Distributed Zscore': 'tab:pink',

}

# Antenna numbers
antennas = [3, 4, 5, 6, 8, 10, 12, 15, 17, 19, 21, 25, 30, 35, 40, 45, 50]

centralised_rss = [
    1561.6627, 1561.0479, 1561.5101, 1559.8311, 1558.3192, 1556.4082, 1558.6975,
    1557.9748, 1563.3138, 1559.6990, 1554.6534, 1562.5564, 1560.8379, 1557.9924,
    1562.8178, 1565.2289, 1558.5472
]

centralised_aoa = [
    2133.8276, 1801.7773, 1588.3021, 1442.1122, 1272.7654, 1135.3571, 1036.8017,
    958.5679, 904.9014, 869.5939, 848.8472, 801.6391, 756.5806, 737.3459, 724.2714,
    695.8741, 690.2187
]

centralised_hybrid = [
    1449.8962, 1232.3250, 1085.4959, 984.8125, 873.2705, 782.3099, 715.8262,
    660.6409, 627.4929, 607.0640, 589.4853, 559.4591, 528.5258, 515.3087, 505.7269,
    485.5603, 481.9775
]

distributed_median = [
    6295.4739, 6273.3537, 6250.9692, 6246.1864, 6230.1560, 6231.0322, 6219.5267,
    6232.6267, 6236.4823, 6242.7825, 6243.7959, 6257.9727, 6246.7758, 6248.6076,
    6246.2263, 6252.6987, 6255.4195
]

distributed_mean = [
    342.6563, 334.2278, 329.4517, 325.3443, 321.3055, 318.0320, 316.2871,
    314.0033, 312.4718, 312.2500, 311.5312, 310.1729, 309.2882, 308.7666, 308.1899,
    307.8287, 307.8016
]
distributed_bayesian = [
    244.8340, 244.1447, 243.9524, 243.5138, 243.3335, 243.0125, 242.9883, 
    242.8006, 242.7055, 242.7524, 242.6358, 242.7065, 242.4597, 242.4903, 
    242.5814, 242.3883, 242.4579
]

distributed_zscore = [
    309.5185, 310.3237, 311.2927, 311.5519, 312.6158, 313.0130, 313.6768, 
    314.1892, 314.5131, 315.0611, 315.2544, 315.8543, 315.9571, 316.2788, 
    316.7245, 316.7827, 316.9773
]


# Create figure and axes
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6), 
                               gridspec_kw={'height_ratios': [1.5, 8.5]})
fig.subplots_adjust(hspace=0.05)

# Plot on both axes
for ax in [ax1, ax2]:
    ax.plot(antennas, centralised_rss, marker='^', label='Centralized-RSS', markersize=4,
            color=color_map['Centralized RSS'])
    ax.plot(antennas, centralised_hybrid, marker='^', label='Centralized-Hybrid',
            color=color_map['Centralized Hybrid'])
    ax.plot(antennas, centralised_aoa, marker='^', label='Centralized-AOA', markersize=5,
            color=color_map['Centralized AOA'])
    ax.plot(antennas, distributed_median, marker='o', label='Distributed-Median', markersize=4.5,
            color=color_map['Distributed Median'])
    ax.plot(antennas, distributed_bayesian, marker='o', label='Distributed-Bayesian', markersize=5,
            color=color_map['Distributed Bayesian'])
    ax.plot(antennas, distributed_mean, marker='o', label='Distributed-Mean',
        color=color_map['Distributed Mean'])  # 0.0 = fully transparent, 1.0 = fully opaque
    ax.plot(antennas, distributed_zscore, marker='o', label='Distributed-Z-score', markersize=7,
            color=color_map['Distributed Zscore'],alpha=0.5)


# Grid and axis adjustments
ax1.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
ax2.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)

ax1.set_ylim(5800, 7000)
ax2.set_ylim(0, 2500)
ax2.set_yticks([0, 500, 1000, 1500, 2000,2500])

ax2.set_xlim(3, 50)
ax2.set_xticks([3] + list(range(5, 55, 5)))

ax1.xaxis.set_tick_params(labelsize=15)
ax1.yaxis.set_tick_params(labelsize=15)

ax2.xaxis.set_tick_params(labelsize=15)
ax2.yaxis.set_tick_params(labelsize=15)

ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Axis breaks
d = 0.5
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels and legend
fig.text(0.5, 0, 'Number of Antennas per AP (N)', ha='center', fontsize=23)
fig.text(0.02, 0.5, r'Area of 95% Error Ellipse (m$^\boldsymbol{2}$)', va='center', rotation='vertical', fontsize=23)

ax2.legend(fontsize=15, loc='upper right',bbox_to_anchor=(0.98, 1.13),framealpha=0.4)

# Save and show
plt.savefig('Fig5a_err_ellips_area_vs_ant_L25_SF8dB_Tz1_K225.png', dpi=300, bbox_inches='tight')
plt.show()
