import matplotlib.pyplot as plt

color_map = {
    'Centralized RSS': 'tab:blue',
    'Centralized Hybrid': 'tab:orange',
    'Centralized AOA': 'tab:green',
    'Distributed Median': 'tab:red',
    'Distributed Mean': 'tab:purple',
    'Distributed Bayesian': 'tab:brown',
    'Distributed Zscore (Tz=1)': 'tab:pink',
    'Centralized RSS KNN': 'tab:gray',
    'Distributed Median KNN': 'tab:olive',
    'Centralized RSS LR': 'tab:cyan',
    'Centr Hybrid (no train noise)': 'gold',
    'Distr Bayesian (no train noise)': 'darkred',
    'Distr Zscore (Tz=1, no train noise)': 'deeppink',
    'FCNN': 'indigo'
}

ant_number = [3, 4, 5, 6, 8, 10, 12, 15, 17, 19, 21, 25, 30, 35, 40, 45, 50]

centralised_rss = [
    14.1259, 14.1379, 14.1322, 14.1353, 14.1482, 14.1315, 14.1255, 
    14.1318, 14.1321, 14.1376, 14.1452, 14.1151, 14.1260, 14.1349, 
    14.1443, 14.1442, 14.1294
]

centralised_hybrid = [
    13.5224, 12.3249, 11.5504, 11.0129, 10.3331, 9.7864, 9.5346, 
    9.1522, 8.9985, 8.9135, 8.8019, 8.6299, 8.5152, 8.3917, 
    8.3918, 8.2641, 8.2184
]

centralised_aoa = [
    14.6627, 13.1957, 12.4632, 11.7752, 10.9834, 10.3962, 10.1261, 
    9.7353, 9.5042, 9.4813, 9.2729, 9.0767, 9.0291, 8.8606, 
    8.8386, 8.7055, 8.6935
]

distributed_median = [
    7.0502, 6.9760, 6.9487, 6.9166, 6.8938, 6.8896, 6.8803, 
    6.8971, 6.9147, 6.9235, 6.9321, 6.9456, 6.9377, 6.9464, 
    6.9411, 6.9493, 6.9385
]

distributed_mean = [
    8.7218, 8.3963, 8.1927, 8.0819, 7.9105, 7.8197, 7.7520, 
    7.6787, 7.6554, 7.6365, 7.6106, 7.5929, 7.5623, 7.5466, 
    7.5323, 7.5255, 7.5135
]

distributed_bayesian =  [
    8.2994, 8.2098, 8.1598, 8.1382, 8.0962, 8.0757, 8.0534, 8.0462,
    8.0352, 8.0443, 8.0471, 8.0411, 8.0260, 8.0288, 8.0418, 8.0357, 8.0302
]

distributed_zscore_tz_1 = [
    7.5074, 7.4449, 7.4231, 7.3932, 7.3630, 7.3461, 7.3316, 7.3357,
    7.3333, 7.3527, 7.3581, 7.3596, 7.3441, 7.3464, 7.3597, 7.3605, 7.3622
]


knn_rss = [
    10.7618, 10.7662, 10.7759, 10.7621, 10.7635, 10.7492, 10.7508, 
    10.7552, 10.7493, 10.7492, 10.7574, 10.7477, 10.7512, 10.7480, 
    10.7456, 10.7423, 10.7528
]

lr_rss = [25.1954, 25.2593, 25.2295, 25.1994, 25.2054, 25.2101, 25.2071, 25.2012, 25.2083, 
               25.1772, 25.2186, 25.1940, 25.1849, 25.1823, 25.2028, 25.1822, 25.1793]


distributed_knn = [7.1564, 7.0779, 7.0601, 7.0286, 7.0245, 7.0027, 7.0080, 7.0521, 7.0497, 
                    7.0867, 7.0945, 7.1084, 7.1114, 7.1301, 7.1320, 7.1499, 7.1344]

distributed_lr = [46.8575, 46.8756, 46.8776, 46.8938, 46.8931, 46.9008, 46.9059, 46.9028, 46.9040, 
                  46.8949, 46.8973, 46.9030, 46.9021, 46.8974, 46.8974, 46.8977, 46.9017]

# Note: The FCNN architecture (i.e. number of hidden layer nodes) was tuned specifically for K=225 training points.
# When same architecture is used for smaller configurations (e.g., 64 RPs or fewer), the model likely overfits,
# causing fluctuations in accuracy as antenna count increases. 
# These irregular points are excluded for clarity and the results are shown only for selected antenna counts

ant_number_fcnn = [3, 4, 5, 6, 8, 10, 12, 17, 19, 25, 30, 40, 45, 50]
fcnn = [
    18.970834, 17.409723, 16.641902, 16.047993, 15.157581, 14.457430,
    14.118685, 13.381954, 13.320687, 12.898046, 12.840221,
    12.383619, 12.612727, 12.724741
]

# Create the figure and axes
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 8),
                               gridspec_kw={'height_ratios': [0.3, 2.7]})
fig.subplots_adjust(hspace=0.07)  # Adjust vertical space

# Top section (zoomed in for high error values)
ax1.plot(ant_number, centralised_rss, '-^', label='Centralized-RSS', color=color_map['Centralized RSS'])
ax1.plot(ant_number, centralised_hybrid, '-^', label='Centralized-Hybrid', color=color_map['Centralized Hybrid'])
ax1.plot(ant_number, centralised_aoa, '-^', label='Centralized-AOA', color=color_map['Centralized AOA'])
ax1.plot(ant_number, distributed_median, '-o', label='Distributed-Median', color=color_map['Distributed Median'])
ax1.plot(ant_number, distributed_mean, '-o', label='Distributed-Mean', color=color_map['Distributed Mean'])
ax1.plot(ant_number, distributed_bayesian, '-o', label='Distributed-Bayesian', color=color_map['Distributed Bayesian'])
ax1.plot(ant_number, distributed_zscore_tz_1, '-o', label='Distributed-Z-score (Tz=1)', color=color_map['Distributed Zscore (Tz=1)'])
ax1.plot(ant_number, knn_rss, '-^', label='Centralized-RSS (KNN)', color=color_map['Centralized RSS KNN'])
#ax1.plot(ant_number, lr_rss, '-^', label='Centralized-RSS (LR)', color=color_map['Centralized RSS LR'])
ax1.plot(ant_number, distributed_knn, '-o', label='Distributed-Median (KNN)', color=color_map['Distributed Median KNN'])
ax1.plot(ant_number, distributed_lr, '-o', label='Distributed-Median (LR)', color=color_map['Centralized RSS LR'])

# Bottom section (main error range)
ax2.plot(ant_number, centralised_rss, '-^', color=color_map['Centralized RSS'])
ax2.plot(ant_number, centralised_hybrid, '-^', color=color_map['Centralized Hybrid'])
ax2.plot(ant_number, centralised_aoa, '-^', color=color_map['Centralized AOA'])
ax2.plot(ant_number, distributed_median, '-o', color=color_map['Distributed Median'])
ax2.plot(ant_number, distributed_mean, '-o', color=color_map['Distributed Mean'])
ax2.plot(ant_number, distributed_bayesian, '-o', color=color_map['Distributed Bayesian'])
ax2.plot(ant_number, distributed_zscore_tz_1, '-o', color=color_map['Distributed Zscore (Tz=1)'])
ax2.plot(ant_number, knn_rss, '-^', color=color_map['Centralized RSS KNN'])
#ax2.plot(ant_number, lr_rss, '-^', color=color_map['Centralized RSS LR'])
ax2.plot(ant_number, distributed_knn, '-o', color=color_map['Distributed Median KNN'])
ax2.plot(ant_number, distributed_lr, '-o', color=color_map['Centralized RSS LR'])
ax2.plot(ant_number_fcnn, fcnn, '-^', label='FCNN', color=color_map['FCNN'])

# Grids
ax1.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
ax2.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)

# Y-axis limits and ticks
ax1.set_ylim(45, 48)
ax1.set_yticks([45, 47])
ax2.set_ylim(6, 19.5)
ax2.set_yticks(range(6, 19, 2))

# X-axis limits and ticks
ax2.set_xlim(3, 50)
ax2.set_xticks([3] + list(range(5, 55, 5)))

# Tick label sizes
ax1.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis='both', which='major', labelsize=18)

# Hide spines between axes
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Diagonal break lines
d = 0.5  # Diagonal line size
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels
fig.text(0.5, 0.02, 'Number of Antennas per AP (N)', ha='center', fontsize=25)
fig.text(0.06, 0.5, 'Mean Localization Error (m)', va='center', rotation='vertical', fontsize=25)

# Legend
fig.legend(loc='upper right', bbox_to_anchor=(0.89, 0.86), bbox_transform=fig.transFigure,
           fontsize=12, frameon=True, framealpha=0.4)

# Save
plt.savefig('Fig3a_loc_accuracy_vs_ant_L25_SF8dB_Tz1_K64.png', dpi=300, bbox_inches='tight')
plt.show()
