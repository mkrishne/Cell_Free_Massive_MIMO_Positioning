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
    9.0022, 8.9982, 8.9996, 8.9984, 9.0064, 9.0017, 9.0080, 8.9947, 9.0031,
    9.0036, 9.0116, 8.9987, 9.0084, 9.0056, 9.0030, 9.0049, 9.0025
]

centralised_hybrid = [
    11.1336, 9.5870, 8.6331, 7.9190, 7.0104, 6.3724, 5.9581, 5.5559, 5.3896,
    5.2403, 5.1470, 4.9550, 4.7592, 4.6735, 4.6069, 4.4900, 4.4727
]

centralised_aoa = [
    12.3965, 10.5733, 9.5338, 8.6364, 7.6176, 6.8985, 6.3907, 5.9577, 5.7590,
    5.6159, 5.4792, 5.2896, 5.0292, 4.9572, 4.8688, 4.7479, 4.7284
]

distributed_median = [
    6.2350, 6.1668, 6.1267, 6.0935, 6.0696, 6.0471, 6.0411, 6.0746, 6.0925,
    6.1199, 6.1236, 6.1383, 6.1458, 6.1493, 6.1545, 6.1513, 6.1528
]

distributed_mean = [
    9.0308, 8.7205, 8.5242, 8.3810, 8.2339, 8.1225, 8.0374, 7.9812, 7.9478,
    7.9287, 7.9113, 7.8777, 7.8556, 7.8456, 7.8345, 7.8152, 7.8116
]

distributed_bayesian = [
    8.1366, 8.0421, 7.9905, 7.9486, 7.9036, 7.8678, 7.8518, 7.8358, 
    7.8345, 7.8344, 7.8316, 7.8281, 7.8229, 7.8234, 7.8261, 7.8289, 
    7.8249
]

distributed_zscore_tz_1 = [6.9285, 6.8533, 6.8177, 6.7832, 6.7460, 6.7246, 6.7137, 6.7118, 
                 6.7247, 6.7355, 6.7327, 6.7382, 6.7390, 6.7436, 6.7487, 6.7515, 6.7584]

centralized_hybrid_no_noise = [
    16.6990, 14.5018, 12.7739, 11.6756, 10.2871, 9.3306, 8.6166, 7.9255,
    7.6293, 7.4804, 7.2240, 6.9292, 6.6524, 6.4118, 6.3308, 6.1294, 6.0693
]

knn_rss = [
    5.5672, 5.5646, 5.5614, 5.5649, 5.5627, 5.5608, 5.5556, 5.5635, 5.5586,
    5.5518, 5.5550, 5.5557, 5.5495, 5.5530, 5.5495, 5.5470, 5.5440
]

lr_rss = [
    19.5375, 19.5382, 19.5336, 19.5322, 19.5300, 19.5294, 19.5252, 19.5227, 19.5173,
    19.5249, 19.5196, 19.5197, 19.5175, 19.5187, 19.5167, 19.5176, 19.5161
]

distributed_median_no_noise = [
    17.4812, 17.4757, 17.2091, 16.8359, 17.0391, 16.9105, 17.1079, 16.7953,
    16.9302, 17.1616, 17.1397, 17.0620, 16.9441, 17.0941, 17.0888, 17.1959, 17.3121
]

distributed_bayesian_no_noise = [
    19.2617, 19.2841, 19.2675, 19.2041, 19.2468, 19.1732, 19.2426, 19.2844, 
    19.3029, 19.3043, 19.3539, 19.3624, 19.3869, 19.4750, 19.4592, 19.5656, 
    19.5473
]

distributed_bayesian_zscore_no_noise = [14.4602, 14.3866, 14.3478, 14.3508, 14.3243, 14.3077, 14.3165, 14.3374, 
                   14.3584, 14.4162, 14.4060, 14.4667, 14.4862, 14.5443, 14.6075, 14.6031, 14.5338]


distributed_knn = [5.1439, 5.0611, 5.0193, 4.9999, 4.9742, 4.9703, 4.9726, 5.0019, 5.0298, 
                    5.0558, 5.0707, 5.0910, 5.1108, 5.1220, 5.1254, 5.1293, 5.1455]

fcnn =     [12.972596, 11.395674, 10.183187, 9.421949, 8.345975, 7.762053, 7.368658,
    6.863254, 6.579629, 6.517370, 6.414308, 6.253715, 6.072018, 5.901977, 5.841282,
    5.711648, 5.798684]

# Create the figure and axes with new break
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 8),
                               gridspec_kw={'height_ratios': [0.3, 2.7]})
fig.subplots_adjust(hspace=0.07)

# Top axis (values above 19)
ax1.plot(ant_number, centralised_rss, '-^', label='Centralized-RSS', color=color_map['Centralized RSS'])
ax1.plot(ant_number, centralised_hybrid, '-^', label='Centralized-Hybrid', color=color_map['Centralized Hybrid'])
ax1.plot(ant_number, centralised_aoa, '-^', label='Centralized-AOA', color=color_map['Centralized AOA'])
ax1.plot(ant_number, distributed_median, '-o', label='Distributed-Median', color=color_map['Distributed Median'])
ax1.plot(ant_number, distributed_mean, '-o', label='Distributed-Mean', color=color_map['Distributed Mean'])
ax1.plot(ant_number, distributed_bayesian, '-o', label='Distributed-Bayesian', color=color_map['Distributed Bayesian'])
ax1.plot(ant_number, distributed_zscore_tz_1, '-o', label='Distributed-Z-score (Tz=1)', color=color_map['Distributed Zscore (Tz=1)'])
ax1.plot(ant_number, knn_rss, '-^', label='Centralized-RSS (KNN)', color=color_map['Centralized RSS KNN'])
ax1.plot(ant_number, distributed_knn, '-o', label='Distributed-Median (KNN)', color=color_map['Distributed Median KNN'])
ax1.plot(ant_number, lr_rss, '-^', label='Centralized-RSS (LR)', color=color_map['Centralized RSS LR'])
ax1.plot(ant_number, distributed_bayesian_no_noise, '-o', label='Distr-Bayesian (no train noise)', color=color_map['Distr Bayesian (no train noise)'])
ax1.plot(ant_number, fcnn, '-^', label='FCNN', color=color_map['FCNN'])

# Bottom axis (values below 14)
ax2.plot(ant_number, centralised_rss, '-^', color=color_map['Centralized RSS'])
ax2.plot(ant_number, centralised_hybrid, '-^', color=color_map['Centralized Hybrid'])
ax2.plot(ant_number, centralised_aoa, '-^', color=color_map['Centralized AOA'])
ax2.plot(ant_number, distributed_median, '-o', color=color_map['Distributed Median'])
ax2.plot(ant_number, distributed_mean, '-o', color=color_map['Distributed Mean'])
ax2.plot(ant_number, distributed_bayesian, '-o', color=color_map['Distributed Bayesian'])
ax2.plot(ant_number, distributed_zscore_tz_1, '-o', color=color_map['Distributed Zscore (Tz=1)'])
ax2.plot(ant_number, knn_rss, '-^', color=color_map['Centralized RSS KNN'])
ax2.plot(ant_number, distributed_knn, '-o', color=color_map['Distributed Median KNN'])
ax2.plot(ant_number, lr_rss, '-^', color=color_map['Centralized RSS LR'])
ax2.plot(ant_number, distributed_bayesian_no_noise, '-o', color=color_map['Distr Bayesian (no train noise)'])
ax2.plot(ant_number, fcnn, '-^', color=color_map['FCNN'])

# Grids
ax1.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
ax2.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)

# Set Y-limits and ticks
ax1.set_ylim(18, 20)  # Set this according to your actual data
ax1.set_yticks([18, 20])  # Adjust ticks to suit your data
ax2.set_ylim(4, 13.5)
ax2.set_yticks([4,6,8,10,12,13])

# X-axis
ax2.set_xlim(3, 50)
ax2.set_xticks([3] + list(range(5, 55, 5)))

# Tick labels
ax1.tick_params(axis='both', which='major', labelsize=18)
ax2.tick_params(axis='both', which='major', labelsize=18)

# Hide spines between axes
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Diagonal lines
d = 0.5
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels
fig.text(0.5, 0.02, 'Number of Antennas per AP (N)', ha='center', fontsize=25)
fig.text(0.06, 0.5, 'Mean Localization Error (m)', va='center', rotation='vertical', fontsize=25)

# Legend
fig.legend(loc='upper right', bbox_to_anchor=(0.89, 0.86), bbox_transform=fig.transFigure,
           fontsize=15, frameon=True, framealpha=0.4)

# Save
plt.savefig('Fig3b_loc_accuracy_vs_ant_L25_SF8dB_Tz1_K225.png', dpi=300, bbox_inches='tight')
plt.show()

