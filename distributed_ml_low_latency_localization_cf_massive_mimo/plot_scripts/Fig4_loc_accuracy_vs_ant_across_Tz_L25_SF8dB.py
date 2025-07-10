import matplotlib.pyplot as plt

antennas = [3, 4, 5, 6, 8, 10, 12, 15, 17, 19, 21, 25, 30, 35, 40, 45, 50]
# Extracted values for tz_0p6, tz_1 and tz_1p4

# Given the values for each of the required variables, we can store them in Python lists

tz_0p6 = [7.0844, 7.0019, 6.9619, 6.9296, 6.8864, 6.8715, 6.8617, 6.8688, 
                  6.8854, 6.9036, 6.9018, 6.9201, 6.9107, 6.9185, 6.9206, 6.9214, 6.9279]

tz_1 = [6.9285, 6.8533, 6.8177, 6.7832, 6.7460, 6.7246, 6.7137, 6.7118, 
                 6.7247, 6.7355, 6.7327, 6.7382, 6.7390, 6.7436, 6.7487, 6.7515, 6.7584]

tz_1p4 = [7.3104, 7.2436, 7.2113, 7.1810, 7.1506, 7.1286, 7.1182, 7.1083, 
                  7.1164, 7.1176, 7.1179, 7.1192, 7.1209, 7.1257, 7.1256, 7.1297, 7.1324]


tz_0p6_k_64 = [7.3565, 7.2698, 7.2264, 7.1929, 7.1467, 7.1354, 7.1209, 7.1230, 7.1256, 7.1438, 7.1532, 7.1587, 7.1553, 7.1414, 7.1527, 7.1577, 7.1505]


tz_1_k_64 = [7.5074, 7.4449, 7.4231, 7.3932, 7.3630, 7.3461, 7.3316, 7.3357, 7.3333, 7.3527, 7.3581, 7.3596, 7.3441, 7.3464, 7.3597, 7.3605, 7.3622]


tz_1p4_k_64 = [7.9145, 7.8587, 7.8361, 7.8182, 7.7883, 7.7761, 7.7637, 7.7659, 7.7619, 7.7728, 7.7794, 7.7860, 7.7690, 7.7685, 7.7868, 7.7879, 7.7810]


# Data is ready for further processing or visualization!
# Plotting all `tz` variables

tz_0p6_label = r'$\mathcal{T}_z = 0.6\ (K=225)$'
tz_1_label = r'$\mathcal{T}_z = 1$ (K=225)'
tz_1p4_label = r'$\mathcal{T}_z = 1.4$ (K=225)'

tz_0p6_label_k_64 = r'$\mathcal{T}_z = 0.6\ (K=64)$'
tz_1_label_k_64 = r'$\mathcal{T}_z = 1$ (K=64)'
tz_1p4_label_k_64 = r'$\mathcal{T}_z = 1.4$ (K=64)'

# Plotting with corresponding labels
plt.figure(figsize=(10, 6))
plt.plot(antennas, tz_0p6, marker='d', label=tz_0p6_label)
plt.plot(antennas, tz_1, marker='d', label=tz_1_label)
plt.plot(antennas, tz_1p4, marker='d', label=tz_1p4_label)

plt.plot(antennas, tz_0p6_k_64, marker='o', label=tz_0p6_label_k_64)
plt.plot(antennas, tz_1_k_64, marker='o', label=tz_1_label_k_64)
plt.plot(antennas, tz_1p4_k_64, marker='o', label=tz_1p4_label_k_64)

# Adding axis labels and title
plt.xlabel('Number of Antennas per AP (N)', fontsize=22)
plt.ylabel('Mean Localization Error (m)', fontsize=22)

# Setting axis limits and ticks
plt.xlim(3, 50)
plt.xticks([3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], fontsize=16)
plt.yticks(fontsize=17)
plt.ylim(6.5, 8)

# Adding grid and legend
plt.grid(True)
plt.legend(fontsize=15, loc='upper right', framealpha=0.3, bbox_to_anchor=(1, 0.87))

# Saving the plot
plt.savefig('Fig4_loc_accuracy_vs_ant_across_Tz_L25_SF8dB.png', dpi=300, bbox_inches='tight')  # Save as PNG
plt.show()
