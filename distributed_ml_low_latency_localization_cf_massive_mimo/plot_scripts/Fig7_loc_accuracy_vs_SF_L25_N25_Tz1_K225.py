import matplotlib.pyplot as plt
import numpy as np
color_map = {
    'Centralised RSS': 'tab:blue',
    'Centralised Hybrid': 'tab:orange',
    'Centralised AOA': 'tab:green',
    'Distributed Median': 'tab:red',
    'Distributed Mean': 'tab:purple',
    'Distributed Zscore': 'tab:pink',
}

# Arrays
sigma_sf_db  = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14])
import numpy as np

sigma_sf_db = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

centralised_rss = np.array([
    1.957213, 2.268275, 2.682538, 3.214548, 3.928738, 4.787190, 5.924837, 
    7.301879, 9.085506, 11.140955, 13.827751, 16.727453, 20.009599
])

centralised_hybrid = np.array([
    4.468625, 4.478863, 4.518484, 4.574172, 4.655606, 4.702937, 4.785719, 
    4.875013, 4.946961, 5.160180, 5.504218, 6.352655, 8.231510
])

centralised_aoa = np.array([
    5.132362, 5.105897, 5.101275, 5.133147, 5.151586, 5.134354, 5.175420, 
    5.223702, 5.249256, 5.451215, 5.818690, 6.774001, 8.967482
])

distributed_median = np.array([
    1.678945, 1.941745, 2.286986, 2.696421, 3.171613, 3.769398, 4.470429, 
    5.302088, 6.234282, 7.233592, 8.366765, 9.528907, 10.492159
])

distributed_mean = np.array([
    2.502629, 2.851580, 3.314578, 3.875358, 4.511710, 5.268664, 6.070805, 
    6.981081, 7.888605, 8.955533, 10.173478, 11.640289, 13.425115
])

distributed_bayesian = np.array([2.051421, 2.414021, 2.903892, 3.501231, 4.214618, 5.063359, 5.973292, 
    6.943080, 7.907049, 8.816040, 9.738217, 10.730613, 11.815505])

distributed_zscore = np.array([
    1.737249, 2.437944, 2.923422, 3.495353, 4.204584, 5.001194, 5.884809,
    6.489285, 7.132664, 7.826741, 8.381923, 9.148423, 9.948275, 10.831472
])

# Plotting
plt.figure(figsize=(10,6))

# Slicing arrays to get values from 12 down to 0
plt.plot(sigma_sf_db[12::-1], centralised_rss[12::-1],
         label='Centralised-RSS', marker='^', color=color_map['Centralised RSS'])

plt.plot(sigma_sf_db[12::-1], centralised_hybrid[12::-1],
         label='Centralised-Hybrid', marker='^', color=color_map['Centralised Hybrid'])

plt.plot(sigma_sf_db[12::-1], centralised_aoa[12::-1],
         label='Centralised-AOA', marker='^', color=color_map['Centralised AOA'])

plt.plot(sigma_sf_db[12::-1], distributed_median[12::-1],
         label='Distributed-Median', marker='o', color=color_map['Distributed Median'])

plt.plot(sigma_sf_db[12::-1], distributed_mean[12::-1],
         label='Distributed-Mean', marker='o', color=color_map['Distributed Mean'])

plt.plot(sigma_sf_db[12::-1], distributed_zscore[12::-1],
         label='Distributed-Z-score', marker='o', color=color_map['Distributed Zscore'])

plt.xlabel('Shadow Fading (dB)',fontsize=24)
plt.ylabel('Mean Localization Error (m)',fontsize=23)
plt.grid(True)
plt.xticks(np.arange(12, -1, -1))  # Setting x-ticks to display from 12 down to 0
plt.xlim(12, 0)  # Setting the x-axis limits to reverse the direction from 12 to 0
plt.yticks(np.arange(0, 21, 2))
plt.ylim(0, 20)
plt.gca().xaxis.set_tick_params(labelsize=18)
plt.gca().yaxis.set_tick_params(labelsize=18)

plt.legend(loc='upper right', frameon=True, framealpha=0.2, prop={'size': 18})
plt.tight_layout()
plt.savefig('Fig7_loc_accuracy_vs_SF_L25_N25_Tz1_K225.png', dpi=300, bbox_inches='tight')  # Save as PNG
plt.show()
