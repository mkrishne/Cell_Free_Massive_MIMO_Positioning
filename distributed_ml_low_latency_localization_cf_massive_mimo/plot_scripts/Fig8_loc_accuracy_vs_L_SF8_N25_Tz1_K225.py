import matplotlib.pyplot as plt

color_map = {
    'Centralised RSS': 'tab:blue',
    'Centralised Hybrid': 'tab:orange',
    'Centralised AOA': 'tab:green',
    'Distributed Median': 'tab:red',
    'Distributed Mean': 'tab:purple',
    'Distributed Bayesian': 'tab:brown',
    'Distributed Zscore': 'tab:pink',
}

# L values
L_values = [5, 10, 15, 20, 25, 30, 35, 40]

centralised_rss = [
    30.174985, 16.501194, 12.157211, 10.202265, 9.006133, 8.175925, 7.668072, 7.210008
]

centralised_hybrid = [
    5.119359, 5.223964, 5.223749, 4.967229, 4.926112, 4.818557, 4.762998, 4.681160
]

centralised_aoa = [
    5.024539, 5.200207, 5.394134, 5.245903, 5.241528, 5.137127, 5.085731, 5.024922
]

distributed_median = [
    10.755688, 7.971310, 7.036707, 6.414104, 6.102863, 5.875072, 5.741433, 5.576133
]

distributed_mean = [
    11.894948, 9.639391, 8.703668, 8.183461, 7.793657, 7.580800, 7.415917, 7.282402
]

distributed_bayesian = [
    11.048180, 9.105919, 8.375196, 8.014217, 7.757946, 7.633195, 7.537549,
    7.442208
]

distributed_zscore = [
    10.675475, 8.284714, 7.423685, 6.980457, 6.701986, 6.560490, 6.440062,
    6.327643
]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(L_values, centralised_rss, marker='^', label='Centralised-RSS',
         color=color_map['Centralised RSS'])

plt.plot(L_values, centralised_hybrid, marker='^', label='Centralised-Hybrid',
         color=color_map['Centralised Hybrid'])

plt.plot(L_values, centralised_aoa, marker='^', label='Centralised-AOA',
         color=color_map['Centralised AOA'])

plt.plot(L_values, distributed_mean, marker='o', label='Distributed-Mean',
         color=color_map['Distributed Mean'])

plt.plot(L_values, distributed_median, marker='o', label='Distributed-Median',
         color=color_map['Distributed Median'])

plt.plot(L_values, distributed_bayesian, marker='o', label='Distributed-Bayesian',
         color=color_map['Distributed Bayesian'])

plt.plot(L_values, distributed_zscore, marker='o', label='Distributed-Z-score',
         color=color_map['Distributed Zscore'])

plt.xlabel('Number of APs (L)',fontsize=25)
plt.ylabel('Mean Localization Error (m)',fontsize=23)
plt.gca().xaxis.set_tick_params(labelsize=17)
plt.gca().yaxis.set_tick_params(labelsize=17)

plt.legend(frameon=True, framealpha=0.5, prop={'size': 17})
plt.grid(True)
plt.savefig('Fig8_loc_accuracy_vs_L_SF8_N25_Tz1_K225.png', dpi=300, bbox_inches='tight')  # Save as PNG
plt.show()
