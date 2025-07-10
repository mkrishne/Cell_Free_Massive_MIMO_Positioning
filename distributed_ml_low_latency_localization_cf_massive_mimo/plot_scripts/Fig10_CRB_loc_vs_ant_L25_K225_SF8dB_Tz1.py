import matplotlib.pyplot as plt

# Data extracted from the file
antennas = [
    3, 4, 5, 6, 8, 10, 12, 15, 17, 19, 21, 25, 30, 35, 40,45,50
]
centralised_hybrid = [
    11.1336, 9.5870, 8.6331, 7.9190, 7.0104, 6.3724, 5.9581, 5.5559, 5.3896, 
    5.2403, 5.1470, 4.9550, 4.7592, 4.6735, 4.6069, 4.4900, 4.4727
]

centralised_aoa = [
        12.3965, 10.5733, 9.5338, 8.6364, 7.6176, 6.8985, 6.3907, 5.9577, 5.7590, 
        5.6159, 5.4792, 5.2896, 5.0292, 4.9572, 4.8688, 4.7479, 4.7284
    ]
    
centralised_hybrid_crlb = [
    4.0858, 3.5249, 3.2667, 3.0659, 2.8828, 2.7909, 2.7688, 2.7110, 2.6956, 
    2.6796, 2.7096, 2.6637, 2.6654, 2.6537, 2.6494, 2.6556, 2.6460
]

centralised_aoa_crlb = [
    4.1807, 3.5293, 3.2568, 3.0040, 2.8069, 2.6979, 2.6762, 2.6098, 2.5927, 
    2.5747, 2.5869, 2.5621, 2.5205, 2.5471, 2.5213, 2.5394, 2.5210
]

distributed_median = [
    6.2350, 6.1668, 6.1267, 6.0935, 6.0696, 6.0471, 6.0411, 6.0746, 6.0925, 
    6.1199, 6.1236, 6.1383, 6.1458, 6.1493, 6.1545, 6.1513, 6.1528
]

distributed_median_crlb = [
    6.1719, 6.0994, 6.0548, 6.0282, 5.9928, 5.9743, 5.9573, 5.9554, 5.9546, 
    5.9476, 5.9438, 5.9405, 5.9402, 5.9347, 5.9316, 5.9341, 5.9361
]

distributed_mean = [
        9.0308, 8.7205, 8.5242, 8.3810, 8.2339, 8.1225, 8.0374, 7.9812, 7.9478,
        7.9287, 7.9113, 7.8777, 7.8556, 7.8456, 7.8345, 7.8152, 7.8116
    ]
distributed_mean_crlb = [
        7.8037, 7.6873, 7.6337, 7.6011, 7.5684, 7.5463, 7.5385, 7.5284, 7.5291,
        7.5269, 7.5248, 7.5193, 7.5193, 7.5207, 7.5178, 7.5181, 7.5157
    ]
    
distributed_zscore  = [6.9285, 6.8533, 6.8177, 6.7832, 6.7460, 6.7246, 6.7137, 6.7118, 
                 6.7247, 6.7355, 6.7327, 6.7382, 6.7390, 6.7436, 6.7487, 6.7515, 6.7584]
                 
distributed_zscore_crlb =  [
    6.7625, 6.6951, 6.6685, 6.6454, 6.6221, 6.6117, 6.6024, 6.5966, 6.5966,
    6.5979, 6.5914, 6.5920, 6.5868, 6.5938, 6.5871, 6.5937, 6.5931 ]

color_map = {
    'Centralised Hybrid (MUSIC)': 'tab:orange',
    'Centralised AOA (MUSIC)': 'tab:green',
    'Distributed Median (MUSIC)': 'tab:red',
    'Distributed Mean (MUSIC)': 'tab:purple',
    'Distributed Bayesian (MUSIC)': 'tab:brown',
    'Distributed Zscore (MUSIC)': 'tab:pink',
}
crlb_color_map = {
    'Centralised Hybrid (CRLB)': 'mediumorchid',       # Softer purple-pink
    'Centralised AOA (CRLB)': 'cadetblue',             # Muted teal/blue
    'Distributed Median (CRLB)': 'lightskyblue',       # Gentle, light blue
    'Distributed Mean (CRLB)': 'slategray',            # Neutral, soft gray-blue
    'Distributed Zscore (CRLB)': 'rosybrown',          # Warm muted rose
}

plt.figure(figsize=(10, 6))

# Centralised
plt.plot(antennas, centralised_hybrid, marker='o', label='Centralised-Hybrid (MUSIC)', 
         color=color_map['Centralised Hybrid (MUSIC)'])
plt.plot(antennas, centralised_hybrid_crlb, marker='o', label='Centralised-Hybrid (CRLB)', 
         color=crlb_color_map['Centralised Hybrid (CRLB)'])

plt.plot(antennas, centralised_aoa, marker='o', label='Centralised-AOA (MUSIC)', 
         color=color_map['Centralised AOA (MUSIC)'])
plt.plot(antennas, centralised_aoa_crlb, marker='o', label='Centralised-AOA (CRLB)', 
         color=crlb_color_map['Centralised AOA (CRLB)'])

# Distributed
plt.plot(antennas, distributed_median, marker='^', label='Distributed-Median (MUSIC)', 
         color=color_map['Distributed Median (MUSIC)'])
plt.plot(antennas, distributed_median_crlb, marker='^', label='Distributed-Median (CRLB)', 
         color=crlb_color_map['Distributed Median (CRLB)'])

plt.plot(antennas, distributed_mean, marker='^', label='Distributed-Mean (MUSIC)', 
         color=color_map['Distributed Mean (MUSIC)'])
plt.plot(antennas, distributed_mean_crlb, marker='^', label='Distributed-Mean (CRLB)', 
         color=crlb_color_map['Distributed Mean (CRLB)'])

plt.plot(antennas, distributed_zscore, marker='^', label='Distributed-Z-score (MUSIC)', 
         color=color_map['Distributed Zscore (MUSIC)'])
plt.plot(antennas, distributed_zscore_crlb, marker='^', label='Distributed-Z-score (CRLB)', 
         color=crlb_color_map['Distributed Zscore (CRLB)'])

plt.xlabel('Number of Antennas per AP (N)',fontsize=22)
plt.ylabel('Mean Localization Error (m)',fontsize=22)
plt.xlim(3, 50)
plt.xticks([3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], fontsize=17)
plt.yticks([2,4,6,8,10,12],fontsize=17)
plt.grid(True)
plt.legend(fontsize=9.5)
plt.savefig('Fig10_CRB_loc_vs_ant_L25_K225_SF8dB_Tz1.png', dpi=300, bbox_inches='tight')  # Save as PNG
plt.show()
