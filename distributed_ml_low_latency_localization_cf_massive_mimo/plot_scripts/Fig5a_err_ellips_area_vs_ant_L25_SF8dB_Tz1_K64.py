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

centralised_rss = [4412.5946, 4411.0039, 4410.4181, 4416.3038, 4406.7156, 4409.8389, 4418.8060,
 4408.8669, 4411.0814, 4406.4410, 4406.3794, 4420.1525, 4431.9171, 4404.8143,
 4404.9821, 4402.4499, 4419.2472]

centralised_aoa = [4442.6138, 3999.2355, 3685.7412, 3490.1946, 3213.4737, 3037.2547, 2933.3028,
 2779.2063, 2739.5775, 2705.4928, 2648.1000, 2584.5708, 2535.0146, 2496.8555,
 2473.4484, 2459.4041, 2436.5686]

centralised_hybrid = [3708.9221, 3356.6409, 3101.9084, 2948.3553, 2730.7472, 2576.0851, 2505.7620,
 2389.1722, 2342.9940, 2314.7086, 2278.3877, 2236.0537, 2199.2267, 2164.1999,
 2136.7081, 2130.0801, 2111.1889]

distributed_median = [8385.6792, 8332.8329, 8331.2954, 8295.2890, 8284.4866, 8251.6624, 8269.7265,
 8275.9974, 8275.1611, 8268.2133, 8291.5824, 8306.0711, 8298.8986, 8297.5639,
 8303.6373, 8282.9555, 8310.3127]

distributed_mean = [456.1576, 444.1115, 437.0775, 431.7266, 425.5173, 421.0532, 418.5842,
 415.0087, 413.9155, 413.1236, 411.9401, 410.2213, 409.1260, 408.9626,
 408.2258, 407.6224, 407.1653]

distributed_bayesian = [
    323.7182, 322.8684, 323.1020, 322.4449, 322.0987, 321.5978, 321.0950, 
    321.5036, 320.9649, 321.1005, 321.1361, 320.8318, 320.7128, 320.5976, 
    320.4848, 320.6629, 320.5629
]

distributed_zscore = [
    412.1532, 413.4865, 415.4650, 416.0347, 417.2826, 417.7987, 418.1377, 
    419.6672, 419.3139, 420.1688, 420.6675, 421.0013, 421.4706, 421.6891, 
    421.8241, 422.3835, 422.5378
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

ax1.set_ylim(7800, 9000)
ax2.set_ylim(0, 5000)
ax2.set_yticks([0,1000, 2000,3000,4000,5000])

ax2.set_xlim(3, 50)
ax2.set_xticks([3] + list(range(5, 55, 5)))

ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

ax1.xaxis.set_tick_params(labelsize=13.5)
ax1.yaxis.set_tick_params(labelsize=12)

ax2.xaxis.set_tick_params(labelsize=13.5)
ax2.yaxis.set_tick_params(labelsize=12)

# Axis breaks
d = 0.5
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Labels and legend
fig.text(0.5, 0.02, 'Number of Antennas per AP (N)', ha='center', fontsize=23)
fig.text(0.04, 0.5, r'Area of 95% Error Ellipse (m$^\boldsymbol{2}$)', va='center', rotation='vertical', fontsize=23)

ax2.legend(fontsize=14, loc='upper right',bbox_to_anchor=(1, 1.06),framealpha=0.4)

# Save and show
plt.savefig('Fig5a_err_ellips_area_vs_ant_L25_SF8dB_Tz1_K64.png', dpi=300, bbox_inches='tight')
plt.show()
