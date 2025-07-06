import matplotlib.pyplot as plt

# Set the font globally to Arial
plt.rc('font', family='Arial')

# Data from the provided text
ant_number = [3, 4, 5, 6, 8, 10, 12, 15, 17, 19, 21, 25, 30, 35, 40, 45, 50]
# Data from the provided text
rss_gpr = [
    9.0022, 8.9982, 8.9996, 8.9984, 9.0064, 9.0017, 9.0080, 8.9947, 9.0031,
    9.0036, 9.0116, 8.9987, 9.0084, 9.0056, 9.0030, 9.0049, 9.0025
]

aoa_gpr = [
    12.3965, 10.5733, 9.5338, 8.6364, 7.6176, 6.8985, 6.3907, 5.9577, 5.7590,
    5.6159, 5.4792, 5.2896, 5.0292, 4.9572, 4.8688, 4.7479, 4.7284
]

hybrid_gpr = [
    11.1336, 9.5870, 8.6331, 7.9190, 7.0104, 6.3724, 5.9581, 5.5559, 5.3896,
    5.2403, 5.1470, 4.9550, 4.7592, 4.6735, 4.6069, 4.4900, 4.4727
]
rss_lr = [
    19.5375, 19.5382, 19.5336, 19.5322, 19.5300, 19.5294, 19.5252, 19.5227, 19.5173,
    19.5249, 19.5196, 19.5197, 19.5175, 19.5187, 19.5167, 19.5176, 19.5161
]
rss_knn = [
    5.5672, 5.5646, 5.5614, 5.5649, 5.5627, 5.5608, 5.5556, 5.5635, 5.5586,
    5.5518, 5.5550, 5.5557, 5.5495, 5.5530, 5.5495, 5.5470, 5.5440
]
# Create the figure and axes
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 6),
                               gridspec_kw={'height_ratios': [0.75, 3]})
fig.subplots_adjust(hspace=0.06)  # Adjust space between Axes

# Plot each dataset on each subplot
ax1.plot(ant_number, rss_gpr, '-s', label='RSS GPR',markersize=5, markeredgewidth=2)
ax1.plot(ant_number, rss_knn, '-o', label='RSS WKNN')
ax1.plot(ant_number, aoa_gpr, '-^', label='AOA GPR',markersize=5, markeredgewidth=2)
ax1.plot(ant_number, hybrid_gpr, '-o', label='Hybrid GPR')
ax1.plot(ant_number, rss_lr, '-d', label='RSS LR',markersize=5, markeredgewidth=2)

ax2.plot(ant_number, rss_gpr, '-s',markersize=5, markeredgewidth=2)
ax2.plot(ant_number, rss_knn, '-o')
ax2.plot(ant_number, aoa_gpr, '-^',markersize=5, markeredgewidth=2)
ax2.plot(ant_number, hybrid_gpr, '-o')
ax2.plot(ant_number, rss_lr, '-d',markersize=5, markeredgewidth=2)



# Add light grids to both subplots
ax1.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
ax2.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)

# Set y-limits for each axis to create the break
ax1.set_ylim(18, 21)  # Top section
ax1.set_yticks([19,21])
ax2.set_ylim(3, 13)   # Bottom section
#ax2.set_yticks([4,6,8,10,12,13])
ax2.set_yticks([3,5,7,9,11,13])


# Set x-axis limits to start from 3 and go up to 50
ax2.set_xlim(3, 50)
ax2.set_xticks([3] + list(range(5, 55, 5)))

# Increase the font size of the axis numbers (tick labels)
ax1.tick_params(axis='both', which='major', labelsize=17)  # Adjust font size for ax1
ax2.tick_params(axis='both', which='major', labelsize=17)  # Adjust font size for ax2

# Hide spines between broken axes
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Create diagonal lines to indicate the break
d = 0.5  # Diagonal line size
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle='none', color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)  # Slashes for ax1 and ax2
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)

# Add labels and title
fig.text(0.5, 0, 'Number of Antennas per AP (N)', ha='center', fontsize=19)
fig.text(0.05, 0.5, 'Mean Positioning Error (m)', va='center', rotation='vertical', 
         fontsize=19)

# Unified legend at the top right
fig.legend(loc='upper left', bbox_to_anchor=(0.74, 0.78), fontsize=10, frameon=True, framealpha=0.4)

# Save the figure in high quality for an IEEE paper
plt.savefig('Fig2_positioning_across_N.png', dpi=300, bbox_inches='tight')  # Save as PNG
# plt.savefig('high_quality_plot.pdf', dpi=300, bbox_inches='tight')  # Save as PDF

plt.show()
