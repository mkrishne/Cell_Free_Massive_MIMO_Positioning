import matplotlib.pyplot as plt

# Data for the algorithms and their respective coverage percentages
algorithms = [
    'Centr-Hybrid', 'Centr-AOA', 'Centr-RSS',
    'Distr-Median', 'Distr-Mean', 'Distr-Bayesian', 'Distr-Z-score'
]
coverage_percentages = [89.63, 89.31, 97.44, 99.99, 72.51, 66.90, 80.64]

# Create a bar chart
plt.figure(figsize=(10, 6))
plt.barh(algorithms, coverage_percentages, color='skyblue')

# Add labels and title
plt.xlabel('Percentage of Test Points Covered by 95% Error Ellipse (%)', fontsize=22)
plt.ylabel('Localization Algorithms', fontsize=23)
#plt.title('Coverage of 95% Error Ellipse by Different Algorithms', fontsize=16)

# Set different font sizes for x and y tick labels
plt.gca().xaxis.set_tick_params(labelsize=15)  # x-axis tick font size
plt.gca().yaxis.set_tick_params(labelsize=15)  # y-axis tick font size


# Show grid
plt.grid(axis='x', linestyle='--', alpha=0.6)
plt.savefig('Fig6_coverage_err_ellips_L25_SF8dB_Tz1_K225.png', dpi=300, bbox_inches='tight')  # Save as PNG


# Display the plot
plt.show()
