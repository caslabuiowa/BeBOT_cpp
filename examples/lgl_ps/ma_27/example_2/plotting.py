import pandas as pd
import matplotlib.pyplot as plt

# Function to read data from a CSV file
def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

# Read data from the three CSV files
times1, values1 = read_data('x.csv')
times2, values2 = read_data('x1.csv')
times3, values3 = read_data('u.csv')

# Plotting
plt.figure(figsize=(10, 6))

# Plotting the first trajectory
plt.plot(times1, values1, marker='o', linestyle='-', color='blue', label='x')

# Plotting the second trajectory
plt.plot(times2, values2, marker='x', linestyle='-', color='red', label='x1')

# Plotting the third trajectory
plt.plot(times3, values3, marker='^', linestyle='-', color='green', label='u')

plt.title('Lagrange Results Comparison')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.show()
