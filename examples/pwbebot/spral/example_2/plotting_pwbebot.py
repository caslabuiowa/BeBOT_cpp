import pandas as pd
import matplotlib.pyplot as plt

def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

# Read main data
times1, values1 = read_data('x.csv')
times2, values2 = read_data('x1.csv')
times3, values3 = read_data('u.csv')

# Read control points data
cp_times1, cp_values1 = read_data('x_controlpoints.csv')
cp_times2, cp_values2 = read_data('x1_controlpoints.csv')
cp_times3, cp_values3 = read_data('u_controlpoints.csv')

# Read continuity points data
cont_times1, cont_values1 = read_data('x_continuity.csv')
cont_times2, cont_values2 = read_data('x1_continuity.csv')
cont_times3, cont_values3 = read_data('u_continuity.csv')

plt.figure(figsize=(10, 6))
line_width = 0.2

# Plot main data
plt.plot(times1, values1, marker='.', linestyle=':', color='blue', linewidth=line_width, label='x')
plt.plot(times2, values2, marker='.', linestyle=':', color='red', linewidth=line_width, label='x1')
plt.plot(times3, values3, marker='.', linestyle=':', color='green', linewidth=line_width, label='u')

# Plot control points
plt.scatter(cp_times1, cp_values1, color='white', s=70, label='x Control Points', edgecolor='blue') 
plt.scatter(cp_times2, cp_values2, color='white', s=60, label='x1 Control Points', edgecolor='red')
plt.scatter(cp_times3, cp_values3, color='white', s=50, label='u Control Points', edgecolor='green')  

# Plot continuity points
plt.scatter(cont_times1, cont_values1, color='white', s=150, label='x Continuity Points', edgecolor='black') 
plt.scatter(cont_times2, cont_values2, color='white', s=150, label='x1 Continuity Points', edgecolor='black')
plt.scatter(cont_times3, cont_values3, color='white', s=150, label='u Continuity Points', edgecolor='black')

plt.title('BeBOT Results')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.show()
