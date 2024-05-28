import pandas as pd
import matplotlib.pyplot as plt

def read_data(filename):
    df = pd.read_csv(filename)
    return df['Time'].tolist(), df['Value'].tolist()

times1, values1 = read_data('x1.csv')
times2, values2 = read_data('x2.csv')
times3, values3 = read_data('u.csv')

cp_times1, cp_values1 = read_data('x1_controlpoints.csv')
cp_times2, cp_values2 = read_data('x2_controlpoints.csv')
cp_times3, cp_values3 = read_data('u_controlpoints.csv')

plt.figure(figsize=(10, 6))

line_width = 0.5

plt.plot(times1, values1, marker='.', linestyle=':', color='blue', linewidth=line_width, label='x1')
plt.plot(times2, values2, marker='.', linestyle=':', color='red', linewidth=line_width, label='x2')
plt.plot(times3, values3, marker='.', linestyle=':', color='green', linewidth=line_width, label='u')

plt.scatter(cp_times1, cp_values1, color='white', s=70, label='x Control Points', edgecolor='blue') 
plt.scatter(cp_times2, cp_values2, color='white', s=60, label='x1 Control Points', edgecolor='red')
plt.scatter(cp_times3, cp_values3, color='white', s=50, label='u Control Points', edgecolor='green')  

plt.title('BeBOT Results')
plt.xlabel('Time')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.show()
