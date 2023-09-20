"""
@author: Jalal Alali
"""

import numpy as np
import matplotlib.pyplot as plt

den_water_sc=1E1 #kg/m3

p_water=np.array([230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,
305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,
400])

den_water=np.array([972.85,973.06,973.28,973.5,973.71,973.93,974.15,974.36,
974.58,974.79,975.01,975.22,975.43,975.65,975.86,976.08,976.29,976.5,976.71,
976.93,977.14,977.35,977.56,977.77,977.99,978.2,978.41,978.62,978.83,979.04,
979.25,979.46,979.67,979.87,980.08])

fvf_water=np.array([1.0279077,1.0276859,1.0274536,1.0272214,1.0269998,
1.0267678,1.0265360,1.0263147,1.0260830,1.0258620,1.0256305,1.0254097,
1.0251889,1.0249577,1.0247372,1.0245062,1.0242858,1.0240655,1.0238454,
1.0236148,1.0233948,1.0231749,1.0229551,1.0227354,1.0225053,1.0222858,
1.0220664,1.0218471,1.0216279,1.0214087,1.0211897,1.0209707,1.0207519,
1.0205435,1.0203249])

#Water Density at STD (1.01325 bar, 15 °C)
den_water_sc=1000.00 #kg/m3

# # Create figure and axes
# fig, ax1 = plt.subplots(1)

# # Plot data for "p" on the first y-axis
# ax1.plot(p_water, den_water, 'r-', label='Density')
# ax1.set_xlabel('Pressure (bar)')
# ax1.set_ylabel('Density (kg/m3)', color='r')
# ax1.tick_params('y', colors='r')

# # # Create a second y-axis
# ax2 = ax1.twinx()

# # # Plot data for "fvf" on the second y-axis
# ax2.plot(p_water, fvf_water, 'b-', label='FVF')

# ax2.set_xlabel('Pressure (bar)')
# ax2.set_ylabel('FVF', color='b')
# ax2.tick_params('y', colors='b')

# # # Add legend
# lines = ax1.get_lines() + ax2.get_lines()
# labels = [line.get_label() for line in lines]
# ax1.set_title('Water Density & FVF Vs. Pressure')
# ax1.legend(lines, labels, loc='best')

# #Adjust the figure size to fit the page
# fig=plt.gcf()
# fig.set_size_inches(8, 6)
# plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)

# # Show the plot
# plt.show()
# plt.savefig('Water FVF Plot.pdf')


