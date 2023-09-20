"""
@author: Jalal Alali
"""

import numpy as np
import matplotlib.pyplot as plt

p_spann=np.array([150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,
310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490,
500,510,520,530,540,550,560,570,580,590,600])

den=np.array([339.2598818,371.5183818,403.317185,434.0160957,463.1216168,
490.3323815,515.5301064,538.7369052,560.0634069,579.6650226,597.7122246,
614.3729516,629.8032612,644.1431781,657.5157308,670.0278043,681.7718488,
692.8278309,703.265095,713.1439957,722.5172736,731.431196,739.9265001,
748.0391761,755.8011209,763.2406882,770.3831547,777.2511161,783.864826,
790.2424864,796.4004982,802.3536764,808.1154374,813.6979602,819.1123271,
824.368646,829.4761569,834.4433252,839.2779231,843.9871012,848.5774518,
853.0550638,857.4255722,861.694201,865.8658017,869.9448878])

fvf=np.array([0.00551745,0.005038376,0.004641135,0.004312857,0.004041809,
0.003817511,0.003630921,0.003474515,0.00334221,0.003229191,0.00313169,
0.003046764,0.002972118,0.002905952,0.002846851,0.002793689,0.002745566,
0.002701753,0.002661655,0.002624785,0.002590733,0.00255916,0.002529777,
0.002502341,0.002476643,0.002452502,0.002429764,0.002408294,0.002387975,
0.002368702,0.002350387,0.002332948,0.002316314,0.002300423,0.002285217,
0.002270646,0.002256664,0.002243231,0.002230309,0.002217865,0.002205867,
0.002194289,0.002183104,0.00217229,0.002161824,0.002151687])

#CO2 Density at STD (1.01325 bar, 15 °C)
den_co2_sc=1.87184934007942 #kg/m3

# Define the point where we want to interpolate
p_target=350 #bar

# Perform linear interpolation
fvf_interp = np.interp(p_target, p_spann, fvf)
# Print the interpolated values
print("FVF at "+str(p_target)+" bar= ", fvf_interp)


# Perform linear interpolation
den_interp = np.interp(p_target, p_spann, den)
# Print the interpolated values
print("CO2 Density at Standard Condition (SC)= ", den_co2_sc, "kg/m3")
print("CO2 Density at "+str(p_target)+" bar= ", den_interp, "kg/m3")

# # Create figure and axes
# fig, ax1 = plt.subplots(1)

# # Plot data for "p" on the first y-axis
# ax1.plot(p_spann, den, 'r-', label='Density')
# # ax1.scatter(p_target,den_interp, color='black')
# # ax1.annotate("Density at "+str(p_target)+" bar:\n"+str(den_interp), xy=(p_target,den_interp), xytext=(10,20),
# #              textcoords='offset points', ha='left', va='bottom')
# ax1.set_xlabel('Pressure (bar)')
# ax1.set_ylabel('Density (kg/m3)', color='r')
# ax1.tick_params('y', colors='r')

# # # Create a second y-axis
# ax2 = ax1.twinx()

# # # Plot data for "fvf" on the second y-axis
# ax2.plot(p_spann, fvf, 'b-', label='FVF')
# # ax2.scatter(p_target,fvf_interp, color='black')
# # ax2.annotate("FVF at "+str(p_target)+" bar:\n"+str(fvf_interp), xy=(p_target,fvf_interp), xytext=(10,20),
# #                  textcoords='offset points', ha='left', va='bottom')
# ax2.set_xlabel('Pressure (bar)')
# ax2.set_ylabel('FVF', color='b')
# ax2.tick_params('y', colors='b')

# # # Add legend
# lines = ax1.get_lines() + ax2.get_lines()
# labels = [line.get_label() for line in lines]
# ax1.set_title('CO2 Density & FVF Vs. Pressure')
# ax1.legend(lines, labels, loc='best')

# #Adjust the figure size to fit the page
# fig=plt.gcf()
# fig.set_size_inches(8, 6)
# plt.subplots_adjust(left=0.2, right=0.8, bottom=0.2, top=0.8)

# # # Show the plot
# # plt.show()
# plt.savefig('FVF Original.pdf')









