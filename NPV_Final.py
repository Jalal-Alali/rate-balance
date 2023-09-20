"""
@author: Jalal Alali
"""

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import ecl
from ecl.eclfile import EclFile
from ecl.summary import EclSum
from ecl.eclfile import EclKW, EclFileView
from ecl.grid import EclGrid
import os
import sys
import pandas as pd
from CO2_interpolation import den_co2_sc,p_spann,fvf,den
from Water_interpolation import den_water_sc,p_water,den_water,fvf_water

#simulation=input("Simulation File Name:") 
# Source Files  Naming
file_name='NEW_AURORA_BASECASE_CO2_MULT'
original_file=file_name+'.DATA'
well1='EOS-5-7'
well2='EOS-5-C'

#Create a Figure Counter
n=0

#Rate Unit Convert
r_unit_convert=1465094.13  #Sm3/Day
scf_to_sm3=0.028317

#Well Control Variables
BHP_limit=400 #bar
S_girr=0.254

#Openning Data File For Further Precesses
with open(original_file, 'r') as basecase_file:
    file_contents = basecase_file.read()

#Create a Scenario Cluster to Save the Files
os.makedirs('SCENARIO_CLUSTER', exist_ok=True)
os.chdir('SCENARIO_CLUSTER')

#3D Plotting
#Plotting
x_values_list=[]
y_values_list=[]
condition_list=[]

#Initial Amount (Initial Guess)
#Well EOS-5-7
#Initial Rate
r1_MT=0.5 #Mt/Year
for i in range(10):
    r1=r1_MT
    r1_SM3=r1*r_unit_convert
    #Creating First Directory
    r1_folder=well1+'-RATE_'+str(r1)+'_Mt'
    r1_folder=r1_folder.replace('.','P')
    os.makedirs(r1_folder, exist_ok=True)
    #Openning BaseCase DATA File
    os.chdir(r1_folder)
    #Well EOS-5-C
    #Initial Rate
    r2_MT=0.5 #Mt/Year
    for j in range(10):
        r2=r2_MT
        r2_SM3=r2*r_unit_convert
        total_injection=r1+r2
        # All the actions should be take place here
        #Creating Second Directory
        r2_folder=well1+'-RATE_'+str(r1)+'_Mt'+'_AND_'+well2+'-RATE_'+str(r2)+'_Mt'
        r2_folder=r2_folder.replace('.','P')
        os.makedirs(r2_folder, exist_ok=True)
        os.chdir(r2_folder)
        #Copy Original Data To a New-Established File
        new_file_name=file_name+'_R1-'+str(r1)+'_R2-'+str(r2)
        new_file_name=new_file_name.replace('.','P')
        with open(new_file_name+'.DATA', 'w') as new_file:
            new_file.write(file_contents)
        with open(new_file_name+'.DATA', 'r') as new_file:
            new_file=new_file.read()
        # Replace the target string
        new_file = new_file.replace('<rate1>', str(r1_SM3))
        new_file = new_file.replace('<rate2>', str(r2_SM3))
        new_file = new_file.replace('<total_injection>', str(total_injection))
        new_file = new_file.replace('<rate1m>', str(r1))
        new_file = new_file.replace('<rate2m>', str(r2))
     	# Write the file out in the correct directory
        with open(new_file_name+'.DATA', 'w') as new_file_2:
            new_file_2.write(new_file)
        os.system('nohup flow '+new_file_name+'.DATA')
        
        #sys.stdout=open('Simulation Report.txt','a')
        report_file = open("Simulation Report.txt", "w")
        # Redirect the Standard Output to the Report File
        sys.stdout = report_file
        
        # Create a PDF file
        pdf_file = "Plots.pdf"
        # Create a PdfPages object to save the plots to the PDF file
        pdf_pages = PdfPages(pdf_file)
        
        #simulation=input("Simulation File Name:") 
        # Source Files
        hRestartFile=EclFile(new_file_name+'.UNRST')
        UNS=EclFile(new_file_name+'.UNSMRY')
        Smspec=EclSum(new_file_name+'.SMSPEC')
        Sum=EclSum(new_file_name+'.DATA')
        Grid=EclGrid(new_file_name+'.EGRID')
        Init=EclFile(new_file_name+'.INIT')
        iTimeSteps=hRestartFile.num_report_steps()
        
        # Density Data
        print("CO2 Density at SC= ", den_co2_sc, "Kg/Sm3")
        print("Brine Density at SC= ", den_water_sc, "Kg/Sm3")
        
        #Time Steps
        print("Time Steps= ", iTimeSteps, "Total Time Steps= ", np.sum(iTimeSteps))
        #Time Function
        fTime_list=[]
        for iTimeStep in range (0,iTimeSteps):
            fTime=hRestartFile.iget_restart_sim_days(iTimeStep)/365
            fTime_list.append(fTime)
        fTime=np.array(fTime_list)
        print("fTime Emelements= ", len(fTime))
        fTime2=Sum.numpy_vector('TIME')
        fTime2=fTime2/365
        #Pressure Unit
        pressure_unit=Sum.unit('FPR')
        
        #Listing the Accesible Data
        keys_file=hRestartFile.keys()
        # print("Keys from EclFile= ", keys_file)
        keys_sum=Sum.keys(pattern=None)
        # print("Keys from EclSum= ", keys_sum)
        keys_init=Init.keys()
        # print("Keys from INIT= ", keys_init)
        # index_list=Sum.report_index_list()
        # print("Index List= ", index_list)
        num_rate=hRestartFile.num_named_kw('Rate')
        print("Rate Keyword and Their Related Indexes= ", num_rate)
        num_porv=hRestartFile.num_named_kw('PORV')
        print("PORV Keyword and Their Related Indexes= ", num_porv)
        wells=Sum.wells(pattern=None)
        print("Wells= ", wells)

        #Injection and Monitoring Time Steps
        injection=30
        injection_period=hRestartFile.iget_restart_sim_days(injection)
        print("Injection Period= ", injection_period, "days or ", injection_period/365, "years")
        #injection=input("Injection Period: ") can be used instead in order to give the injection period to the file at the beginning
        monitoring=43
        monitoring_period=hRestartFile.iget_restart_sim_days(iTimeSteps-1)
        print("Monitoring Period ", monitoring_period, "days or ", monitoring_period/365, "years")

        #Time Steps Calling
        time_range=Sum.time_range(start=None, end=None, interval="1Y", num_timestep=None, extend_end=True)
        # print("Time Range= ", time_range)
        print("Time Range Elements= ", len(time_range))

        #report_steps=hRestartFile.report_steps()
        #print("Report Steps:", report_steps)
        #report_dates=hRestartFile.report_dates()
        #print("Report Dates:", report_dates)

        #Data Formatted in the File
        #data_formatted=hRestartFile.fprintf_data('NEW_AURORA_CO2_RW6.UNRST',fmt=None)
        #print("Data Formatted in the File= ", data_formatted)

        #File Information
        file_info=Grid.load_from_file(new_file_name+'.EGRID')
        print("File Info=", file_info)
        Active=Grid.get_num_active()
        print("Number of Active Cells= ", Active)
        #Define NX,NY,NZ
        NX=Grid.get_nx()
        print("NX= ", NX)
        NY=Grid.get_ny()
        print("NY= ", NY)
        NZ=Grid.get_nz()
        print("NZ= ", NZ)

        #Index Import
        index=Grid.export_index(active_only=True)
        # print("Indexes Panada Frame= ", index)
        # i=index['i']
        # print("Pandas i Column= ", i)

        indexf=np.array(index)
        # print("Indexes Array= ", indexf)
        print("Length of indexes: ", len(indexf))
        print("Active ACTNUM (Cells) Indexes:")
        #INDEX I
        #i=[115, 144]
        i=indexf[:,0]
        print("i=", i)
        i_min=np.min(i)
        i_max=np.max(i)
        # print('i_min= ', i_min)
        # print('i_max= ', i_max)
        #INDEX J
        j=indexf[:,1]
        print("j=", j)
        j_min=np.min(j)
        j_max=np.max(j)
        # print('j_min= ', j_min)
        # print('j_max= ', j_max)
        #INDEX K
        k=indexf[:,2]
        print("k=", k)
        k_min=np.min(k)
        k_max=np.max(k)
        # print('k_min= ', k_min)
        # print('k_max= ', k_max)

        #Defining The Penalty Area
        #i=[115, 144]
        i_penalty_min=115
        i_penalty_max=144
        print('i_min_penalty= ', i_penalty_min)
        print('i_max_penalty= ', i_penalty_max)
        #j=[20, 80]
        j_penalty_min=20
        j_penalty_max=80
        print('j_min_penalty= ', j_penalty_min)
        print('j_max_penalty= ', j_penalty_max)
        #j=[9, 13]
        k_penalty_min=9
        k_penalty_max=13
        print('k_min_penalty= ', k_penalty_min)
        print('k_max_penalty= ', k_penalty_max)
        print("Penalty Area Lower Limit= ", [i_penalty_min,j_penalty_min,k_penalty_min])
        print("Penalty Area Top Limit= ", [i_penalty_max,j_penalty_max,k_penalty_max])
        # penalty_area_cells=(indexf[:,0]>=(115))&(indexf[:,0]<=144)&(indexf[:,1]>=20)&(indexf[:,1]<=80)
        penalty_area_cells=(i>=(i_penalty_min-1))&(i<=(i_penalty_max-1))&(j>=(j_penalty_min-1))&(j<=(j_penalty_max-1))&((k>=(k_penalty_min-1))&(k<=(k_penalty_max-1)))
        #ijk_penalty=((indexf[:,0])&(indexf[:,1])&(indexf[:,2]))
        #penalty_area=indexf[penalty_area_indexes, indexf]
        #PP=np.create_3d(penalty_area)
        #print("Penalty Area ijk= ", PP)
        cell_penalty=indexf[penalty_area_cells,3]
        print("Cell Numbers of the Penalty Area= ", cell_penalty)
        min_cell_penalty=np.min(cell_penalty)
        max_cell_penalty=np.max(cell_penalty)
        print("Number of active Cells in the Penalty Area= ", len(cell_penalty))

        #FPR Calling
        FPR=Sum.numpy_vector('FPR')
        # print("FPR= ", FPR)
        Max_FPR=np.max(FPR)
        print("Max FPR: ", Max_FPR, pressure_unit)

        #WBHP: EOS-5-7 Calling
        BHP1_total=Sum.numpy_vector('WBHP:'+well1)
        BHP1=[x for x in BHP1_total if x!=0]
        # print("WBHP:EOS-5-7=", BHP1, pressure_unit)
        print("WBHP:"+well1+" Elements= ", len(BHP1_total))
        BHP1_min=np.min(BHP1_total)
        BHP1_max=np.max(BHP1_total)
        # print("Minimum WBHP:EOS-5-7= ", BHP1_min, pressure_unit)
        print("Maximum WBHP:"+well1+"= ", BHP1_max, pressure_unit)

        #WBHP: EOS-5-C Calling
        BHP2_total=Sum.numpy_vector('WBHP:'+well2)
        BHP2=[x for x in BHP2_total if x!=0]
        # print("WBHP:EOS-5-C=", BHP2, pressure_unit)
        print("WBHP:"+well2+" Elements= ", len(BHP2_total))
        BHP2_min=np.min(BHP2_total)
        BHP2_max=np.max(BHP2_total)
        # print("Minimum WBHP:EOS-5-C= ", BHP2_min, pressure_unit)
        print("Maximum WBHP:"+well2+"= ", BHP2_max, pressure_unit)
        
        #Well Control Loop
        print("BHP Limit=", BHP_limit, pressure_unit)
        if BHP1_max>=BHP_limit:
            condition1=0.2
            print("Pressure Build-Up Happened in Well "+well1)
        else:
            condition1=0
            print("No Pressure Build-up Happened in Well "+well1)
            
        if BHP2_max>=BHP_limit:
            condition2=0.3
            print("Pressure Build-Up Happened in Well "+well2)
        else:
            condition2=0
            print("No Pressure Build-up Happened in Well "+well2)
        if Max_FPR>=BHP_limit:
            condition3=0.4
            print("Pressure Build-Up Due to FPR Happened During the Injection Operation")
        else:
            condition3=0
            print("No Pressure Build-Up Due to FPR Happened During the Injection Operation")    
        condition=condition1+condition2+condition3
        condition_list.append(condition)
        #FGIR Calling
        FGIR=Sum.numpy_vector('FGIR')
        FGIR_unit=Sum.unit('FGIR')
        # print("FGIR: ", FGIR)
        # print("FGIR Elements= ", len(FGIR_total))
        FGIR_min=np.min(FGIR)
        FGIR_max=np.max(FGIR)
        print("Minimum FGIR= ", FGIR_min, FGIR_unit)
        print("Maximum FGIR= ", FGIR_max, FGIR_unit)

        if FGIR_min==FGIR_max:
            print("No Rate Decline Happened During the Injection Operation")
        else:
            print("Rate Decline Happened During the Injection Operation")

        #WGIR: EOS-5-7 Calling
        GIR1=Sum.numpy_vector('WGIR:'+well1)
        rate_unit=Sum.unit("FGIR")
        # print("WGIR of EOS-5-7= ", rate1)
        # print("WGIR:"+well1+" Elements= ", len(rate1_total))
        GIR1_min=np.min(GIR1)
        GIR1_max=np.max(GIR1)
        print("Minimum WGIR:"+well1+" = ", GIR1_min, rate_unit)
        print("Maximum WGIR:"+well1+" = ", GIR1_max, rate_unit)

        if GIR1_min==GIR1_max:
            print("No Rate Decline Happened in Well "+well1)
        else:
            print("Rate Decline Happened in Well "+well1)

        #WGIR: EOS-5-C Calling
        GIR2=Sum.numpy_vector('WGIR:'+well2)
        GIR2_min=np.min(GIR2)
        GIR2_max=np.max(GIR2)
        print("Minimum WGIR:"+well2+"= ", GIR2_min, rate_unit)
        print("Maximum WGIR:"+well2+"= ", GIR2_max, rate_unit)

        if GIR2_min==GIR2_max:
            print("No Rate Decline Happened in Well "+well2)
        else:
            print("Rate Decline Happened in Well "+well2)

        #Calling PORV From the ".INIT" File
        PORV=Init.iget_named_kw('PORV',0)
        PORV3D=Grid.create_3d(PORV)
        # PORV3D_nonzero=[x for x in PORV3D if x!=0]
        PORV3D_penalty=PORV3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
        PORV_unit='Sm3'
        print("Minimum PV in the Penalty Area= ", np.min(PORV3D_penalty), PORV_unit)
        print("Maximum PV in the Penalty Area= ", np.max(PORV3D_penalty), PORV_unit)
        print("Total Number of PORV Elements in the Penalty Area: ", len(PORV3D_penalty))
        print("Total Pore Volume of the Penalty Area= ", np.sum(PORV3D_penalty), PORV_unit)

        #FGIT Calling (Total Injected Gas Calculation)
        FGIT=Sum.numpy_vector('FGIT')
        FGIT_mass=FGIT*den_co2_sc
        FGIT_max=np.max(FGIT)
        FGIT_mass_max=np.max(FGIT_mass)
        FGIT_unit=Sum.unit('FGIT')
        print("Total Injected Gas (FGIT)= ", FGIT_max, FGIT_unit," = ", FGIT_mass_max/1E9, "Mt")
        #print("FGIT Vector= ", FGIT)


        #WGIT:EOS-5-7 Calling
        GIT1=Sum.numpy_vector('WGIT:'+well1)
        GIT1_mass=GIT1*den_co2_sc
        GIT1_max=np.max(GIT1)
        GIT1_mass_max=np.max(GIT1_mass)
        print("Total Injected Gas Through Well "+well1+"= ", GIT1_max, FGIT_unit, " = ", GIT1_mass_max/1E9, "Mt")

        #WGIT:EOS-5-C Calling
        GIT2=Sum.numpy_vector('WGIT:'+well2)
        GIT2_mass=GIT2*den_co2_sc
        GIT2_max=np.max(GIT2)
        GIT2_mass_max=np.max(GIT2_mass)
        print("Total Injected Gas Through Well "+well2+"= ", GIT2_max, FGIT_unit, " = ", GIT2_mass_max/1E9, "Mt")

        #GIT & FGIT Volume Plotting
        FGIT_transposed=FGIT.T
        GIT1_transposed=GIT1.T
        GIT2_transposed=GIT2.T
        fig1=plt.figure(n+1)
        plt.plot(fTime2,FGIT,linestyle='-', label='Total Injected CO2 (FGIT)',color='black')
        plt.plot(fTime2,GIT1,linestyle='--', label=well1,color='blue')
        plt.plot(fTime2,GIT2,linestyle='--', label=well2,color='red')
        plt.title('Total Injected CO2 Volume Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('Injected CO2 Volume [Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('FGIT Profile.pdf')
        pdf_pages.savefig(fig1)
        
        #GIT & FGIT Mass Plotting
        fig2=plt.figure(n+2)
        plt.plot(fTime2,FGIT_mass/1E9,linestyle='-', label='Total Injected CO2 Mass',color='black')
        plt.plot(fTime2,GIT1_mass/1E9,linestyle='--', label=well1,color='blue')
        plt.plot(fTime2,GIT2_mass/1E9,linestyle='--', label=well2,color='red')
        plt.title('Total Injected CO2 Mass Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('Injected CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('FGIT Mass Profile.pdf')
        pdf_pages.savefig(fig2)

        # Penalty for 1 Sm3 of Injected CO2 Into the Penalty Area
        CO2_penalty_price=20 #Dollar($/Ton)
        Dollar_to_NOK=10.05
        CO2_penalty_price_NOK=CO2_penalty_price*Dollar_to_NOK
        print("Penalty Cost for each ton of Leaked CO2 to the Penalty Area= ", CO2_penalty_price, "Dollars ($) = ", CO2_penalty_price_NOK, "NOK")

        #Penalty Area Free Gas (CO2) for the Injection Time Step
        SGASi_list=[]
        SGASi_penalty_max_list=[]
        GVi_penalty_list=[]
        GVi_penalty_mass_list=[]
        for iTimeStep in range(0,injection+1):
            SGAS=hRestartFile.iget_named_kw('SGAS',iTimeStep)
            SGAS3D=Grid.create_3d(SGAS)
            SGAS3D_max=np.max(SGAS3D)
            SGASi_list.append(SGAS3D_max)
            SGAS3D_penalty=SGAS3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            SGASi_penalty_max=np.max(SGAS3D_penalty)
            SGASi_penalty_max_list.append(SGASi_penalty_max)
            P=hRestartFile.iget_named_kw('PRESSURE',iTimeStep)
            P3D=Grid.create_3d(P)
            #FVF Interpolation
            fvf_interp_3D = np.interp(P3D, p_spann, fvf)
            fvf_interp_penalty_3D=fvf_interp_3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            GVi_penalty=SGAS3D_penalty*PORV3D_penalty/fvf_interp_penalty_3D
            GVi_penalty_mass=GVi_penalty*den_co2_sc
            GVi_penalty_sum=np.sum(GVi_penalty)
            GVi_penalty_mass_sum=np.sum(GVi_penalty_mass)
            GVi_penalty_list.append(GVi_penalty_sum)
            GVi_penalty_mass_list.append(GVi_penalty_mass_sum)
            #A grapg of SGAS vs. Time is needed
        GVi_penalty_injection=GVi_penalty_list[-1]
        GVi_penalty_mass_injection=GVi_penalty_mass_list[-1]
        # print("GVi Penalty List= ", GVi_penalty_list)
        # print("SGAS in the Penalty Area after the Injection Period= ", SGASi_list)
        print("Accumulated Gas Volume in the Penalty Area After the Injection Period= ", GVi_penalty_injection, FGIT_unit, " = ", GVi_penalty_mass_injection/1E9, "Mt")
        GVi_penalty_max=np.max(GVi_penalty_list)
        GVi_penalty_mass_max=np.max(GVi_penalty_mass_list)
        print("Maximum Accumulated Gas Volume in the Penalty Area During the Injection Period= ", GVi_penalty_max, FGIT_unit, " = ", GVi_penalty_mass_max/1E9, "Mt")


        #Penalty Area Free Gas (CO2) Loop for the Monitoring Time Step
        SGASm_list=[]
        SGASm_penalty_max_list=[]
        GVm_penalty_list=[]
        GVm_penalty_mass_list=[]
        for iTimeStep in range(injection+1,monitoring+1):
            SGAS=hRestartFile.iget_named_kw('SGAS',iTimeStep)
            SGAS3D=Grid.create_3d(SGAS)
            SGAS3D_max=np.max(SGAS3D)
            SGASm_list.append(SGAS3D_max)
            SGAS3D_penalty=SGAS3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            SGASm_penalty_max=np.max(SGAS3D_penalty)
            SGASm_penalty_max_list.append(SGASm_penalty_max)
            P=hRestartFile.iget_named_kw('PRESSURE',iTimeStep)
            P3D=Grid.create_3d(P)
            #FVF Interpolation
            fvf_interp_3D = np.interp(P3D, p_spann, fvf)
            fvf_interp_penalty_3D=fvf_interp_3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            GVm_penalty=SGAS3D_penalty*PORV3D_penalty/fvf_interp_penalty_3D
            GVm_penalty_sum=np.sum(GVm_penalty)
            GVm_penalty_mass=GVm_penalty*den_co2_sc
            GVm_penalty_mass_sum=np.sum(GVm_penalty_mass)
            GVm_penalty_list.append(GVm_penalty_sum)
            GVm_penalty_mass_list.append(GVm_penalty_mass_sum)
            #A grapgh of SGAS vs. Time is needed
        GVm_penalty_monitoring=GVm_penalty_list[-1]
        GVm_penalty_mass_monitoring=GVm_penalty_mass_list[-1]
        # print("GVm Penalty List= ", GVm_penalty_list)
        # print("SGAS in the Penalty Area after the Monitoring Period= ", SGASm_list)
        print("Accumulated Gas Volume in the Penalty Area After the End of the Monitoring Period= ", GVm_penalty_monitoring, FGIT_unit, " = ", GVm_penalty_mass_monitoring/1E9, "Mt")

        GVm_penalty_max=np.max(GVm_penalty_list)
        GVm_penalty_mass_max=np.max(GVm_penalty_mass_list)
        print("Maximum Accumulated Gas Volume in the Penalty Area During the Monitoring Period= ", GVm_penalty_max, FGIT_unit, " = ", GVm_penalty_mass_max/1E9, "Mt")

        #Contribute a Total GV List
        GV_penalty_list=[]
        GV_penalty_mass_list=[]
        GV_penalty_list.append(GVi_penalty_list+GVm_penalty_list)
        GV_penalty_mass_list.append(GVi_penalty_mass_list+GVm_penalty_mass_list)
        GV_penalty=np.array(GV_penalty_list)
        GV_penalty_mass=np.array(GV_penalty_mass_list)
        GV_penalty_max=np.max(GV_penalty)
        GV_penalty_mass_max=np.max(GV_penalty_mass)
        GV_penalty_transposed=GV_penalty.T
        GV_penalty_mass_transposed=GV_penalty_mass.T
        # print("GV Penalty List= ",GV_penalty_list)
        # print("GV Penalty List Arrays= ", GV_penalty_transposed)
        # print("GV Penalty Elements= ", len(GV_penalty_transposed))
        
        #Maximum Gas Volume Status Loop
        if GVi_penalty_max>=GVm_penalty_max:
            GV_max=GVi_penalty_max
            GV_mass_max=GVi_penalty_mass_max
        else:
            GV_max=GVm_penalty_max
            GV_mass_max=GVm_penalty_mass_max
        print("Maximum Accumulated Gas Mass in the Penalty Area During the Hole Project Period= ", GV_mass_max/1E3, 'tons')

        #Penalty Amount Settlement
        Penalty_cost=(GV_mass_max/1E3)*CO2_penalty_price
        Penalty_cost_NOK=Penalty_cost*Dollar_to_NOK
        print("Penalty Cost That the Contractor Company Should Pay to Compensate= ", Penalty_cost/1e+6, "Million Dollars (m$)","= ", Penalty_cost_NOK/1e+6, "Million NOK (mNOK)")
        
        if (Penalty_cost/1e+6)>1:
            print("The Project has no economic justification!")
        else:
            print("The Project is economically justified!")
        
        #Contribute a Total SGAS List
        SGAS_list=[]
        SGAS_list.append(SGASi_list+SGASm_list)
        SGAS=np.array(SGAS_list)
        SGAS_max=np.max(SGAS)
        SGAS_transposed=SGAS.T
        print("Maximum SGAS in the Reservoir= ", SGAS_max)
        SGAS_penalty_list=[]
        SGAS_penalty_list.append(SGASi_penalty_max_list+SGASm_penalty_max_list)
        SGAS_penalty=np.array(SGAS_penalty_list)
        SGAS_penalty_max=np.max(SGAS_penalty)
        SGAS_penalty_transposed=SGAS_penalty.T
        print("Maximum SGAS in the Penalty Area= ", SGAS_penalty_max)

        #SGAS Plotting
        fig3=plt.figure(n+3)
        plt.plot(fTime,SGAS_transposed,linestyle='-',label='SGAS in the Field',color='black')
        plt.plot(fTime,SGAS_penalty_transposed,linestyle='--',label='SGAS in the Penalty Area',color='red')
        plt.title('CO2 Saturation (SGAS) Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('SGAS [-]')
        plt.legend()
        plt.show()
        plt.savefig('SGAS Profile.pdf')
        pdf_pages.savefig(fig3)

        #Rs Calling and Dissolved Gas Calculating
        Rs_max_list=[]
        DISGAS_list=[]
        DISGAS_penalty_list=[]
        SWAT_list=[]
        SWAT_penalty_list=[]
        GV_list=[]
        GV_mass_list=[]
        GV_free_list=[]
        GV_trapped_list=[]
        DISGAS_mass_list=[]
        DISGAS_penalty_mass_list=[]
        GV_free_mass_list=[]
        GV_trapped_mass_list=[]
        GV_free_penalty_mass_sum_list=[]
        for iTimeStep in range(0,iTimeSteps):
            Rs=hRestartFile.iget_named_kw('RS', iTimeStep)
            Rs3D=Grid.create_3d(Rs)
            Rs3D_penalty=Rs3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            Rs3D_max=np.max(Rs3D)
            Rs_max_list.append(Rs3D_max)
            P=hRestartFile.iget_named_kw('PRESSURE',iTimeStep)
            P3D=Grid.create_3d(P)
            #FVF Interpolation
            fvf_interp_3D = np.interp(P3D, p_spann, fvf)
            fvf_interp_penalty_3D=fvf_interp_3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            fvf_water_interp_3D = np.interp(P3D, p_water, fvf_water)
            fvf_water_penalty_interp_3D=fvf_water_interp_3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            SGAS=hRestartFile.iget_named_kw('SGAS',iTimeStep)
            SGAS3D=Grid.create_3d(SGAS)
            SGAS3D_array=np.array(SGAS3D)
            GV=SGAS3D*PORV3D
            GV_mass=GV*den_co2_sc
            GV_sum=np.sum(GV)
            GV_mass_sum=np.sum(GV_mass)
            GV_list.append(GV_sum)
            GV_mass_list.append(GV_mass_sum)
            SGAS3D_free=np.where(SGAS3D_array > S_girr, SGAS3D_array-S_girr, 0)
            SGAS3D_free_penalty=SGAS3D_free[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            GV_free=SGAS3D_free*PORV3D/fvf_interp_3D
            GV_free_mass=GV_free*den_co2_sc
            GV_free_sum=np.sum(GV_free)
            GV_free_mass_sum=np.sum(GV_free_mass)
            GV_free_list.append(GV_free_sum)
            GV_free_mass_list.append(GV_free_mass_sum)
            GV_free_penalty=SGAS3D_free_penalty*PORV3D_penalty/fvf_interp_penalty_3D
            GV_free_penalty_mass=GV_free_penalty*den_co2_sc
            GV_free_penalty_mass_sum=np.sum(GV_free_penalty_mass)
            GV_free_penalty_mass_sum_list.append(GV_free_penalty_mass_sum)
            SGAS3D_trapped=np.where(SGAS3D_array > S_girr, S_girr, SGAS3D_array)
            GV_trapped=SGAS3D_trapped*PORV3D/fvf_interp_3D
            GV_trapped_mass=GV_trapped*den_co2_sc
            GV_trapped_sum=np.sum(GV_trapped)
            GV_trapped_mass_sum=np.sum(GV_trapped_mass)
            GV_trapped_list.append(GV_trapped_sum)
            GV_trapped_mass_list.append(GV_trapped_mass_sum)
            SWAT3D=1-SGAS3D
            SWAT3D_min=np.min(SWAT3D)
            SWAT_list.append(SWAT3D_min)
            SWAT3D_penalty=SWAT3D[(i_penalty_min-1):(i_penalty_max-1),(j_penalty_min-1):(j_penalty_max-1),(k_penalty_min-1):(k_penalty_max-1)]
            SWAT3D_penalty_min=np.min(SWAT3D_penalty)
            SWAT_penalty_list.append(SWAT3D_penalty_min)
            V_water_penalty_rc=SWAT3D_penalty*PORV3D_penalty
            V_water_penalty_sc=V_water_penalty_rc/fvf_water_penalty_interp_3D
            DISGAS_penalty=Rs3D_penalty*V_water_penalty_sc
            DISGAS_penalty_mass=DISGAS_penalty*den_co2_sc
            DISGAS_penalty_sum=np.sum(DISGAS_penalty)
            DISGAS_penalty_mass_sum=np.sum(DISGAS_penalty_mass)
            DISGAS_penalty_list.append(DISGAS_penalty_sum)
            DISGAS_penalty_mass_list.append(DISGAS_penalty_mass_sum)
            V_water_rc=SWAT3D*PORV3D
            V_water_sc=V_water_rc/fvf_water_interp_3D
            DISGAS=Rs3D*V_water_sc
            DISGAS_mass=DISGAS*den_co2_sc
            DISGAS_sum=np.sum(DISGAS)
            DISGAS_mass_sum=np.sum(DISGAS_mass)
            DISGAS_list.append(DISGAS_sum)
            DISGAS_mass_list.append(DISGAS_mass_sum)
        # print("Rs Max List= ", Rs_max_list)
        Rs=np.array(Rs_max_list)
        Rs_max=np.max(Rs)
        Rs_transposed=Rs.T
        print("Maximum Rs= ", Rs_max)

        #Dissolved Gas (CO2) of the Field and the Penalty Area
        DISGAS_penalty=np.array(DISGAS_penalty_list)
        DISGAS_penalty_mass=np.array(DISGAS_penalty_mass_list)
        DISGAS_penalty_max=np.max(DISGAS_penalty)
        DISGAS_penalty_mass_max=np.max(DISGAS_penalty_mass)
        DISGAS_penalty_transposed=DISGAS_penalty.T
        DISGAS_penalty_mass_transposed=DISGAS_penalty_mass.T
        # print("Penalty Area Dissolved Gas= ", DISGAS_penalty_list)

        # Field Dissolved Gas Amount
        DISGAS=np.array(DISGAS_list)
        DISGAS_mass=np.array(DISGAS_mass_list)
        DISGAS_max=np.max(DISGAS)
        DISGAS_mass_max=np.max(DISGAS_mass)
        DISGAS_transposed=DISGAS.T
        DISGAS_mass_transposed=DISGAS_mass.T
        # print("Dissolved Gas in the Reservoir= ", DISGAS_list)
        print("Dissolved Gas in the Penalty Area Till the End of the Monitoring Period= ", DISGAS_penalty[-1], "Sm3", " = ", DISGAS_penalty_mass[-1]/1E9, "Mt")
        print("Dissolved Gas in the Reservoir Till the End of the Monitoring Period= ", DISGAS[-1], "Sm3", " = ", DISGAS_mass[-1]/1E9, "Mt")
        
        #Free Gas in the Penalty Area
        GV_free_penalty_mass=np.array(GV_free_penalty_mass_sum_list)
        # GV_free_penalty_mass_transposed=GV_free_penalty_mass.T
        print("Mobile CO2 Mass in The Penalty Area= ", np.max(GV_free_penalty_mass)/1E9, "Mt")
        
        #Total Leaked Gas
        GV_penalty_area_mass=GV_penalty_mass+DISGAS_penalty_mass
        GV_penalty_area_mass_transposed=GV_penalty_area_mass.T

        #CO2 Volume Leakage into the Penalty Area Plotting
        fig4=plt.figure(n+4)
        plt.plot(fTime,GV_penalty_mass_transposed/1E9,label='CO2 Mass',color='red')
        plt.title('Total CO2 Mass Leaked Into The Penalty Area')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('CO2 Mass in the Penalty Area.pdf')
        pdf_pages.savefig(fig4)    

        #CO2 Mass Leakage into the Penalty Area Plotting
        fig5=plt.figure(n+5)
        plt.plot(fTime,GV_penalty_mass_transposed/1E9,label='CO2 Mass in the Penalty Area',color='red')
        plt.title('Total CO2 Mass Leaked Into The Penalty Area')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('CO2 Mass in the Penalty Area.pdf')
        pdf_pages.savefig(fig5)

        #Rs Plotting
        fig6=plt.figure(n+6)
        plt.plot(fTime,Rs_transposed,linestyle='-', label='Rs',color='blue')
        plt.title('Gas Solution Ratio in the Reservoir (Rs) Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('Rs [Sm3/Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('Rs Profile.pdf')
        pdf_pages.savefig(fig6)

        #Dissolved CO2 Volume Plotting
        fig7=plt.figure(n+7)
        plt.plot(fTime,DISGAS_transposed,linestyle='-',label='Dissolved CO2 Volume in The Reservoir',color='black')
        plt.plot(fTime,DISGAS_penalty_transposed,linestyle='--',label='Dissolved CO2 Volume in the Penalty Area',color='red')
        plt.title('Dissolved CO2 Volume Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('Dissolved CO2 Volume [Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('DISGAS Volume Profile.pdf')
        pdf_pages.savefig(fig7)
        
        #Dissolved CO2 Mass Plotting
        fig8=plt.figure(n+8)
        plt.plot(fTime,DISGAS_mass_transposed/1E9,linestyle='-',label='Dissolved CO2 Mass in The Reservoir',color='black')
        plt.plot(fTime,DISGAS_penalty_mass_transposed/1E9,linestyle='--',label='Dissolved CO2 Mass in the Penalty Area',color='red')
        plt.title('Dissolved CO2 Mass Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('Dissolved CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('DISGAS Mass Profile.pdf')
        pdf_pages.savefig(fig8)
        
        #CO2 Leakage Mass into the Penalty Area & in the Reservoir Plotting
        GV_mass=np.array(GV_mass_list)
        GV_mass_max=np.max(GV_mass)
        GV_mass_transposed=GV_mass.T
        print("Maximum GV in the Reservoir= ", GV_max, "Sm3", " = ", GV_mass_max/1E9, "Mt")
        fig9=plt.figure(n+9)
        plt.plot(fTime, GV_mass_transposed/1E9,linestyle='-', label='CO2 Mass in The Reservoir',color='blue')
        plt.plot(fTime, GV_penalty_mass_transposed/1E9, label='CO2 Mass in the Penalty Area',color='red')
        plt.title('CO2 Mass in The Reservoir and the in The Penalty Area')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('CO2 Mass in the Penalty Area and in the Reservoir.pdf')
        pdf_pages.savefig(fig9)

        #Volume of Free CO2 Captured in the Reservoir (Free CO2 Volume Plotting)
        # print("Free SGAS List= ", SGAS3D_free)
        # print("Free Gas Volume= ", GV_free)
        GV_free=np.array(GV_free_list)
        GV_free_max=np.max(GV_free)
        GV_free_transposed=GV_free.T
        fig10=plt.figure(n+10)
        plt.plot(fTime,GV_free_transposed,linestyle='-', label='Free CO2 Volume',color='blue')
        plt.title('Free CO2 Volume in the Reservoir Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Volume [Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('Free CO2 Volume Profile.pdf')
        pdf_pages.savefig(fig10)
        
        #Mass of Free CO2 Captured in the Reservoir (Free CO2 Mass Plotting)
        # print("Free SGAS List= ", SGAS3D_free)
        # print("Free Gas Volume= ", GV_free)
        GV_free_mass=np.array(GV_free_mass_list)
        GV_free_mass_max=np.max(GV_free_mass)
        GV_free_mass_transposed=GV_free_mass.T
        fig11=plt.figure(n+11)
        plt.plot(fTime,GV_free_mass_transposed/1E9,linestyle='-', label='Free CO2 Mass',color='blue')
        plt.title('Free CO2 Mass in the Reservoir Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('Free CO2 Mass Profile.pdf')
        pdf_pages.savefig(fig11)
        print("Free CO2 After the Monitoring Period=", GV_free[-1], "Sm3", " = ", GV_free_mass[-1]/1E9, "Mt")

        #Volume of the CO2 Captured in the Reservoir By the Trapping Mechanism Plotting (Trapped CO2 Volume Plotting)
        # print("Trapped SGAS List= ", SGAS3D_trapped)
        # print("Trapped Gas Volume= ", GV_trapped)
        GV_trapped=np.array(GV_trapped_list)
        GV_trapped_max=np.max(GV_trapped)
        GV_trapped_transposed=GV_trapped.T
        fig12=plt.figure(n+12)
        plt.plot(fTime,GV_trapped_transposed,linestyle='-', label='Trapped CO2 Volume',color='black')
        plt.title('Trapped CO2 Volume in the Reservoir Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Volume [Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('Trapped CO2 Volume Profile.pdf')
        pdf_pages.savefig(fig12)
        
        #Mass of the CO2 Captured in the Reservoir By the Trapping Mechanism Plotting (Trapped CO2 Mass Plotting)
        # print("Trapped SGAS List= ", SGAS3D_trapped)
        # print("Trapped Gas Volume= ", GV_trapped)
        GV_trapped_mass=np.array(GV_trapped_mass_list)
        GV_trapped_mass_max=np.max(GV_trapped_mass)
        GV_trapped_mass_transposed=GV_trapped_mass.T
        fig13=plt.figure(n+13)
        plt.plot(fTime,GV_trapped_mass_transposed/1E9,linestyle='-', label='Trapped CO2 Mass',color='black')
        plt.title('Trapped CO2 Mass in the Reservoir Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('Trapped CO2 Mass Profile.pdf')
        pdf_pages.savefig(fig13)
        print("Trapped CO2 After the Monitoring Period=", GV_trapped[-1], "Sm3", " = ", GV_trapped_mass[-1]/1E9, "Mt")

        #A comprehensive Plot Comparing Volume of CO2 Captured in the Reservoir By Different Mechanisms
        fig14=plt.figure(n+14)
        plt.plot(fTime2,FGIT,linestyle='-', label='Total CO2 Injected',color='black')
        plt.plot(fTime,GV_free_transposed,linestyle='-', label='Free CO2',color='blue')
        plt.plot(fTime,DISGAS_transposed,linestyle='-',label='Dissolved CO2',color='red')
        plt.plot(fTime,GV_trapped_transposed,linestyle='-', label='Trapped CO2',color='orange')
        plt.title('Captured CO2 Volume By Different Mechanisms Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Volume [Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('Captured CO2 Volume by Different Mechanisms.pdf')
        pdf_pages.savefig(fig14)
        
        #A comprehensive Plot Comparing Mass of CO2 Captured in the Reservoir By Different Mechanisms
        fig15=plt.figure(n+15)
        plt.plot(fTime2,FGIT_mass/1E9,linestyle='-', label='Total CO2 Injected',color='black')
        plt.plot(fTime,GV_free_mass_transposed/1E9,linestyle='-', label='Free CO2',color='blue')
        plt.plot(fTime,DISGAS_mass_transposed/1E9,linestyle='-',label='Dissolved CO2',color='red')
        plt.plot(fTime,GV_trapped_mass_transposed/1E9,linestyle='-', label='Trapped CO2',color='orange')
        plt.title('Captured CO2 Mass By Different Mechanisms Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('Captured CO2 Mass by Different Mechanisms.pdf')
        pdf_pages.savefig(fig15)
        
        #Perform a Compatibilty Test by the Total Injected Volume
        CO2_sum=DISGAS_transposed+GV_free_transposed+GV_trapped_transposed
        fig16=plt.figure(n+16)
        plt.plot(fTime2,FGIT,linestyle='-', label='FGIT',color='red')
        plt.plot(fTime,CO2_sum,linestyle='-', label='Total CO2 Volume Calculated',color='blue')
        plt.title('FGIT Volume By Simulator and Calculations Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Volume [Sm3]')
        plt.legend()
        plt.show()
        plt.savefig('Calculated Injected CO2 Volume Compatibility Test.pdf')
        pdf_pages.savefig(fig16)
        
        #Perform a Compatibilty Test by the Total Injected Mass
        CO2_mass_sum=DISGAS_mass_transposed+GV_free_mass_transposed+GV_trapped_mass_transposed
        fig17=plt.figure(n+17)
        plt.plot(fTime2,FGIT_mass/1E9,linestyle='-', label='FGIT Mass',color='red')
        plt.plot(fTime,CO2_mass_sum/1E9,linestyle='-', label='Total CO2 Mass Calculated',color='blue')
        plt.title('FGIT Mass By Simulator and Calculations Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('CO2 Mass [Mt]')
        plt.legend()
        plt.show()
        plt.savefig('Calculated Injected CO2 Mass Compatibility Test.pdf')
        pdf_pages.savefig(fig17)

        #Contribute a Total SWAT List
        SWAT_penalty=np.array(SWAT_penalty_list)
        SWAT_penalty_min=np.min(SWAT_penalty)
        SWAT_penalty_transposed=SWAT_penalty.T
        SWAT=np.array(SWAT_list)
        SWAT_min=np.min(SWAT)
        SWAT_transposed=SWAT.T
        print("Minimum SWAT in the Penalty Area= ", SWAT_penalty_min)
        print("Minimum SWAT in the Reservoir= ", SWAT_min)

        #SWAT Plotting
        fig18=plt.figure(n+18)
        plt.plot(fTime,SWAT_transposed,linestyle='-',label='SWAT in the Reservoir',color='black')
        plt.plot(fTime,SWAT_penalty_transposed,linestyle='--',label='SWAT in the Penalty Area',color='blue')
        plt.title('Water Saturation Profile Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('SWAT [-]')
        plt.legend()
        plt.show()
        plt.savefig('SWAT Profile.pdf')
        pdf_pages.savefig(fig18)

        # Maximum Pressure Below the Cap Rock Till the End of the Injection Period
        P1_list=[]
        P1_max_list=[]
        for iTimeStep in range(0,injection+1):
            P1=hRestartFile.iget_named_kw('PRESSURE',iTimeStep)
            P13D=Grid.create_3d(P1)
            P13D_top_layer=P13D[:,:,(k_penalty_min-1):(k_penalty_min)]
            P13D_top_layer_max=np.max(P13D_top_layer)
            P1_list.append(P13D_top_layer)
            P1_max_list.append(P13D_top_layer_max)
        Pi=np.max(P1_list)
        print('Time Step=', injection, 'Maximum Pressure Below the Cap Rock Till the End of the', "Injection", "Period= ", Pi, pressure_unit)
        # Maximum Pressure Below the Cap Rock at the End of the Monitoring Period
        P2_list=[]
        P2_max_list=[]
        for iTimeStep in range(injection+1,monitoring+1):
            P2=hRestartFile.iget_named_kw('PRESSURE',iTimeStep)
            P23D=Grid.create_3d(P2)
            P23D_top_layer=P23D[:,:,(k_penalty_min-1):(k_penalty_min)]
            P23D_top_layer_max=np.max(P23D_top_layer)
            P2_list.append(P23D_top_layer)
            P2_max_list.append(P23D_top_layer_max)
        Pm=np.max(P2_list)
        print('Time Step=', monitoring, 'Maximum Pressure Below the Cap Rock Till the End of the', "Monitoring", "Period= ", Pm, pressure_unit)

        #Pressure Status Loop (Pc)
        if Pi>=Pm:
            Pc=Pi
            print("The highest pressure applies to the cap rock during the", "injection", "operation")
        else:
            Pc=Pm
            print("The highest pressure applies to the cap rock after the", "monitoring", "period")

        #Maximum Pressure Definition
        #Choose an array between the two values
        #Maximum Pressure Toward Cap Rock
        Pc_list=[P1_max_list+P2_max_list]
        Pc_top_layer=np.array(Pc_list)
        Pc_top_layer_transposed=Pc_top_layer.T
        # print("Pc List= ", Pc_top_layer_transposed)
        # print("Pc List Element= ", len(Pc_top_layer_transposed))
        print("Maximum Pressure enforced to the Cap Rock= ", Pc, pressure_unit)

        #Cep Rock Fracture Pressure
        CRFP=350
        print("Cap Rock Fracture Pressure= ", CRFP, pressure_unit)

        #Cap Rock Pressure Status Loop
        print("Cap Rock Pressure Status:")
        if Pc>=CRFP:
            print("Maximum Pressure enforced to the Cap Rock is", "higher than the cap rock fracture pressure")
        else:
            print("Maximum Pressure enforced to the Cap Rock is", "less than the cap rock fracture pressure")

        #Cap Rock Pressure Profile Plotting
        fig19=plt.figure(n+19)
        plt.plot(fTime,Pc_top_layer_transposed,label='Maximum pressure enforced to the cap rock',color='blue')
        plt.title('Pressure Enforced to Cap Rock Vs. Time')
        plt.xlabel('Time [years]')
        plt.ylabel('Pressure [BARSA]')
        plt.legend()
        plt.show()
        plt.savefig('Cap Rock Pressure Profile.pdf')
        pdf_pages.savefig(fig19)


        #Reservoir Formation Fracture Pressure
        RFFP=500
        print("Reservoir Formation Fracture Pressure= ", RFFP, pressure_unit)

        #FPR Calling
        FPR=Sum.numpy_vector('FPR')
        FPR_max=np.max(FPR)
        print("Maximum FPR= ", np.max(FPR), pressure_unit)
        BHP1=Sum.numpy_vector('WBHP:'+well1)
        BHP1_max=np.max(BHP1)
        print("Maximum BHP of "+well1+"=", BHP1_max, pressure_unit)
        BHP2=Sum.numpy_vector('WBHP:'+well2)
        BHP2_max=np.max(BHP2)
        print("Maximum BHP of "+well2+"=", BHP2_max, pressure_unit)
        # FGIR=Smspec.numpy_vector('FGIR')
        # print("FGIR: ", FGIR)

        #Formation Pressure Status Loop
        print("Reservoir Formation Pressure Status:")
        if BHP1_max>=RFFP:
            print("Injection pressure in well "+well1+" is higher than the formation fracture presure")
        else:
            print("Injection pressure in well "+well1+" is less than the formation fracture presure")

        if BHP2_max>=RFFP:
            print("Injection pressure in well "+well2+" is higher than the formation fracture presure")
        else:
            print("Injection pressure in well "+well2+" is less than the formation fracture presure")

        if FPR_max>=RFFP:
            print("Field Average Pressure is higher than the formation fracture presure")
        else:
            print("Field Average Pressure is less than the formation fracture presure")

        # Close the Report File
        report_file.close()
        
        #3Dplotting
        x_values_list.append(r1_MT)
        y_values_list.append(r2_MT)
        condition_list.append(condition)
        
        # Close the Pdf Pages object
        n=n+20
        pdf_pages.close()
        r2_MT=r2_MT+(0.5)
        os.chdir("..")
    r1_MT=r1_MT+(0.5)
    os.chdir("..")
#3D Plotting
x_values=np.array(x_values_list)
y_values=np.array(y_values_list)
condition=np.array(condition_list)
fig=plt.figure(n+1)
color_mapping = {
    0.0: 'blue',
    0.2: 'yellow',
    0.3: 'orange',
    0.4: 'pink',
    0.5: 'magenta',
    0.6: 'purple',
    0.7: 'maroon',
    0.9: 'red'}

# Create a new array by replacing values with corresponding colors
new_condition = np.array([color_mapping[val] if val in color_mapping else val for val in condition])
# print("Colors= ", new_condition)

# Iterate over the points and plot them
for k in range(len(x_values)):
    x = x_values[k]
    y = y_values[k]
    color = new_condition[k]

    # Plot the point
    plt.scatter(x, y, s=40, c=color)

    # Annotate the point with x,y values
    # plt.annotate(f'({x},{y})', (x, y), textcoords="offset points", xytext=(5, 5), ha='center',fontsize=7)

#  Set axis range limits
plt.xlim(x_values[0]-0.5, x_values[-1]+0.5)  # X-axis range: 2 to 5
plt.ylim(y_values[0]-0.5, y_values[-1]+0.5)  # Y-axis range: 2 to 5
# Set labels for the axes
plt.title('Wells Injection Rates vs. Build-up/No-Build-up')
plt.xlabel('EOS-5-7 Rate (Mt/Year)')
plt.ylabel('EOS-5-C Rate (Mt/Year)')
plt.legend()
# Show the plot
plt.show()
plt.savefig('Rates Balance 2D Plot.pdf')
os.chdir("..")
basecase_file.close()

#Finished




