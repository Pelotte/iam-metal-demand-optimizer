'''
This code aims to calculate metal demand linked to Integrated Assessment Models (IAM) energy projections,
in addition to the metal demand for the rest of the transition and broader economy. Then, it generates
an optimisation of the IAM power mix, minimizing the variation in the energy sources used to meet IAM
energy demand projections, without exceeding metal constraints (reserve, resources, mining production).

Author : Pénélope Bieuville
        Master student at CIRAIG
Supervisor : Guillaume Majeau-Bettez
Co-supervisor : Anne De Bortoli
Created at Montréal, Canada, in 2024
'''


# Import usual libraries
import pandas as pd
import numpy as np
import os
from pyomo.environ import *
import pyomo as po
# Import ipopt solver
import cyipopt as ipopt

# Choose the solver for the optimization
#opt = SolverFactory('ipopt', executable="C:/Users/Penel/Ipopt-3.11.1-win64-intel13.1/bin/ipopt.exe")
#opt.options['max_iter'] = 10000  # define maximum number of optimization iterations

opt = SolverFactory('cplex')
opt.options['qpmethod'] = 2
# Define mipgap, as relative difference to optimal solution
opt.options['mip_tolerances_mipgap'] = 0.05  # Écart relatif de 5 %
opt.options['optimality'] = 1
# Adjust convergence parameters
opt.options['barrier convergetol'] = 1e-2  # Set barrier feasibility tolerance

class IAM_Metal_Optimisation_EV :
    '''
    Init folder :
    If the Modelisation Type is Init, this code calculates the cumulated and annual demand of 29 metals
    necessary to meet initial IAM power projections until 2050. The demand is divided in three sectors :
    Energy sector demand, Other transition demand and Other sector demand.
    Opti folder :
    If the Modelisation Type is Opti, the code optimise the power mix, by calculating the minimum variation necessary
    to meet metal constraints (reserve and production) without increasing fossil fuel technologies. It also calculates
    cumulated and annual demand of the new mix.
    '''

    def __init__(self, folder_path, result_path, model_s0, scenario, ModelisationType, Alpha, ResLimit):
        '''
        This function aims to initialize the opti code variables and parameters.

        :param folder_path: [string] path to external datas used in the code, finish by /
        :param result_path: [string] path to stock results, finish by /
        :param model_s0: [string] Integrated Assessment Model (IAM) that wants to be studied
        (AIM..CGE, GCAM4, IMAGE, MESSAGE-GLOBIOM, REMIND-MAGPIE, WITCH-GLOBIOM)
        :param scenario: [string] scenario of Shared Socioeconomic Pathway and  Representative Concentration Pathway
        (SSP-RCP) that wants to be studied
        :param ModelisationType: [string] 'Init' to study initial IAM power mix metal demand
                                 'Opti' to generate an optimisation of IAM power mix
        :param Alpha: [int] percentage value of acceptable share of metal energy demand in comparison to OSD and Storage_Demand
        for the optimisation metal production by year constraint (See CstrMetalMining) ; original value is 5
        :param ResLimit : [string] Limit for the optimisation of cumulated metal demand : "Reserves" or "Resources" ;
        original value is "Resources"
        '''
        self.folder_path = folder_path # Folder path for external datas used
        self.result_path = result_path # Folder path to save results
        self.model_s0 = model_s0 # Chosen IAM to study
        self.scenario = scenario # Chosen scenario ssp-rcp to study
        self.ResLimit = ResLimit # Chosen limit for cumulated metal demand in 2050 ("Reserves" or "Resources")
        # Choose to study the demand in metal of the initial IAM scenario ("Init") or with an optimisation ("Opti")
        self.ModelisationType = ModelisationType
        # Acceptable share of metal energy demand in comparison to OSD and Storage_Demand for the optimisation
        self.Alpha = Alpha

        # Lifetime of a vehicle
        self.Lifetime_V = 12

        # Penalisation by M of the relaxation variable in the objective function
        self.M = 10**6
        # Folder path for results according to the modelisation type chosen
        self.Res_folder = self.result_path+ self.ModelisationType

        # Definition of the scope of this study in terms of models, scenarios, regions, decades, metals and techno
        # External excel defining the scope of the study
        studyScope = self.folder_path + 'Scope of the study.xlsx'

        with pd.ExcelFile(studyScope) as file:

            # Power capacities

            # list of the IAM models you can optimise
            self.listModels = pd.read_excel(file, 'model', index_col=0).squeeze().tolist()
            # list of the SSP - RCP scenarios you can optimise (SSP3 excluded)
            self.listScenarios = pd.read_excel(file, 'scenario', index_col=0).squeeze().tolist()
            # List of the 5 regions studied by the IAM
            self.listRegions = pd.read_excel(file, 'region', index_col=0).squeeze().tolist()
            # List of the decades studied
            self.listDecades = pd.read_excel(file, 'decades', index_col=0).squeeze().astype(str).tolist()
            # List of the metals studied
            self.listMetals = pd.read_excel(file, 'metals', index_col=0).squeeze().tolist()
            # List of the sub-technologies studied
            self.listTechno = pd.read_excel(file, 'technologies precises', index_col=0).squeeze().tolist()
            # List of the initial disaggregated technologies studied
            self.listTechnoAgg= pd.read_excel(file, 'technologies flex', index_col=0).squeeze().tolist()

            # Vehicles

            # List of the aggregated vehicles studied
            self.listVehicleType = ['ICEV', 'PHEV', 'BEV']
            # List of the disaggregated batteries studied
            self.listBatteryDisag = pd.read_excel(file, 'EV battery disag', index_col=0).squeeze().tolist()
            # List of the aggregated batteries studied
            self.listBatteryAgg = pd.read_excel(file, 'EV battery ag', index_col=0).squeeze().tolist()
            # List of the motors studied
            self.listMotor = pd.read_excel(file, 'EV motor', index_col=0).squeeze().tolist()

        # Create a list of Years from the first decade to the last one
        listYears = [str(year) for year in range(int(self.listDecades[0]), int(self.listDecades[-1]) + 1)]
        self.listYears = listYears

        self.list_i = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        # List of the decades with an additional past one (ex : from 2020, add 2010)
        listDecadeswPast = self.listDecades.copy()
        listDecadeswPast.insert(0, str(int(self.listDecades[0]) - 10))
        self.listDecadeswPast = listDecadeswPast

        # Create a list to store fossil technologies
        self.listTechno_Foss = []
        # Iterate through the list of technologies and append the indices of those containing "Foss"
        for techno in self.listTechno:
            if "Foss" in techno:
                self.listTechno_Foss.append(techno)

        # Create a list to store renewable and carbon neutral technologies
        self.listTechno_Ren = []
        # Iterate through the list of technologies and append the names of those not containing "Foss"
        for techno in self.listTechno:
            if "Foss" not in techno:
                self.listTechno_Ren.append(techno)

        # Import general datas on metal : resources, reserves, 2020 production, recovery rates
        self.MetalData = pd.read_excel(self.folder_path + 'Metal data.xlsx', index_col=0)

        # Creation of a list of metals with known ResLimit (exclusion of metal with unknown reserve or resources)
        self.listMetals_knownRes = self.listMetals.copy() # Create a copy of all the metals studied
        for metal in self.listMetals:
            if pd.isna(self.MetalData[self.ResLimit].loc[metal]): # If the data for ResLimit is missing
                self.listMetals_knownRes.remove(metal)  # Remove the metal

        # Import datas for total 2050 biomass availability [GW], according to various scenario
        self.Biomass_availability = pd.read_excel(self.folder_path + 'Future techno availability.xlsx', sheet_name = 'Biomass' , index_col=0)

        # Import datas for remaining hydroelectricity availability, for various two economic and ecologic scenarios, for each regions [GW]
        self.Hydro_availability = pd.read_excel(self.folder_path + 'Future techno availability.xlsx', sheet_name = 'Hydro',index_col=0)

        # Import CF from the literature
        self.CF_from_litt = pd.read_excel(self.folder_path + 'CF_from_litt.xlsx', index_col=0)

        # Import matrix for corresponding CF of aggregated and disaggregated sub-technologies
        Matrix_CF_Init_disag = pd.read_excel(self.folder_path + 'Matrix_CF_disagg.xlsx',
                                             index_col=0)
        self.Matrix_CF_disag = Matrix_CF_Init_disag.sort_values(by='Techno ') # Rearrange matrix in the alphabetic order of technologies

        # Import metal consumption by decade for storage in [t/decade] (from IEA)
        self.Storage_Demand_init = pd.read_excel(self.folder_path + 'Demand_Storage&Network_IEA.xlsx', index_col=0)

        # Import initial metal consumption by decade for network in [t/decade] (from IEA)
        self.Network_Demand_init = pd.read_excel(self.folder_path + 'Demand_Storage&Network_IEA.xlsx',
                                            sheet_name='Demand_Network', index_col=0)

        # Tab of future metal production according to different scenarios from 2020 to 2050 [t/yr]
        #self.Prod = pd.read_excel(self.folder_path + 'Future metal prod.xlsx', index_col=0)

        # Importation of GDP projections for a specific model, scenario ssp-rcp, at a specific year (from IAM)
        self.GDP_folder = self.folder_path + 'GDP IAM/GDP ' + self.model_s0 + '_' + self.scenario

        # Importation of GDP datas for a specific model, ssp, at a specific year
        self.Population_folder = self.folder_path + 'Population IAM/Population ' + model_s0 + '_' + self.scenario

        # Folder used to stock initial Power capacity projections for specific model_s0 and scenario, for various regions
        self.folderIAM = self.folder_path + 'Power Capacity IAM/Dossier s0 ' + self.model_s0 + '_' + self.scenario

        # Import the matrix of correspondence between Techno and Techno Init
        self.TechnoMatrix = pd.read_excel(self.folder_path + 'Techno Matrix.xlsx', index_col=0)

        # Import the matrix of correspondence between TechnoEV agg and Techno EV disag
        self.TechnoMatrixEV = pd.read_excel(self.folder_path + 'Techno Matrix.xlsx', sheet_name='Techno matrix EV', index_col=0)

        self.MS_Vehicle_Type = pd.read_excel(self.folder_path + 'Market Share in time EV.xlsx', sheet_name='MS_Vehicle',
                                        index_col=0)
        # Import market share of EV battery type, in percentage, by battery, by year
        self.MS_Battery = pd.read_excel(self.folder_path + 'Market Share in time EV.xlsx', sheet_name='MS_Battery', index_col=0)
        # Import market share of EV motor type, in percentage, fixed in time
        self.MS_Motor = pd.read_excel(self.folder_path + 'Market Share in time EV.xlsx', sheet_name='MS_Motor', index_col=0)

        # Import metal intensity in g of metal per vehicle type
        self.MI_VehicleType = pd.read_excel(self.folder_path + 'Metal Intensity EV.xlsx', sheet_name='MI Vehicle', index_col=0)
        # Import metal intensity in g of metal per kWh of battery aggregated
        self.MI_BatteryAgg = pd.read_excel(self.folder_path + 'Metal Intensity EV.xlsx', sheet_name='MI Battery pack ag',
                                      index_col=0)
        # Import metal intensity in g of metal per kWh of battery disaggregated
        self.MI_BatteryDisag = pd.read_excel(self.folder_path + 'Metal Intensity EV.xlsx', sheet_name='MI Battery pack disag',
                                         index_col=0)
        # Import metal intensity in g of metal per kW of motor
        self.MI_Motor = pd.read_excel(self.folder_path + 'Metal Intensity EV.xlsx', sheet_name='MI Motor', index_col=0)

        # Import statistics on vehicles : kWh of battery by vehicle type, kW of motor by vehicle type
        self.Vehicle_stat = pd.read_excel(self.folder_path + 'Metal Intensity EV.xlsx', sheet_name='Vehicle caracteristics',
                                     index_col=0)

        # Read initial recycling rates for EV metals
        self.Initial_Recycling_Rate = pd.read_excel(self.folder_path + 'Recycling_EV.xlsx', index_col=0)

        # Read scenario of evolution of recycling rates in time for EV metals
        self.Recycling_Evolution = pd.read_excel(self.folder_path + 'Recycling_EV.xlsx', sheet_name='Recycling_Evolution',
                                            index_col=0)

    # Calculate Metal Intensity (MI) evolution in time
    def Power_Metal_Intensity(self):
        '''
        Calculation of the evolution of the metal intensity of energy sources technology, according
        to the scenario SSP chosen.
        :return: Metal intensity in tons of metal per GW of energy sources, with evolution in time [t/GW]
        '''

        # Import table of metal intensities for specific sub-technologies [t/GW] according to the MI scenario chosen
        MetalIntensity = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx', sheet_name = 'MI',index_col=0)
        # Import the table of metal intensities for specific sub-technologies
        self.MetalIntensity_Agg = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx', sheet_name='Aggregated MI 2010',
                                           index_col=0)

        # Creation of a dictionary for the Metal Intensity of technologies over time
        MetalIntensity_doc = {}
        # Take the initial data for the first decade metal intensity
        MetalIntensity_doc[self.listDecades[0]] = MetalIntensity

        # Tab of different scenario of metal intensity reduction
        ReductionScenario = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx', sheet_name='MI reduction scenario',
                                          index_col=0)
        # Load categorization of metals : bulk vs techno-specific
        MetalCategorisation = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx',
                                            sheet_name='Metal categorisation', index_col=0)
        # Create sub-list of metals for each category according to the table categorization
        listBulkMetals = [m for m in self.listMetals if MetalCategorisation['Bulk metals'].loc[m] == 'x']
        listTechnoSpecificMetals = [m for m in self.listMetals if MetalCategorisation['Technology-specific'].loc[m] == 'x']

        # Create a variable to stock the scenario of metal intensity reduction, according to the scenario ssp
        MI_reduc_scenario = str()
        if self.scenario.startswith('SSP1'):
            MI_reduc_scenario = 'Optimistic'
        elif self.scenario.startswith('SSP2') or self.scenario.startswith('SSP5'):
            MI_reduc_scenario = 'Neutral'
        elif self.scenario.startswith('SSP3') or self.scenario.startswith('SSP4'):
            MI_reduc_scenario = 'Conservative'


        # Reductions for each category of metals, according to the scenario ssp
        bulk_reduction = (100 - ReductionScenario[MI_reduc_scenario].loc['Bulk metals']) / 100
        techno_specific_reduction = (100 - ReductionScenario[MI_reduc_scenario].loc['Technology-specific']) / 100
        # Iterate over the years starting from the first consecutive year
        for y in range(1, len(self.listYears)):
            current_year = self.listYears[y]
            previous_year = self.listYears[y - 1]

            # Create a DataFrame for the current year
            MetalIntensity_doc[current_year] = pd.DataFrame(index=self.listMetals, columns=self.listTechno)

            for t in self.listTechno:
                # Apply the MI reduction in time for bulk metals
                MetalIntensity_doc[current_year].loc[listBulkMetals, t] = (
                        MetalIntensity_doc[previous_year].loc[listBulkMetals, t] * bulk_reduction
                )
                # Apply the MI reduction in time for technology-specific metals
                MetalIntensity_doc[current_year].loc[listTechnoSpecificMetals, t] = (
                        MetalIntensity_doc[previous_year].loc[listTechnoSpecificMetals, t] * techno_specific_reduction
                )
        # Final need in metal m for the techno t at the decade d
        # in [t/GW] with reduction in time according to ssp
        self.MetalIntensity_doc = MetalIntensity_doc

    def Power_Charging_Factor(self):
        '''
        Calculate Charging Factor (CF) by sub-technology, based on initial IAM CF
        '''

        # Define the folder with CF data from IAM
        CF_folder = self.folder_path + 'Capacity Factor IAM/FC ' + self.model_s0 + '_' + self.scenario

        # Reads the folder with Excel files of CF from IAM dataset
        excelCF = [f for f in os.listdir(CF_folder) if f.endswith('xlsx')]

        # Creation of a dictionary of dataFrame to stock final CF
        CF = {}
        # Copies initially calculated CF by IAM in the dictionary
        for region_name, file in zip(self.listRegions, excelCF):
            fileCFbyRegion = os.path.join(CF_folder, file)  # copies datas from files in dataFrames
            CF[region_name] = pd.read_excel(fileCFbyRegion,index_col=0)  # copies the dataFrames in dic
            CF[region_name] = CF[region_name].drop('Total', axis=0) # Drop the total column

        # List of techno included in the initial CF
        list_CF_techno = list(CF[self.listRegions[0]].index)
        self.list_CF_techno =list_CF_techno

        # Correct the CF for missing technology from IAM data
        for region in self.listRegions:
            # If Hydro not in IAM data :
            if "Hydro" not in list_CF_techno:
                # Add CF from the literature
                CF[region].loc["Hydro"] = self.CF_from_litt["CF from literature"].loc["Hydro"] * 100  # Add a line for hydro when not included, datas from the literature
            # If Oil not in IAM data :
            if "Oil" not in list_CF_techno:
                # Add CF from the literature
                CF[region].loc["Oil"] = self.CF_from_litt["CF from literature"].loc["Oil"] * 100  # Add a line for oil when not included, datas from the literaturz
            # Correction of unlikely Biomass CF from IAM
            for decade in self.listDecades:
                # If CF for Biomass is more than 85%
                if CF[region][decade].loc["Biomass"] >= 0.85:
                    # Replace CF from IAM by CF from the literature
                    CF[region][decade].loc["Biomass"] = self.CF_from_litt["CF from literature"].loc["Biomass"] * 100
            # If Geothermal not in IAM data :
            if "Geothermal" not in list_CF_techno:
               # Hypothesis of Geothermal having the same CF from IAM data than Biomass
               CF[region].loc["Geothermal"] = CF[region].loc["Biomass"]

            if "Wind Offshore" not in list_CF_techno:
                CF[region].loc["Wind Offshore"] = self.CF_from_litt["CF from literature"].loc[
                                                      "Wind Offshore"] * 100  # Add a line for wind offshore from the literature

            CF[region] = CF[region].sort_index(axis=0) # Sort CF techno

        # Add CF from literature for missing data
        for region in self.listRegions:
            for techno in list_CF_techno : # Loop through years included in the initial CF
                for year in list(CF[self.listRegions[0]].columns): # Loop through techno included in the initial CF
                    if np.isnan(CF[region][year].loc[techno]):  # Replace if the value is missing
                        CF[region][year].loc[techno] = self.CF_from_litt["CF from literature"].loc[
                                                           techno] * 100  # Take values from the litt (*100 bc percentage)

            # Rename the line "Wind" in "Wind Onshore"
            CF[region] = CF[region].rename(index={
                "Wind": "Wind Onshore"})  # Wind Onshore is 90% of Wind. We suppose IAM data are applied to Wind Onshore

        # Initialise dictionary to store disaggregated CF datas
        CF_disag = {}
        for region in self.listRegions:
            # Transpose the matrix with the corresponding CF techno / sub-techno matrix
            CF_disag[region] = pd.DataFrame(self.Matrix_CF_disag.transpose() @ CF[region])

        # Change the CF of CSP for those of the literature, for the decades studied
        for region in self.listRegions:
            for years in self.listDecades:
                for techno in self.listTechnoAgg:
                    # if the techno is a sub-techno of solar CSP
                    if "CSP" in techno:
                        # replace the CF by the one from the literature
                        CF_disag[region][years][techno] = self.CF_from_litt["CF from literature"].loc["CSP"] * 100

        # Change the CF of solar to be constant in time, for the decades studied
        for region in self.listRegions:
            for years in self.listDecades[1:]: # exclude the first decade
                for techno in self.listTechnoAgg:
                    # if the techno is a sub-techno of solar PV
                    if "Sol" in techno:
                        # replace the CF by the one from IAM at the first decade
                        CF_disag[region][years][techno] = CF_disag[region][self.listDecades[0]][techno]
        # Final regionalized CF for the disaggregated sub-technologies
        self.CF_disag = CF_disag

    def Power_Capacities(self):
        '''
        Collect and disaggregate initial IAM power capacity
        '''

        # Loop through all files with IAM power capacity data
        for file in os.listdir(self.folderIAM):
            if file.endswith('.xlsx'):  # Check if the file is an Excel file
                fileIAMbyRegion = os.path.join(self.folderIAM, file)  # Complete file path
                # Read the Excel file into a DataFrame with the first column as index
                df = pd.read_excel(fileIAMbyRegion)
        # Creation of the initial IAM dictionary
        excelRegionsIAM = [f for f in os.listdir(self.folderIAM) if
                           f.endswith('xlsx')]  # Reads the folder with excel files from IAM dataset
        # Creation of a dictionary of dataFrame to stock datas from Excel files
        IAM_initial_data = {}
        for region_name, file in zip(self.listRegions, excelRegionsIAM):
            fileIAMbyRegion = os.path.join(self.folderIAM, file)  # copies datas from files in dataFrames
            IAM_initial_data[region_name] = pd.read_excel(fileIAMbyRegion,
                                                          index_col=0)  # copies the dataFrames in dic, with first col as index

        # Delete lines with power capacity repetition (ex. Capacity total, Capacity wind and onshore/offshore)
        # Identify the lines to delete from the IAM Power Capacity, according to the IAM chosen
        if (self.model_s0 == 'AIM..CGE'):
            excelLinesToDelete = [2, 10, 12]
        if (self.model_s0 == 'GCAM4'):
            excelLinesToDelete = [2, 9]
        if (self.model_s0 == 'IMAGE'):
            excelLinesToDelete = [2, 9, 10, 13]
        if (self.model_s0 == 'MESSAGE-GLOBIOM'):
            excelLinesToDelete = [2, 10, 13]
        if (self.model_s0 == 'REMIND-MAGPIE'):
            excelLinesToDelete = [2, 9]
        if (self.model_s0 == 'WITCH-GLOBIOM'):
            excelLinesToDelete = [2, 9, 10, 13]

        # List of lines to delete from the excel files (total of capacity, "other", repetitions...)
        LinesToDelete = [line - 2 for line in excelLinesToDelete]
        for region in self.listRegions:
            # Delete excessive lines in the dictionary
            IAM_initial_data[region] = IAM_initial_data[region].drop(IAM_initial_data[region].index[LinesToDelete],
                                                                     axis=0)

        # Stock Market Share (MS) matrix of the IAM techno and the corresponding disaggregated sub-techno
        MarketShareInit_df = {}
        # Import MS data for each decade
        for d in self.listDecadeswPast:
            name_sheet = d
            MarketShareInit_df[d] = pd.read_excel(f'{self.folder_path}Market Share in time - {self.ModelisationType}.xlsx', sheet_name=name_sheet,
                                                      index_col=0)

        self.MarketShareInit_df = MarketShareInit_df

        # Align IAM initial techno and matrix of MS, according to the techno represented by the model_s0 studied
        MarketShare_df = {}
        for y in self.listDecadeswPast: # Loop through decades studied + first previous decade
            # Filter the rows of MarketShareInit_df based on the list of indexes from the initial IAM datas
            MarketShare_df[y] = MarketShareInit_df[y][MarketShareInit_df[y].index.isin(IAM_initial_data[self.listRegions[0]].index)]

        # Initialise dictionary to store disaggregated IAM power capacity, for each region, techno and decade [GW]
        self.s0 = {}
        for r in self.listRegions:
            # Initialize a dictionary for the current region
            self.s0[r] = pd.DataFrame()
            for y in self.listDecadeswPast:
                # Transpose initial IAM power cap and MS matrix of sub-techno correspondence
                self.s0[r][y] = pd.DataFrame(MarketShare_df[y].transpose() @ IAM_initial_data[r][y])

        if self.ModelisationType =='Opti':
            self.listTechnoIAM = self.listTechnoAgg
        if self.ModelisationType =='Init':
            self.listTechnoIAM = self.listTechno

    def Demand_Network_Storage(self):
        '''
        Calculate the demand in metal for the rest of the transition : electric vehicle, grid, storage
        Based on different International Energy Agency (IEA) scenarios
        '''

        # 0. Match IEA scenario with SSP-RCP ones

        # Addition of the data from the IEA in this DataFrame. The IEA scenario considered depends on the SSP
        # Net Zero
        if self.scenario == 'SSP1-19' or self.scenario == 'SSP4-19':
            scenarioIEA = 'NZE'
            self.scenarioIEA = scenarioIEA
        # ABS
        elif self.scenario == 'SSP1-26' or self.scenario == 'SSP4-26':
            scenarioIEA = 'APS'
            self.scenarioIEA = scenarioIEA
        # Stated policy
        else:
            scenarioIEA = 'SPS'
            self.scenarioIEA = scenarioIEA

        # 1. Demand for Storage of energy, from IEA data

        # Select only the data of the scenario studied, thanks to the column scenario
        self.Storage_Demand_init = self.Storage_Demand_init[self.Storage_Demand_init["Scenario"] == scenarioIEA]

        # Add missing years to the initial dataFrame
        all_years = np.arange(2020, 2051)
        self.Storage_Demand_init = self.Storage_Demand_init.reindex(columns=all_years)

        # Linear Interpolation to fill intermediate years
        for m in self.Storage_Demand_init.index:
            self.Storage_Demand_init.loc[m] = self.Storage_Demand_init.loc[m].interpolate(method='linear')
        # Linear extrapolation for 2020 and 2021
        self.Storage_Demand_init[2020] = 2 * self.Storage_Demand_init[2022] - self.Storage_Demand_init[2025]
        self.Storage_Demand_init[2021] = (self.Storage_Demand_init[2022] + self.Storage_Demand_init[2020]) / 2  # Mean value
        # Change columns as strings
        self.Storage_Demand_init.columns = self.Storage_Demand_init.columns.astype(str)

        # Create the final Storage_Demand dataFrame
        self.Storage_Demand = pd.DataFrame(0.0, index=self.listMetals, columns=self.listYears)
        for y in self.listYears:
            for m in self.Storage_Demand_init.index:
                if self.Storage_Demand_init[y].loc[m] > 0.0:
                    self.Storage_Demand[y].loc[m] = self.Storage_Demand_init[y].loc[m]

        # 2. Demand for Grid network, from IEA data

        # Select only the data of the scenario studied, thanks to the column scenario
        self.Network_Demand_init = self.Network_Demand_init[self.Network_Demand_init["Scenario"] == scenarioIEA]

        # Linearization of metal demand for missing years
        self.Network_Demand_init = self.Network_Demand_init.reindex(columns=all_years)

        # Linear Interpolation to fill intermediate years
        for m in self.Network_Demand_init.index:
            self.Network_Demand_init.loc[m] = self.Network_Demand_init.loc[m].interpolate(method='linear')

        # Linear extrapolation for 2020 and 2021
        self.Network_Demand_init[2020] = 2 * self.Network_Demand_init[2022] - self.Network_Demand_init[2025]
        self.Network_Demand_init[2021] = (self.Network_Demand_init[2022] + self.Network_Demand_init[2020]) / 2  # Mean value

        # Change columns as strings
        self.Network_Demand_init.columns = self.Network_Demand_init.columns.astype(str)

        self.Network_Demand = pd.DataFrame(0.0, index=self.listMetals, columns=self.listYears)

        for y in self.listYears:
            for m in self.Network_Demand_init.index:
                if self.Network_Demand_init[y].loc[m] > 0.0:
                    self.Network_Demand[y].loc[m] = self.Network_Demand_init[y].loc[m]


    def Future_Production(self):

        '''
        Estimation of future metal production in tons by year, until 2050
        '''

        # Initialisation of a dataframe with metals in lines and years in columns
        Prod = pd.DataFrame(index=self.listMetals, columns=self.listYears)

        # Initialise the datas from 2020 with actual known 2020 metal production
        Prod["2020"] = self.MetalData['Production 2020 [t/yr]']
        # Historical annual growth production rate by metal
        Beta = pd.read_excel(self.folder_path + 'Future Metal Production.xlsx', sheet_name = 'Beta historic', index_col=0)
        # Initialise a final demand dF to stock future demand in metals, with GDP increase
        ProdFinalbyYear = {year: Prod[year].copy() for year in self.listYears}

        for m in self.listMetals:
            for y in range(1, len(self.listYears)):
                # Loop through list of years, excluding the first year 2020
                previous_prod = ProdFinalbyYear[self.listYears[y - 1]].loc[
                    m]  # Create a variable for the production at the previous year
                ProdFinalbyYear[self.listYears[y]].loc[m] = previous_prod + Beta['Mean annual growth'].loc[
                    m] * previous_prod  # Calculate production with beta increase by year

        # Conversion of the dic in dF
        ProdFinalbyYear = pd.DataFrame.from_dict(ProdFinalbyYear, orient='index')

        # Transpose the DataFrame to have metals as columns and years as index
        ProdFinalbyYear = ProdFinalbyYear.transpose()

        # Create the dataFrame with metal production other each decade
        Prod = pd.DataFrame({'2020': self.MetalData['Production 2020 [t/yr]'], '2030': ProdFinalbyYear['2030'],
                             '2040': ProdFinalbyYear['2040'], '2050': ProdFinalbyYear['2050']})

        # Correction of future prod with IEA High production case projections for Cobalt, Copper, Nickel
        Prod_litt = pd.read_excel(self.folder_path + 'Future metal production.xlsx', sheet_name='Other prod estimates', index_col=0)

        for m in Prod_litt.index:
            for d in self.listDecades:
                Prod[d].loc[m] = Prod_litt[d].loc[m]

        self.Prod = Prod

    def EV_Market_Share(self):

        # list Years from 2005 to 2050 with a step of 5
        self.listYearsTot = [str(year) for year in range(2005, 2051, 1)]

        self.Vehicle_Stock = pd.Series()

        if self.model_s0 in ['IMAGE', 'GCAM4']:
            # LDV Ratio
            Ratio_LDV = pd.read_excel(self.folder_path + 'Market Share in time EV.xlsx', sheet_name='Ratio_LDV', index_col=0)

            # initialize Stock_IAM
            Stock_IAM = pd.Series(dtype=float)

            # Construction du chemin vers le fichier Transport
            path = f"{self.folder_path}Transport IAM/Transport_{self.model_s0}_{self.scenario}.xlsx"

            # Vérification de l'existence du fichier et traitement
            if os.path.exists(path):
                table = pd.read_excel(path)
                table.drop(['VARIABLE', 'REGION'], axis=1, inplace=True)
                Stock_IAM = table.sum()

            # Linear interpolation for missing years between decades
            Stock_IAM = Stock_IAM.reindex(self.listYearsTot)
            Stock_IAM = Stock_IAM.interpolate(method='linear')
            # In million of vehicles, to match MI in g/vehicles
            Stock_IAM = Stock_IAM * 10 ** 3 / 19300

            for y in self.listYearsTot:
                self.Vehicle_Stock[y] = Stock_IAM[y] * Ratio_LDV.loc[self.model_s0, 'LDV_Share']

        else:

            # Import data of stock of passenger vehicles from MATILDA
            Stock_MATILDA = pd.read_excel(self.folder_path + 'Market Share in time EV.xlsx', sheet_name='Vehicle_Stock',
                                          index_col=0)

            # Take only Vehicle Demand for the global region
            Stock_MATILDA = Stock_MATILDA[Stock_MATILDA["Region"] == 'GLOBAL']
            # The low scenario is chosen, because it is the closest to IAM estimations
            Stock_MATILDA = Stock_MATILDA[Stock_MATILDA["Stock_scenario"] == 'Low']

            # Change time values to string
            Stock_MATILDA['Time'] = Stock_MATILDA['Time'].astype(str)

            for y in self.listYearsTot:
                # In million of vehicles, to match MI in g/vehicles
                self.Vehicle_Stock[y] = Stock_MATILDA.loc[(Stock_MATILDA["Time"] == y), "value"].values[0] / 10 ** 6

        scenarioIEA_VE = str()

        # Scenario for MS of vehicles from IEA
        if self.scenarioIEA == 'NZE':
            scenarioIEA_VE = 'Net Zero'
        elif self.scenarioIEA == 'APS':
            scenarioIEA_VE = 'SD'
        else:
            scenarioIEA_VE = 'STEP'

        # Select only the data of the scenario studied, thanks to the column scenario
        self.MS_Vehicle_Type = self.MS_Vehicle_Type[self.MS_Vehicle_Type["scenario"] == scenarioIEA_VE]
        # Change time values to string
        self.MS_Vehicle_Type['year'] = self.MS_Vehicle_Type['year'].astype(str)

        # New installation of vehicles vT at the year y because of growth
        VehicleType_Growth = pd.DataFrame(index=self.listVehicleType, columns=self.listYearsTot[1:])
        for y in range(1, len(self.listYearsTot)):
            for vT in self.listVehicleType:
                # Take the maximum value between the additional stock and 0, if there is a decrease in stock
                VehicleType_Growth.loc[vT, self.listYearsTot[y]] = max(self.Vehicle_Stock[self.listYearsTot[y]] * self.MS_Vehicle_Type.loc[
                    (self.MS_Vehicle_Type["drive_train"] == vT) & (
                                self.MS_Vehicle_Type["year"] == self.listYearsTot[y]), "value"].values[0]
                                                                  - (self.Vehicle_Stock[self.listYearsTot[y - 1]] *
                                                                     self.MS_Vehicle_Type.loc[
                                                                         (self.MS_Vehicle_Type["drive_train"] == vT) & (
                                                                                     self.MS_Vehicle_Type["year"] ==
                                                                                     self.listYearsTot[
                                                                                         y - 1]), "value"].values[0])
                                                                  , 0)
        self.VehicleType_Growth = VehicleType_Growth

        # New installation of vehicles to maintain the stock
        VehicleType_Maintain = pd.DataFrame(index=self.listVehicleType, columns=self.listYearsTot[13:])
        for y in range(13, len(self.listYearsTot)):
            for vT in self.listVehicleType:
                VehicleType_Maintain.loc[vT, self.listYearsTot[y]] = self.VehicleType_Growth.loc[vT, self.listYearsTot[y - self.Lifetime_V]]

        # Create the dataFrame for sales of vehicles by vehicle type by year, from 2018 to 2050
        VehicleType_Sales = pd.DataFrame(index=self.listVehicleType, columns=self.listYearsTot[13:])

        for vT in self.listVehicleType:
            for y in self.listYearsTot[13:]:
                # Stock of total vehicles at the year y * Market share by vehicle type v at the year y
                VehicleType_Sales.loc[vT, y] = self.VehicleType_Growth.loc[vT, y] + VehicleType_Maintain.loc[vT, y]


        # Years from 2018 to 2050
        self.listYearsVehicle = self.listYearsTot[13:]


        # Scenario for MS of batteries from IRENA
        scenarioIRENA_battery = str()
        if self.scenario.startswith('SSP1'):
            scenarioIRENA_battery = 'Increased Innovation Scenario'
        elif self.scenario.startswith('SSP2') or self.scenario.startswith('SSP5'):
            scenarioIRENA_battery = 'Current Trend Scenario'
        else:
            scenarioIRENA_battery = 'Technology Stagnation Scenario'

        # Select only the data of the scenario studied, thanks to the column scenario
        self.MS_Battery = self.MS_Battery[self.MS_Battery["Scenario"] == scenarioIRENA_battery]
        '''
        # Create an empty DataFrame with defined columns
        Vehicle_Sales = pd.DataFrame(0, index=[], columns=self.listYearsVehicle)

        for v in self.listVehicleType:
            if v == 'ICEV':
                # Directly copy data for ICEV (Internal Combustion Engine Vehicle)
                Vehicle_Sales.loc[v] = VehicleType_Sales.loc[v]
            else:
                for b in self.listBatteryAgg:
                    # Create names for the PM (Permanent Magnet) and Induction battery configurations
                    v_PM_Batt = f"{v}_PM_{b}"
                    v_ind_Batt = f"{v}_Ind_{b}"

                    # Add new rows with initial values (set to 0)
                    Vehicle_Sales.loc[v_PM_Batt] = 0
                    Vehicle_Sales.loc[v_ind_Batt] = 0

                    for y in self.listYearsVehicle:
                        # Calculate and assign values based on market share and battery data
                        Vehicle_Sales.loc[v_PM_Batt, y] = (
                                VehicleType_Sales.loc[v, y]
                                * self.MS_Motor.loc['MS', 'Motor_PM']
                                * self.MS_Battery.loc[b, y]
                        )

                        Vehicle_Sales.loc[v_ind_Batt, y] = (
                                VehicleType_Sales.loc[v, y]
                                * self.MS_Motor.loc['MS', 'Motor_induction']
                                * self.MS_Battery.loc[b, y]
                        )
        '''

        # Create a dF with the demand in stock of vehicle, for each years
        Vehicle_Growth = pd.DataFrame(0, index=[], columns=self.listYearsTot[1:])

        for v in self.listVehicleType:
            if v == 'ICEV':
                for y in self.listYearsTot[1:]:
                    # Copier directement les données pour ICEV
                    Vehicle_Growth.loc[v, y] = VehicleType_Growth.loc[v, y]
            else:
                for b in self.listBatteryAgg:
                    # Créer les noms pour les véhicules PM et Ind
                    v_PM_Batt = f"{v}_PM_{b}"
                    v_ind_Batt = f"{v}_Ind_{b}"

                    # Ajouter les nouvelles lignes avec des valeurs initiales
                    Vehicle_Growth.loc[v_PM_Batt] = 0
                    Vehicle_Growth.loc[v_ind_Batt] = 0

                    for y in self.listYearsVehicle:
                        # Calculer et assigner les valeurs directement
                        Vehicle_Growth.loc[v_PM_Batt, y] = (
                                VehicleType_Growth.loc[v, y]
                                * self.MS_Motor.loc['MS', 'Motor_PM']
                                * self.MS_Battery.loc[b, y]
                        )

                        Vehicle_Growth.loc[v_ind_Batt, y] = (
                                VehicleType_Growth.loc[v, y]
                                * self.MS_Motor.loc['MS', 'Motor_induction']
                                * self.MS_Battery.loc[b, y]
                        )
        self.Vehicle_Growth = Vehicle_Growth

        # list of vehicles, with vehicle type, motor type, and aggregated battery types
        self.listVehicleAgg = self.Vehicle_Growth.index.tolist()

        # New installation of vehicles to maintain the stock
        self.Vehicle_Maintain = pd.DataFrame(index=self.listVehicleAgg, columns=self.listYearsTot[13:])
        for y in range(13, len(self.listYearsTot)):
            for v in self.listVehicleAgg:
                self.Vehicle_Maintain.loc[v, self.listYearsTot[y]] = Vehicle_Growth.loc[
                    v, self.listYearsTot[y - self.Lifetime_V]]

        # Create the dataFrame for sales of vehicles type by year, from 2018 to 2050
        self.Vehicle_Sales = pd.DataFrame(index=self.listVehicleAgg, columns=self.listYearsVehicle)

        for v in self.listVehicleAgg:
            for y in self.listYearsVehicle:
                # Stock of total vehicles at the year y * Market share by vehicle type v at the year y
                self.Vehicle_Sales.loc[v, y] = self.Vehicle_Growth.loc[v, y] + self.Vehicle_Maintain.loc[v, y]

        # Create x0, the initial scenario for EV - motor - battery demand
        self.x0 = self.Vehicle_Sales

        # list of vehicles, with vehicle type, motor type, and precise battery types
        listVehicleDisag = []
        for v in self.listVehicleType:
            if v == 'ICEV':
                # Copier directement les données pour ICEV
                listVehicleDisag.append('ICEV')
            else:
                for b in self.listBatteryDisag:
                    listVehicleDisag.append(f"{v}_PM_{b}")
                    listVehicleDisag.append(f"{v}_Ind_{b}")
        self.listVehicleDisag = listVehicleDisag

    def EV_Metal_Intensity(self):

        # Initial scenario uses aggregated battery data
        if self.ModelisationType == 'Init':
            self.MI_Battery = self.MI_BatteryAgg
            self.listBattery = self.listBatteryAgg
            self.listVehicle = self.listVehicleAgg
        # Optimised scenario uses disaggregated battery data
        elif self.ModelisationType == 'Opti':
            self.MI_Battery = self.MI_BatteryDisag
            self.listBattery = self.listBatteryDisag
            self.listVehicle = self.listVehicleDisag

        # Metal intensity in g of metal per vehicle type with agg or disag battery type
        # According to the Modelisation Type chosen
        MI_Vehicle = pd.DataFrame(0.0, columns=self.listVehicle, index=self.listMetals)

        for v in self.listVehicle:
            # Add metal intensities of vehicles by vehicle type
            for v_agg in self.listVehicleType:
                if v_agg in v:
                    for m in self.MI_VehicleType.index:
                        MI_Vehicle.loc[m][v] += self.MI_VehicleType.loc[m][v_agg]
            # Add metal intensities of vehicles by battery type
            for b in self.listBattery:
                if b in v:
                    for v_agg in self.listVehicleType:
                        if v_agg in v:
                            for m in self.MI_Battery.index:
                                MI_Vehicle.loc[m][v] += self.MI_Battery.loc[m][b] * self.Vehicle_stat.loc['Battery'][v_agg]
            # Add metal intensities of vehicles by motor type
            for mo in self.listMotor:
                if mo in v:
                    for v_agg in self.listVehicleType:
                        if v_agg in v:
                            for m in self.MI_Motor.index:
                                MI_Vehicle.loc[m][v] += self.MI_Motor.loc[m][mo]
        self.MI_Vehicle = MI_Vehicle

        # Metal intensity in g of metal per vehicle type, with aggregated battery type
        MI_VehicleAgg = pd.DataFrame(0.0, columns=self.listVehicleAgg, index=self.listMetals)

        for v in self.listVehicleAgg:
            # Add metal intensities of vehicles by vehicle type
            for v_agg in self.listVehicleType:
                if v_agg in v:
                    for m in self.MI_VehicleType.index:
                        MI_VehicleAgg.loc[m][v] += self.MI_VehicleType.loc[m][v_agg]
            # Add metal intensities of vehicles by battery type
            for b in self.listBatteryAgg:
                if b in v:
                    for v_agg in self.listVehicleType:
                        if v_agg in v:
                            for m in self.MI_BatteryAgg.index:
                                MI_VehicleAgg.loc[m][v] += self.MI_BatteryAgg.loc[m][b] * self.Vehicle_stat.loc['Battery'][
                                    v_agg]
            # Add metal intensities of vehicles by motor type
            for mo in self.listMotor:
                if mo in v:
                    for v_agg in self.listVehicleType:
                        if v_agg in v:
                            for m in self.MI_Motor.index:
                                MI_VehicleAgg.loc[m][v] += self.MI_Motor.loc[m][mo]
        self.MI_VehicleAgg = MI_VehicleAgg

    def EV_Recycling(self):
        '''
        Estimation of the recycling rates from 2020 to 2050 of metals used in electric vehicles,
        according to the scenario SSP considered
        '''

        # DataFrame for recycling rates of metals used in EV
        self.Recycling_Rate = pd.DataFrame(columns=self.listYears, index=self.listMetals)

        # Change recycling_Evolution df columns in strings
        self.Recycling_Evolution.columns = self.Recycling_Evolution.columns.astype(str)

        # Add initial recycling rates for the 2020 year
        for m in self.listMetals:
            if m in self.Initial_Recycling_Rate.index:
                self.Recycling_Rate['2020'].loc[m] = self.Initial_Recycling_Rate['Initial_Recycling_Rate'].loc[m]
            else:
                self.Recycling_Rate['2020'].loc[m] = 0.0

        # Add the assumed recycling rate in 2050
        if 'SSP4' in self.scenario or 'SSP5' in self.scenario:
            for m in self.listMetals:
                if m in self.Initial_Recycling_Rate.index:
                    # Constant recycling rate
                    self.Recycling_Rate['2050'].loc[m] = self.Initial_Recycling_Rate['Initial_Recycling_Rate'].loc[m]
                else:
                    self.Recycling_Rate['2050'].loc[m] = 0.0
        else:
            for SSPs in ['SSP1', 'SSP2', 'SSP3']:
                if SSPs in self.scenario:
                    # Select data for the concerned scenario
                    Recycling_EvolutionSSP = self.Recycling_Evolution.loc[self.Recycling_Evolution['SSPs'] == SSPs]
                    for m in self.listMetals:
                        # Add recycling rates for concerned metals only
                        if m in self.Initial_Recycling_Rate.index:
                            # Improved recycling rate according to the scenario
                            self.Recycling_Rate['2028'].loc[m] = Recycling_EvolutionSSP['2028'].loc[m]
                            self.Recycling_Rate['2032'].loc[m] = Recycling_EvolutionSSP['2032'].loc[m]
                            self.Recycling_Rate['2050'].loc[m] = Recycling_EvolutionSSP['2050'].loc[m]
                        else:
                            self.Recycling_Rate['2050'].loc[m] = 0.0

        # Conversion of years in int for the interpolation
        self.Recycling_Rate.columns = self.Recycling_Rate.columns.astype(int)
        # Convert all columns to numeric (even if they are of object type)
        self.Recycling_Rate = self.Recycling_Rate.apply(pd.to_numeric, errors='coerce')
        # Linear interpolation between the initial 2020 value, and the 2050 assumed one
        self.Recycling_Rate = self.Recycling_Rate.interpolate(method='linear', axis=1)
        # Conversion of years in string
        self.Recycling_Rate.columns = self.Recycling_Rate.columns.astype(str)

    def Demand_Other_Sectors(self, growth):
        '''
        Estimation of the metal demand of the rest of the economy
        Based on literature and Gross Domestic Product (GDP) growth from IAM and SSP-RCP
        '''

        # Initialisation of a dF for other sector demand in metals by year, with metals in lines and years in columns
        OSD = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)

        # 1. Estimate OSD by year in 2020

        # New capacity installed between 2010 and 2020, by techno
        new2020Capacity_d = pd.DataFrame(0.0, index=self.listTechnoAgg, columns=['2020'])
        for t in self.listTechnoIAM:
            for r in self.listRegions:
                # Only account if there is an increase in capacity
                if self.s0[r].loc[t, '2020'] - self.s0[r].loc[t, '2010'] > 0:
                    # The new installed capacity in the decade between 2020 and 2010
                    new2020Capacity_d['2020'].loc[t] = sum(
                        self.s0[r]['2020'].loc[t] - self.s0[r]['2010'].loc[t] for r in self.listRegions)

        # Metal demand for the new installed capacity at the year 2020
        PowerMetalDemand_y = pd.DataFrame(0.0, index=self.listMetals, columns=['2020'])
        for m in self.listMetals:
            for t in self.listTechnoAgg:
                if t == 'Sol_Thin_Film' or t == 'Wind_Onshore':
                    # Aggregated MI data, with 2010 market share
                    PowerMetalDemand_y['2020'].loc[m] = PowerMetalDemand_y['2020'].loc[m] + \
                                                        new2020Capacity_d['2020'].loc[t] / 10 * \
                                                        self.MetalIntensity_Agg[t].loc[m]
                elif t == 'Sol_C-si':  # All Sol_C-Si techno is based on silver before 2020
                    PowerMetalDemand_y['2020'].loc[m] = PowerMetalDemand_y['2020'].loc[m] + \
                                                        new2020Capacity_d['2020'].loc[t] / 10 * \
                                                        self.MetalIntensity_doc['2020']['Sol_C-si_Silver'].loc[m]
                else:
                    PowerMetalDemand_y['2020'].loc[m] = PowerMetalDemand_y['2020'].loc[m] + \
                                                        new2020Capacity_d['2020'].loc[t] / 10 * \
                                                        self.MetalIntensity_doc['2020'][t].loc[m]
        self.PowerMetalDemand_y = PowerMetalDemand_y

        # Initialise the datas from 2020 with actual known metal prod in 2020 - metal demand for the transition [t/y]
        for m in self.listMetals:
            # Total prod in 2020 - Metal for grid and storage - Metal for powCap installation in 2020 - Metal for EV
            InitialOSD = ((self.Prod['2020'].loc[m] - (self.Storage_Demand['2020'].loc[m]+self.Network_Demand['2020'].loc[m]) / self.MetalData['RR (%) Prod'].loc[m]
                          -self.PowerMetalDemand_y['2020'].loc[m] / self.MetalData['RR (%) Prod'].loc[m])
                          - sum(self.x0.loc[v, '2020'] * self.MI_VehicleAgg.loc[m][v] / self.MetalData['RR (%) Prod'].loc[m] for v in self.listVehicleAgg))
            # Change OSD value only if it is above 0, if not : = 0.
            if InitialOSD > 0:
                OSD["2020"].loc[m] = InitialOSD

        if growth == 'GDP':
            # 2.1 Estimate OSD by year until 2050, based on GDP growth

            # Creation of the GDP dictionary
            excelGDP = [f for f in os.listdir(self.GDP_folder) if
                        f.endswith('xlsx')]  # Reads the folder with excel files of GDP from IAM dataset
            # Creation of a dictionary of dataFrame to stock GDP datas from IAM Excel files
            self.GDP = {}
            for region_name, file in zip(self.listRegions, excelGDP):
                fileGDPbyRegion = os.path.join(self.GDP_folder, file)  # copies datas from GDP files by region in dF
                self.GDP[region_name] = pd.read_excel(fileGDPbyRegion,
                                                      index_col=0)  # copies the dataFrames in dic, with first col as index

            # Initialise a final demand dF to stock future demand in metals, with GDP increase
            for m in self.listMetals:
                # Loop through list of years, excluding the first year 2020
                for d in range(1, len(self.listDecades)):
                    # Create a variable for the demand at the previous year
                    previous_demand = OSD[self.listDecades[d - 1]].loc[m]
                    # Calculate the World GDP increase with a sum for every region
                    GDP_growth = sum(self.GDP[r][self.listDecades[d]].loc['GDP..PPP'] for r in self.listRegions) / sum(
                        self.GDP[r][self.listDecades[d - 1]].loc['GDP..PPP'] for r in self.listRegions) - 1.3/100*10
                    # Calculate Final demand with the previous one and GDP increase
                    if GDP_growth > 0 :
                        OSD[self.listDecades[d]].loc[m] = previous_demand * GDP_growth
                    else :
                        OSD[self.listDecades[d]].loc[m] = previous_demand

        elif growth == 'Pop':

            # 2.2 Estimate OSD by year until 2050, based on Population growth

            # Creation of the Population dictionary
            excelPop = [f for f in os.listdir(self.Population_folder) if
                        f.endswith('xlsx')]  # Reads the folder with excel files of Pop from IAM dataset

            # Creation of a dictionary of dataFrame to stock datas from Excel files
            Pop = {}

            for region_name, file in zip(self.listRegions, excelPop):
                filePopbyRegion = os.path.join(self.Population_folder, file)  # copies datas from files in dataFrames
                Pop[region_name] = pd.read_excel(filePopbyRegion,
                                                 index_col=0)  # copies the dataFrames in dic, with first col as index

            # Initialise a final demand dF to stock future demand in metals, with Pop increase
            # Demand_byYear = {year: Demand[year].copy() for year in listDecades}

            for m in self.listMetals:
                # Loop through list of years, excluding the first year 2020
                for d in range(1, len(self.listDecades)):
                    # Create a variable for the demand at the previous year
                    previous_demand = OSD[self.listDecades[d - 1]].loc[m]
                    # Calculate the World Pop increase with a sum for every region
                    Pop_growth = sum(Pop[r][self.listDecades[d]].loc['Population'] for r in self.listRegions) / sum(
                        Pop[r][self.listDecades[d - 1]].loc['Population'] for r in self.listRegions)
                    # Calculate Final demand with the previous one and Pop increase
                    OSD[self.listDecades[d]].loc[m] = previous_demand * Pop_growth

        # 3. For available data, replace OSD by data from the literature

        OSD_litt = pd.read_excel(self.folder_path + 'OSD litt.xlsx', index_col=0)

        # Correction with available data for the metal m in the literature of future OSD
        for m in OSD_litt.index:
            for d in self.listDecades[1:]:
                OSD[d].loc[m] = OSD_litt[d].loc[m]

        OSD_litt_byScenario = pd.read_excel(self.folder_path + 'OSD litt.xlsx', sheet_name='OSD by scenario', index_col=0)

        # Correction with projections from IEA, matching scenario
        OSD_ScenarioIEA = OSD_litt_byScenario.loc[OSD_litt_byScenario['Scenario'] == self.scenarioIEA]
        for m in OSD_ScenarioIEA.index:
            for d in self.listDecades[1:]:
                OSD[d].loc[m] = OSD_ScenarioIEA[d].loc[m]

        # Correction with projections for different SSP, matching scenario
        OSD_ScenarioSSP = OSD_litt_byScenario.loc[
            OSD_litt_byScenario['Scenario'] == self.scenario[:4]]  # Match with the corresponding ssp
        for m in OSD_ScenarioSSP.index:
            for d in self.listDecades[1:]:
                OSD[d].loc[m] = OSD_ScenarioSSP[d].loc[m]

        self.OSD = OSD

    def modelDef(self):
        '''
        Definition of the optimisation model and the variables
        '''

        self.model = ConcreteModel()

        # Declaration of the optimization variables

        # Cumulated powder capacities at the region r, of the energy source t, at the decade d
        self.model.s = Var(self.listRegions, self.listTechno, self.listDecades, initialize=0, domain=NonNegativeReals)
        # Millions of specific vehicles v sold at the year y (from 2020 to 2050)
        self.model.x = Var(self.listVehicle, self.listYearsVehicle, initialize=0, domain=NonNegativeReals)
        # Millions of vehicles sold to answer the growth in demand, by vehicle type, at the year y
        self.model.x_growth = Var(self.listVehicleType, self.listYearsTot[1:], initialize=0, domain=NonNegativeReals)
        # Millions of vehicles sold to maintain stock of previous years, by vehicle type, at the year y
        self.model.x_maintain = Var(self.listVehicleType, self.listYearsVehicle, initialize=0, domain=NonNegativeReals)
        # Million tons of Copper and Aluminium demand for grid network, at the year y
        self.model.n = Var(self.listMetals, self.listYears, initialize=0, domain=NonNegativeReals)

        # Declaration of the relaxation variables

        # Potential requirement for additional ResLimit of metal m at the decade d
        self.model.Res_relax = Var(self.listMetals_knownRes, initialize=0, within=NonNegativeReals)
        # Potential requirement for additional mining production of metal m at the decade d
        self.model.Mining_relax = Var(self.listMetals, self.listDecades, initialize=0, within=NonNegativeReals)
        # In case the initial techno mix has a decrease for some technoRen
        self.model.CoherentMix_relax = Var(self.listRegions, self.listTechno_Ren, self.listDecades, initialize=0, within=NonNegativeReals)

    def modelObj(self):
        '''
        Objective of the optimisation model : to minimize the difference with IAM scenario
        '''

        epsilon = 1e-6

        self.model.objective = Objective(
            expr=
            sum(((((sum(self.TechnoMatrix[t].loc[T] * self.model.s[r,t,d] for t in self.listTechno))
                   -self.s0[r].loc[T][d])/(self.s0[r].loc[T][d]+epsilon))**2)
                for r in self.listRegions for T in self.listTechnoAgg for d in self.listDecades)


            + sum((((sum(self.TechnoMatrixEV[v].loc[V] * self.model.x[v, y] for v in self.listVehicle)
                     - self.x0.loc[V, y]) / (self.x0.loc[V, y]+epsilon)) ** 2)
                  for V in self.listVehicleAgg for y in self.listYearsVehicle)

            + sum(((self.model.n[m, y]*10**6 - self.Network_Demand.loc[m, y]) / (self.Network_Demand.loc[m, y] + epsilon)) ** 2 for m in
                  self.listMetals for y in self.listYears)

            + sum(self.M * (self.model.Res_relax[m] / self.MetalData[self.ResLimit].loc[m]) for m in self.listMetals_knownRes)
            + sum(self.M * (self.model.Mining_relax[m, d] / self.Prod[d].loc[m]) for m in self.listMetals for d in self.listDecades)
            + sum(self.M * (self.model.CoherentMix_relax[r, t_r, d]) for r in self.listRegions for t_r in self.listTechno_Ren for d in self.listDecades)
            #+ sum(self.M * (self.model.Stock_relax[vT, y]) ** 2 for vT in self.listVehicleType for y in self.listYearsVehicle)
            , sense=minimize)

    def CstrNetworkSubstitution(self):
        '''
        Constraint the model to meet electric grid requirement,
        with a potential mass substitution of aluminium to copper as 1:2
        '''

        self.model.ConstraintNetworkSubstitution = ConstraintList()
        for y in self.listYears:
            self.model.ConstraintNetworkSubstitution.add(
                (self.model.n['Aluminium', y] * 2 + self.model.n['Copper', y])*10**6
                == self.Network_Demand.loc['Aluminium', y] * 2 + self.Network_Demand.loc['Copper', y]
            )

    def CstrEnergyDemand(self):
        '''
        Constraint the model to meet the energy demand of the IAM scenario,
        for every year and every region
        '''

        # Creation of a list of k constraints, added using ".add" in ConstraintList
        # The power constraint must be respected for each year, and each region, on the total of technologies
        self.model.ConstraintEnergyDemand = ConstraintList()

        for r in self.listRegions:
            for d in self.listDecades:
                self.model.ConstraintEnergyDemand.add(
                    sum(sum(self.TechnoMatrix[t].loc[T] * self.model.s[r, t, d] for t in self.listTechno)
                        * self.CF_disag[r].loc[T][d] for T in self.listTechnoAgg)
                    >= sum(self.s0[r].loc[T][d] * self.CF_disag[r].loc[T][d] for T in self.listTechnoAgg)
                )

    def CstrVehicleDemand(self):
        '''
           Constraint the model to meet the predicted vehicle demand, for every year
        '''

        self.model.ConstraintVehicleGrowth = ConstraintList()
        for y in self.listYearsTot[1:]:
            for vT in self.listVehicleType:
                self.model.ConstraintVehicleGrowth.add(self.model.x_growth[vT, y] == self.VehicleType_Growth.loc[vT, y])

        self.model.ConstraintVehicleMaintain = ConstraintList()
        for y in range(1 + self.Lifetime_V, len(self.listYearsTot)):
            for vT in self.listVehicleType:
                self.model.ConstraintVehicleMaintain.add(self.model.x_maintain[vT, self.listYearsTot[y]]
                                                         == self.model.x_growth[vT, self.listYearsTot[y - self.Lifetime_V]])

        self.model.ConstraintVehicleSold = ConstraintList()
        for y in range(3, len(self.listYearsVehicle)):
            for vT in self.listVehicleType:
                self.model.ConstraintVehicleSold.add(
                    sum(self.model.x[v, self.listYearsVehicle[y]] for v in self.listVehicle if vT in v)
                    == self.model.x_growth[vT, self.listYearsVehicle[y]]
                    + self.model.x_maintain[vT, self.listYearsVehicle[y]]
                    )


    def CstrTechnoCoherence(self):
        '''
        Constraint model to keep the initial EV mix and IAM mix for 2020,
        and to keep the installed capacity of the previous decade during the following years
        '''

        # Creation of a list of constraint to add coherence in the optimised EV mix of 2020
        self.model.constraintEVCoherence2020 = ConstraintList()
        # The initial EV mix of 2020 cannot be changed
        for V in self.listVehicleAgg:
            for y in self.listYearsVehicle[:3]:
                self.model.constraintEVCoherence2020.add(sum(self.TechnoMatrixEV[v].loc[V] * self.model.x[v, y]
                                                             for v in self.listVehicle)
                                                         == self.x0.loc[V, y])

        # Creation of a list of constraint to add coherence in the optimised technological mix of 2020
        self.model.constraintMixCoherence2020 = ConstraintList()
        # The initial technological mix of 2020 cannot be changed
        for r in self.listRegions:
            for T in self.listTechnoAgg:
                self.model.constraintMixCoherence2020.add(sum(self.TechnoMatrix[t].loc[T] *self.model.s[r,t,'2020'] for t in self.listTechno)
                                                          ==self.s0[r].loc[T]['2020']
                                                          )

        # Creation of a list of constraint to add coherence in the evolution of the technological mix
        self.model.constraintMixCoherence = ConstraintList()
        # If a capacity is installed in d-1, it is still installed in the future cumulated installed capacity in d
        for r in self.listRegions:
            for t in self.listTechno_Ren:
                for d in range(1, len(self.listDecades)):
                    self.model.constraintMixCoherence.add(
                        self.model.s[r, t, self.listDecades[d]] + self.model.CoherentMix_relax[r, t, self.listDecades[d]] >= self.model.s[r, t, self.listDecades[d - 1]])

    def CstrCappedCC(self):
        '''
        Constrain the model to not increase fossil fuel technologies to limit climate change
        '''

        self.model.constraintCC = ConstraintList()

        for r in self.listRegions:
            for t in self.listTechno_Foss:
                for d in self.listDecades:
                    self.model.constraintCC.add(self.model.s[r, t, d] - self.s0[r].loc[t][d] == 0)

    def CstrMetal_Res(self):
        '''
        Constrain the model to limit cumulated metal demand in 2050 under the ResLimit restriction
        '''

        self.model.constraint_DemandCumulated = ConstraintList()

        for m in self.listMetals_knownRes:
            self.model.constraint_DemandCumulated.add(sum((self.model.s[r,t,self.listDecades[d]]-self.model.s[r,t,self.listDecades[d-1]])
                                          *((self.MetalIntensity_doc[self.listDecades[d]][t].loc[m]+self.MetalIntensity_doc[self.listDecades[d-1]][t].loc[m])/2)
                                          /self.MetalData['RR (%) Res'].loc[m] for t in self.listTechno_Ren for r in self.listRegions for d in range(1, len(self.listDecades)))
                                                      # Metal demand to create the vehicle stock with recycling between 2020 and 2030
                                                      + sum(
                sum(self.model.x[v, self.listYearsVehicle[y]] * self.MI_Vehicle.loc[m][v] for v in self.listVehicle)
                - self.model.x_growth['ICEV', self.listYearsTot[y + 1]] * self.Recycling_Rate.loc[m,self.listYearsVehicle[y]] * self.MI_Vehicle.loc[m]['ICEV'] for y in
                range(2, 12))/self.MetalData['RR (%) Res'].loc[m]
                                                      # Metal demand to create the vehicle stock with recycling between 2030 and 2040
                                                      + sum(
                (self.model.x[v, self.listYearsVehicle[y]] - self.model.x[v, self.listYearsTot[y + 1]] * self.Recycling_Rate.loc[m,self.listYearsVehicle[y]]) * self.MI_Vehicle.loc[m][v] /
                self.MetalData['RR (%) Res'].loc[m] for v in self.listVehicle for y in range(12, len(self.listYearsVehicle)))

                                                      + sum(self.Storage_Demand.loc[m][y] / self.MetalData['RR (%) Res'].loc[m]
                                                            for y in self.listYears)
                                                      # Metal demand cumulated for network grid
                                                      + sum(
                self.model.n[m, y]*10**6 / self.MetalData['RR (%) Res'].loc[m] for y in self.listYears)
                                                      + sum(
                (self.OSD[self.listDecades[d]].loc[m] + self.OSD[self.listDecades[d - 1]].loc[m]) * 10 / 2 for d in
                range(1, len(self.listDecades)))
                                                      - self.model.Res_relax[m]
                                                      <= self.MetalData[self.ResLimit].loc[m])

        # model.constraint_DemandCumulated.pprint()

    def CstrMetal_Mining(self):
        '''
        Constrain the model to limit annual metal demand under the estimated metal production
        or to limit metal demand for energy under an Alpha percentage of the demand for
        the rest of the transition and economy
        '''

        def renewable_demand(m, d):
            return sum(
                self.MetalIntensity_doc[self.listDecades[d]][t].loc[m]
                * (self.model.s[r, t, self.listDecades[d]] - self.model.s[r, t, self.listDecades[d - 1]])
                / self.MetalData['RR (%) Prod'].loc[m]
                for t in self.listTechno_Ren
                for r in self.listRegions
            ) / 10

        def network_demand(m, d):
            return (sum(
                self.model.n[m, self.listYearsTot[5 + i + 10 * d]]*10**6 for i in self.list_i)
                    / (self.MetalData['RR (%) Prod'].loc[m] * 10))

        def vehicle_demand(m, d):
            if d == 1:
                return sum(
                    (
                        sum(
                            self.model.x[v, self.listYearsTot[5 + i + 10 * d]] * self.MI_Vehicle.loc[m, v]
                            for v in self.listVehicle
                        )
                        - self.model.x_growth['ICEV', self.listYearsTot[5 + i + 10 * d - self.Lifetime_V]]
                        * self.MI_Vehicle.loc[m, 'ICEV'] * self.Recycling_Rate.loc[m,self.listYearsTot[5 + i + 10 * d]]
                    )
                    / self.MetalData['RR (%) Prod'].loc[m]
                    for i in self.list_i
                ) / 10
            else:
                return sum(
                    (
                            sum(
                                self.model.x[v, self.listYearsTot[5 + i + 10 * d]]
                                - self.model.x[v, self.listYearsTot[5 + i + 10 * d - self.Lifetime_V]] * self.Recycling_Rate.loc[m,self.listYearsTot[5 + i + 10 * d]]
                                for i in self.list_i
                            )
                            / 10
                    )
                    * self.MI_Vehicle.loc[m][v]
                    / self.MetalData['RR (%) Prod'].loc[m]
                    for v in self.listVehicle
                )

        def max_production(m, d):
            return max(
                self.Prod[self.listDecades[d]].loc[m]
                - self.Storage_Demand[self.listDecades[d]].loc[m] / self.MetalData['RR (%) Prod'].loc[m]
                - self.OSD[self.listDecades[d]].loc[m],
                self.Alpha / 100
                * (
                        self.OSD[self.listDecades[d]].loc[m]
                        + self.Storage_Demand[self.listDecades[d]].loc[m] / self.MetalData['RR (%) Prod'].loc[m]
                ),
            )

        self.model.Constraint_DemandByYear = ConstraintList()

        for m in self.listMetals:
            for d in range(1, len(self.listDecades)):
                demand_renewable = renewable_demand(m, d)
                demand_network = network_demand(m, d)
                demand_vehicle = vehicle_demand(m, d)
                production_capacity = max_production(m, d)

                # Add constraint
                self.model.Constraint_DemandByYear.add(
                    demand_renewable + demand_network + demand_vehicle
                    <= production_capacity + self.model.Mining_relax[m, self.listDecades[d]]
                )

    def CstrNewTechnoEV(self):

        '''
        Constraint to model the limit of the penetration of the new sodium-ion battery
        Supposition of a price advantage for Na-ion battery from 2035
        '''

        self.model.Constraint_NewTechnoEV = ConstraintList()
        # Before 2035, cannot add more than 10% of Na-ion already presumed
        for y in self.listYearsVehicle[:23]:
            self.model.Constraint_NewTechnoEV.add(sum(
                sum(self.TechnoMatrixEV[v].loc[V] * self.model.x[v, y] for v in self.listVehicle)
                for V in self.listVehicleAgg if "Na-ion" in V)
                <= 1.1 * sum(self.x0.loc[V, y] for V in self.listVehicleAgg if "Na-ion" in V))

    def CstrBiomassAvail(self):
        '''
        Constraint the model to limit biomass power capacity under
        a potential of future biomass availability
        '''

        self.model.Constraint_Biomass = ConstraintList()

        for d in self.listDecades:
            self.model.Constraint_Biomass.add(
                sum(self.model.s[r, 'Biomass', d] for r in self.listRegions)
                <= self.Biomass_availability.loc['Business as usual'][d]
            )

    def CstrHydroAvail(self):
        '''
        Constraint the model to limit hydraulic power capacity under
        a potential of future hydraulic availability by region
        '''

        self.model.Constraint_Hydro = ConstraintList()

        for r in self.listRegions[:3]:  # Remaining hydro availability for ASIA, MAF and LAM
            self.model.Constraint_Hydro.add(
                self.model.s[r, 'Hydro', '2050'] - self.model.s[r, 'Hydro', '2020'] <= self.Hydro_availability[r].loc[
                    'Technical -Ecological']
            )

        # Remaining hydro availability for combined Europe and North-Central America
        self.model.Constraint_Hydro.add(
            sum(self.model.s[r, 'Hydro', '2050'] - self.model.s[r, 'Hydro', '2020'] for r in self.listRegions[3:])
            <= self.Hydro_availability['Europe'].loc['Technical -Ecological'] + self.Hydro_availability['North America'].loc[
                'Technical -Ecological']
        )

    def resOpti(self):
        '''
        Results of the optimisation
        '''
        opt.solve(self.model, tee=True)

    def Save_Opti_Characteristics(self):
        '''
        Saves the optimisation characteristics : objective and relaxation variables
        for the ResLimit and Mining production.
        '''

        # Objective
        self.Obj = pd.Series(self.model.objective())

        # Relaxation variable for the additional resources needed for each metal
        self.Relax_Var_Res = pd.Series(self.model.Res_relax.get_values())
        # Create an empty dataFrame to stock the variable
        Relax_Var_Mining = pd.DataFrame(index=self.listMetals, columns=self.listDecades)
        for m in self.listMetals:
            for y in self.listDecades:
                # Stock the variable for the consumption of each metals by year in 2050 in a dataFrame
                Relax_Var_Mining.loc[m, y] = self.model.Mining_relax.get_values()[(m, y)]
        self.Relax_Var_Mining = Relax_Var_Mining

        '''
        # Create an empty dataFrame to stock the variable
        Relax_Var_Stock = pd.DataFrame(index=self.listVehicleType, columns=self.listDecades)
        for vT in self.listVehicleType:
            for y in self.listDecades:
                # Stock the variable for the consumption of each metals by year in 2050 in a dataFrame
                Relax_Var_Stock.loc[vT, y] = self.model.Stock_relax.get_values()[(vT, y)]
        self.Relax_Var_Stock = Relax_Var_Stock
        '''

        # Save the Relaxation variables and the Objective result
        # Generate the full path for the Excel file for this scenario and model
        folderRelVar = self.Res_folder + '/Relax Var'
        if not os.path.exists(folderRelVar):
            os.makedirs(folderRelVar)
        excel_path = folderRelVar +'/RelVar_' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path)
        self.Relax_Var_Res.to_excel(excel, sheet_name='RelVarRes')
        self.Relax_Var_Mining.to_excel(excel, sheet_name='RelVarMining')
        #self.Relax_Var_Stock.to_excel(excel, sheet_name='RelVarStock')
        self.Obj.to_excel(excel, sheet_name='Obj')
        excel.close()

    def CalculResults (self):
        '''
        Calcul results for initial and optimised power mix, cumulated metal demand and metal demand by year
        '''

        # Stock optimisation results

        # Creation of a dictionary res with the results
        if self.ModelisationType == 'Init':
            # PowerCap initial results
            self.res_dic = self.s0
            # Network initial results
            self.res_Network = self.Network_Demand
            # EV initial results
            self.res_EV = self.x0
        if self.ModelisationType == 'Opti':
            # PowerCap optimised results
            result = self.model.s.get_values()
            res_dic = {}
            for key, value in result.items():
                r, t, d = key
                if r not in res_dic:
                    res_dic[r] = pd.DataFrame()
                # Update the value in the DataFrame at the specified index and column
                res_dic[r].loc[t, d] = value
            self.res_dic = res_dic

            # Network optimised results
            Network_opti_dic = self.model.n.get_values()
            res_Network = pd.DataFrame(index=self.listMetals, columns=self.listYears)
            for m in self.listMetals:
                for y in self.listYears:
                    res_Network.loc[m][y] = Network_opti_dic[m, y]*10**6
            self.res_Network = res_Network

            # EV optimised results
            EV_opti_dic = self.model.x.get_values()
            res_EV = pd.DataFrame(columns= self.listYearsVehicle, index= self.listVehicle)
            for v in self.listVehicle:
                for d in self.listYearsVehicle:
                    res_EV.loc[v][d] = EV_opti_dic[v, d]
            self.res_EV = res_EV

        # Calcul cumulated demand
        # Create a list with metals needed for the energy sources
        EnergySectorDemand = pd.Series()
        # Create a list with metals needed for the storage
        NetworkDemand = pd.Series()
        # Create a list with metals needed for the electric vehicles
        EVSectorDemand = pd.Series()
        # Create a list with metals needed for storage and electric grid
        StorageDemand = pd.Series()
        # Create a list with metals needed for the other sector of society
        OtherSectorDemand = pd.Series()
        # Add datas for each Series
        for m in self.listMetals:
            EnergySectorDemand[m] = sum(
                (self.res_dic[r].loc[t][self.listDecades[d]] - self.res_dic[r].loc[t][self.listDecades[d - 1]])
                * ((self.MetalIntensity_doc[self.listDecades[d]][t].loc[m] +
                    self.MetalIntensity_doc[self.listDecades[d - 1]][t].loc[m]) / 2)
                / self.MetalData['RR (%) Res'].loc[m] for t in self.listTechno_Ren for r in self.listRegions for d
                in range(1, len(self.listDecades)))
            NetworkDemand[m] = sum(self.res_Network.loc[m][y] / self.MetalData['RR (%) Res'].loc[m] for y in self.listYears)
            EVSectorDemand[m] = (sum(sum(self.res_EV.loc[v, self.listYearsVehicle[y]]* self.MI_Vehicle.loc[m][v] for v in self.listVehicle)
                                         - self.VehicleType_Growth.loc['ICEV', self.listYearsTot[y + 1]] * self.Recycling_Rate.loc[m,self.listYearsVehicle[y]] * self.MI_Vehicle.loc[m]['ICEV']
                                         for y in range(2, 12))/ self.MetalData['RR (%) Res'].loc[m]
                                 + sum((self.res_EV.loc[v, self.listYearsVehicle[y]]
                                        - self.res_EV.loc[v, self.listYearsTot[y + 1]] * self.Recycling_Rate.loc[m,self.listYearsVehicle[y]])
                                       * self.MI_Vehicle.loc[m][v] / self.MetalData['RR (%) Res'].loc[m]
                                       for v in self.listVehicle for y in range(12, len(self.listYearsVehicle))))

            StorageDemand[m] = sum(
                self.Storage_Demand.loc[m][y] / self.MetalData['RR (%) Res'].loc[m] for y in self.listYears)
            OtherSectorDemand[m] = + self.OSD['2020'].loc[m] + sum(
                (self.OSD[self.listDecades[d]].loc[m] + self.OSD[self.listDecades[d - 1]].loc[m]) * 10 / 2 for d in
                range(1, len(self.listDecades)))
            # Concatenate the Series along axis 1 (columns), specifying keys for the resulting DataFrame
        DemandCumulated = pd.concat([OtherSectorDemand, StorageDemand, NetworkDemand, EVSectorDemand, EnergySectorDemand], axis=1,
                                    keys=['Other Sector Demand', 'Storage Demand', 'Network Demand', 'EV Sector Demand','Energy Sector Demand'])
        self.DemandCumulated = DemandCumulated

        # Calcul demand by year
        # Create a dF for newly installed capacity by decades
        newCapacityD = {}
        for r in self.listRegions:
            newCapacityD[r] = pd.DataFrame(0.0, index=self.listTechno, columns=self.listDecades[1:])
            for t in self.listTechno:
                for d in range(1, len(self.listDecades)):
                    if self.res_dic[r].loc[t][self.listDecades[d]] - self.res_dic[r].loc[t][
                        self.listDecades[d - 1]] > 0:  # Changes the value only for an addition of capacity
                        newCapacityD[r][self.listDecades[d]].loc[t] = self.res_dic[r].loc[t][self.listDecades[d]] - \
                                                                      self.res_dic[r].loc[t][
                                                                          self.listDecades[d - 1]]
        self.newCapacityD = newCapacityD

        # Create a dF with metals needed by sector
        EnergySectorDemandy = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)
        NetworkDemandy = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)
        EVSectorDemand = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)

        # Add initial data for 2020
        EnergySectorDemandy['2020'] = self.PowerMetalDemand_y['2020']
        NetworkDemandy['2020'] = self.Network_Demand['2020']
        EVSectorDemand['2020'] = sum(self.x0.loc[v, '2020'] * self.MI_VehicleAgg[v] for v in self.listVehicleAgg)

        # Add data for next years
        for m in self.listMetals:
            for d in range(1, len(self.listDecades)):

                # Add metals need for new installed capacity between the actual decade and the previous one
                EnergySectorDemandy[self.listDecades[d]].loc[m] = sum(
                    self.MetalIntensity_doc[self.listDecades[d]][t].loc[m]
                    * newCapacityD[r][self.listDecades[d]].loc[t] / 10
                    for t in self.listTechno_Ren for r in self.listRegions)

                NetworkDemandy[self.listDecades[d]].loc[m] = (sum(
                    self.res_Network.loc[m][self.listYearsTot[5 + i + 10 * d]] for i in self.list_i)
                                                               /(self.MetalData['RR (%) Prod'].loc[m]*10))

                EVSectorDemand[self.listDecades[d]].loc[m] = sum(
                    self.res_EV.loc[v, self.listYearsTot[5 + i + 10 * d]] * self.MI_Vehicle.loc[m, v]
                    for v in self.listVehicle for i in self.list_i
                    )/(10*self.MetalData['RR (%) Prod'].loc[m])
        self.EVSectorDemand = EVSectorDemand

        # Create a three dimensions dataFrame to add all datas of consumption by decade
        DemandByYear_df = {}
        for m in self.listMetals:
            # Index for different datasSets, columns for each decades
            DemandByYear_df[m] = pd.DataFrame(
                index=['Other Sector Demand', 'Storage Demand', 'Network Demand', 'EV Sector Demand','Energy Sector Demand'],columns=self.listDecades)
            for i in DemandByYear_df[m].index:
                # Add the values for each dataSets
                DemandByYear_df[m].loc['Energy Sector Demand'] = EnergySectorDemandy.loc[m] / \
                                                                   self.MetalData['RR (%) Prod'].loc[m]
                DemandByYear_df[m].loc['Other Sector Demand'] = self.OSD.loc[m]
                DemandByYear_df[m].loc['Network Demand'] = NetworkDemandy.loc[m] / \
                                                           self.MetalData['RR (%) Prod'].loc[m]
                DemandByYear_df[m].loc['Storage Demand'] = self.Storage_Demand.loc[m] / self.MetalData['RR (%) Prod'].loc[m]
                DemandByYear_df[m].loc['EV Sector Demand'] = self.EVSectorDemand.loc[m]
        self.DemandByYear_df = DemandByYear_df

        # Create a dF for secondary production for EV
        # Create a dF for secondary production for EV
        Secondary_Prod_df = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)
        for m in self.listMetals:
            for d in range(1, len(self.listDecades)):
                if d==1:
                    Secondary_Prod_df.loc[m, self.listDecades[d]] = (sum(
                        self.VehicleType_Growth.loc['ICEV', self.listYearsTot[5 + i + 10 * d - self.Lifetime_V]]
                        * self.Recycling_Rate.loc[m,self.listYearsTot[5 + i + 10 * d]] for i in self.list_i)
                                                                      * self.MI_Vehicle.loc[m]['ICEV']
                                                                     /(10*self.MetalData['RR (%) Prod'].loc[m]))
                elif d>1:
                    Secondary_Prod_df.loc[m, self.listDecades[d]] = (sum(self.res_EV.loc[v, self.listYearsTot[5 + i + 10 * d - self.Lifetime_V]]
                                * self.MI_Vehicle.loc[m, v] * self.Recycling_Rate.loc[m,self.listYearsTot[5 + i + 10 * d]] for v in self.listVehicle for i in self.list_i)
                                                                     /(10 * self.MetalData['RR (%) Prod'].loc[m]))

        self.Secondary_Prod_df = Secondary_Prod_df


    def SaveResults(self):
        '''
        Save results of power capacity, cumulated and annual metal demand in the path chosen for results
        '''

        if not os.path.exists(self.Res_folder):
            os.makedirs(self.Res_folder)

        # Save the power capacity of the model
        for r in self.listRegions:
            table = self.res_dic[r]

            # Create a file to put the excels of the regions in it
            folderCap = self.Res_folder+'/Power Capacity/PowCap' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario
            if not os.path.exists(folderCap):
                os.makedirs(folderCap)
            # Generate the full path for the Excel file for this metal
            excel_pathCap = folderCap + '/PowCap_' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario + '_' + r + '.xlsx'
            # Create the Excel file and write the data
            excel = pd.ExcelWriter(excel_pathCap)
            table.to_excel(excel, sheet_name='PowCap')
            excel.close()

        # Save the EV of the model
        # Generate the full path for the Excel file for this scenario and model
        folderEV = self.Res_folder + '/Vehicle'
        if not os.path.exists(folderEV):
            os.makedirs(folderEV)
        excel_path = folderEV + '/EV_' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path)
        self.res_EV.to_excel(excel, sheet_name='EV')
        excel.close()

        # Save the secondary production of the model
        # Generate the full path for the Excel file for this scenario and model
        folderRecycling = self.Res_folder + '/SecondaryProdEV'
        if not os.path.exists(folderRecycling):
            os.makedirs(folderRecycling)
        excel_path = folderRecycling + '/SecondaryProdEV_' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path)
        self.Secondary_Prod_df.to_excel(excel, sheet_name='SecondaryProdfromEV')
        excel.close()

        # Save the Cumulated Demand in the computer
        # Generate the full path for the Excel file for this scenario and model
        folderDemCum = self.Res_folder + '/Metal Demand Cumulated'
        if not os.path.exists(folderDemCum):
            os.makedirs(folderDemCum)
        excel_path = folderDemCum + '/DemCum_' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path)
        self.DemandCumulated.to_excel(excel, sheet_name='DemCum')
        excel.close()

        # Save DemandByYear Results
        for m in self.listMetals:
            table = self.DemandByYear_df[m]
            # Create a file to put the excels of the metals in it
            folderIAM = self.Res_folder + '/Metal Demand ByYear/DemByYear' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario
            if not os.path.exists(folderIAM):
                os.makedirs(folderIAM)
            # Generate the full path for the Excel file for this metal
            excel_path = folderIAM + '/DemByYear_' + self.ModelisationType + '_' + self.model_s0 + '_' + self.scenario + '_' + m + '.xlsx'
            # Create the Excel file and write the data
            excel = pd.ExcelWriter(excel_path)
            table.to_excel(excel, sheet_name='DemByYear')
            excel.close()

    def Save_OSD_Prod (self):

        '''
        Save datas for other sector demand and other transition demand by year and cumulated by decade
        Save estimated production by year
        '''

        Res_Folder_BroaderEconomy = self.result_path+ 'BroaderEconomy'
        if not os.path.exists(Res_Folder_BroaderEconomy):
            os.makedirs(Res_Folder_BroaderEconomy)

        # Save OSD
        # Generate the full path for the Excel file for this scenario and model
        OSD_Folder = Res_Folder_BroaderEconomy + '/OSD'
        if not os.path.exists(OSD_Folder):
            os.makedirs(OSD_Folder)
        excel_path_OSD = OSD_Folder + '/OSD_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path_OSD)
        self.OSD.to_excel(excel, sheet_name='OSD')
        excel.close()

    def SaveEffortBySociety(self):
        '''
        Saves an Excel file, that shows for each decade and each metal, if the mining constraint have been
        respected because total metal demand is under production ('Prod')
        or if the demand of metal for energy is less than Alpha percent of the rest of economy ('Society')
        It would then suppose a requirement of effort from society to decrease consumption
        '''

        ConstraintDemandByYear_file = self.result_path + 'BroaderEconomy/EffortBySociety' + '_' + self.model_s0 + '_' + self.scenario + '.xlsx'

        if not os.path.exists(ConstraintDemandByYear_file):

            ConstraintDemandByYear = pd.DataFrame(index=self.listMetals, columns=self.listDecades[1:])

            for m in self.listMetals:
                for d in self.listDecades[1:]:
                    # If total metal demand is under the estimated production
                    if (self.Prod[d].loc[m] - self.Storage_Demand[d].loc[m] / self.MetalData['RR (%) Prod'].loc[m] - self.OSD[d].loc[m]) >= (
                            self.Alpha / 100 * self.OSD[d].loc[m]):
                        ConstraintDemandByYear[d].loc[m] = 'Prod'
                    # If there is a need for an effort from society
                    else:
                        ConstraintDemandByYear[d].loc[m] = 'Society'

            # Function to apply red color
            def color_red(val):
                color = 'red' if val == 'Society' else 'black'
                return f'color: {color}'

            # Apply conditional style
            Styled_ConstraintDemandByYear = ConstraintDemandByYear.style.map(color_red)
            Styled_ConstraintDemandByYear

            # Create the Excel file and write the data
            excel = pd.ExcelWriter(ConstraintDemandByYear_file)
            Styled_ConstraintDemandByYear.to_excel(excel)
            excel.close()