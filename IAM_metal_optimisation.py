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

# Choose solver
opt = SolverFactory('ipopt', executable="C:/Users/Penel/Ipopt-3.11.1-win64-intel13.1/bin/ipopt.exe")

# Choose the solver for the optimization
opt = SolverFactory('ipopt', executable="C:/Users/Penel/Ipopt-3.11.1-win64-intel13.1/bin/ipopt.exe")
opt.options['max_iter'] = 10000  # define maximum number of optimization iterations

class IAM_Metal_Optimisation :
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

        :param folder_path: [string] path to external datas used in the code
        :param result_path: [string] path to stock results
        :param model_s0: [string] Integrated Assessment Model (IAM) that wants to be studied
        (AIM..CGE, GCAM4, IMAGE, MESSAGE-GLOBIOM, REMIND-MAGPIE, WITCH-GLOBIOM)
        :param scenario: [string] scenario of Shared Socioeconomic Pathway and  Representative Concentration Pathway
        (SSP-RCP) that wants to be studied
        :param ModelisationType: [string] 'Init' to study initial IAM power mix metal demand
                                 'Opti' to generate an optimisation of IAM power mix
        :param Alpha: [int] percentage value of acceptable share of metal energy demand in comparison to OSD and OTD
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
        # Acceptable share of metal energy demand in comparison to OSD and OTD for the optimisation
        self.Alpha = Alpha

        # Penalisation by M of the relaxation variable in the objective function
        self.M = 10**15
        # Folder path for results according to the modelisation type chosen
        self.Res_folder = self.result_path + '/'+ self.ModelisationType

        # Definition of the scope of this study in terms of models, scenarios, regions, decades, metals and techno
        # External excel defining the scope of the study
        studyScope = self.folder_path + 'Scope of the study.xlsx'

        with pd.ExcelFile(studyScope) as file:
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
            self.listTechnoFlex= pd.read_excel(file, 'technologies flex', index_col=0).squeeze().tolist()

        # Create a list of Years from the first decade to the last one
        listYears = [str(year) for year in range(int(self.listDecades[0]), int(self.listDecades[-1]) + 1)]
        self.listYears = listYears

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
        self.MetalDatas = pd.read_excel(self.folder_path + 'Metal datas.xlsx', index_col=0)

        # Creation of a list of metals with known ResLimit (exclusion of metal with unknown reserve or resources)
        self.listMetals_knownRes = self.listMetals.copy() # Create a copy of all the metals studied
        for metal in self.listMetals:
            if pd.isna(self.MetalDatas[self.ResLimit].loc[metal]): # If the data for ResLimit is missing
                self.listMetals_knownRes.remove(metal)  # Remove the metal

        # Import datas for total 2050 biomass availability [GW], according to various scenario
        self.Biomass_availability = pd.read_excel(self.folder_path + 'Future techno availability.xlsx', sheet_name = 'Biomass' , index_col=0)

        # Import datas for remaining hydroelectricity availability, for various two economic and ecologic scenarios, for each regions [GW]
        self.Hydro_availability = pd.read_excel(self.folder_path + 'Future techno availability.xlsx', sheet_name = 'Hydro',index_col=0)

        # Import CF from the literature
        self.CF_from_litt = pd.read_excel(self.folder_path + 'Capacity Factor IAM/1_CF_from_litt.xlsx', index_col=0)

        # Import matrix for corresponding CF of aggregated and disaggregated sub-technologies
        Matrix_CF_Init_disag = pd.read_excel(self.folder_path + 'Capacity Factor IAM/0_Matrix_CF_disagg.xlsx',
                                             index_col=0)
        self.Matrix_CF_disag = Matrix_CF_Init_disag.sort_values(by='Techno ') # Rearrange matrix in the alphabetic order of technologies

        # Import metal consumption by decade for EV, storage in [t/decade] (from IEA)
        self.OTDdec = pd.read_excel(self.folder_path + 'Datas IEA/IEA_ByDecade_Demand_Metal.xlsx')

        # Import metal consumption by year for EV, storage in [t/decade] (from IEA)
        self.OTDyear = pd.read_excel(self.folder_path + 'Datas IEA/IEA_ByYear_Demand_Metal.xlsx')

        # Tab of future metal production according to different scenarios from 2020 to 2050 [t/yr]
        #self.Prod = pd.read_excel(self.folder_path + 'Future metal prod.xlsx', index_col=0)

        # Importation of GDP projections for a specific model, scenario ssp-rcp, at a specific year (from IAM)
        self.GDP_folder = self.folder_path + 'GDP IAM/GDP ' + self.model_s0 + '_' + self.scenario

        # Folder used to stock initial Power capacity projections for specific model_s0 and scenario, for various regions
        self.folderIAM = self.folder_path + 'Power Capacity IAM/Dossier s0 ' + self.model_s0 + '_' + self.scenario

        # Import the matrix of correspondence between Techno and Techno Init
        self.TechnoMatrix = pd.read_excel(self.folder_path + 'Techno Matrix.xlsx', index_col=0)

        # Metal intensity of technologies 
        self.MetalIntensity_doc = MetalIntensity_doc

    # Calculate Metal Intensity (MI) evolution in time
    def MI(self):
        '''
        Calculation of the evolution of the metal intensity of energy sources technology, according
        to the scenario SSP chosen.
        :return: Metal intensity in tons of metal per GW of energy sources, with evolution in time [t/GW]
        '''

        # Import table of metal intensities for specific sub-technologies [t/GW] according to the MI scenario chosen
        MetalIntensity = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx', sheet_name = 'MI base',index_col=0)
        # Import the table of metal intensities for specific sub-technologies
        self.MetalIntensity_Agg = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx', sheet_name='Aggregated MI 2010',
                                           index_col=0)

        # Creation of a dictionary for the Metal Intensity of technologies over time
        MetalIntensity_doc = {}
        # Take the initial data for the first decade metal intensity
        MetalIntensity_doc[self.listDecades[0]] = MetalIntensity

        # Tab of different scenario of metal intensity reduction
        ReductionScenario = pd.read_excel(self.folder_path + 'Metal Intensity Tech.xlsx', sheet_name='IM reduction scenario',
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

    def CF(self):
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
                for techno in self.listTechnoFlex:
                    # if the techno is a sub-techno of solar CSP
                    if "CSP" in techno:
                        # replace the CF by the one from the literature
                        CF_disag[region][years][techno] = self.CF_from_litt["CF from literature"].loc["CSP"] * 100

        # Change the CF of solar to be constant in time, for the decades studied
        for region in self.listRegions:
            for years in self.listDecades[1:]: # exclude the first decade
                for techno in self.listTechnoFlex:
                    # if the techno is a sub-techno of solar PV
                    if "Sol" in techno:
                        # replace the CF by the one from IAM at the first decade
                        CF_disag[region][years][techno] = CF_disag[region][self.listDecades[0]][techno]
        # Final regionalized CF for the disaggregated sub-technologies
        self.CF_disag = CF_disag

    def IAM_power(self):
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
            self.listTechnoIAM = self.listTechnoFlex
        if self.ModelisationType =='Init':
            self.listTechnoIAM = self.listTechno

        # Replace in the initial scenario s0 the 0 power capacity by a very little value, to avoid the 0 division
        for r in self.listRegions:
            for t in self.listTechnoIAM:
                for d in self.listDecades:
                    # if the power capacity is 0
                    if self.s0[r].loc[t, d] == 0:
                        # replace it by a very little value, in comparison to the total capacity
                        self.s0[r].loc[t, d] = 1


    def OTD(self):
        '''
        Calculate the demand in metal for the rest of the transition : electric vehicle, grid, storage
        Based on different International Energy Agency (IEA) scenarios
        '''

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



        # Select only the data of the scenario studied, thanks to the column scenario
        self.OTDdec = self.OTDdec[self.OTDdec["Scenario"] == scenarioIEA]
        # Drop the columns for technology and scenario, to only have the metal information
        self.OTDdec = self.OTDdec.drop(columns=['Technology', 'Scenario', 'Metal'])
        # Group IEA data by metals
        self.OTDdec = self.OTDdec.groupby('Metal [t/decade]')
        # Sum for every IEA studied technology, by year, for each metal
        self.OTDdec = self.OTDdec.sum()

        self.OTDd = pd.DataFrame(index=self.listMetals, columns=self.listDecades)

        for d in self.listDecades:
            for m in self.listMetals:
                if m in self.OTDdec.index:
                    self.OTDd[d].loc[m] = self.OTDdec[d].loc[m]
                else:
                    self.OTDd[d].loc[m] = 0

        # Importation of IEA data, in t/decade, non cumulated
        self.OTDyear = pd.read_excel(self.folder_path + 'Datas IEA/IEA_ByYear_Demand_Metal.xlsx')
        # Select only the data of the scenario studied, thanks to the column scenario
        self.OTDyear = self.OTDyear[self.OTDyear["Scenario"] == scenarioIEA]
        # Drop the columns for technology and scenario, to only have the metal information
        self.OTDyear = self.OTDyear.drop(columns=['Technology [t/yr]', 'Detail', 'Scenario'])
        # Group IEA data by metals
        self.OTDyear = self.OTDyear.groupby('Metal')
        # Sum for every IEA studied technology, by year, for each metal
        self.OTDyear = self.OTDyear.sum()
        self.OTDy = pd.DataFrame(index=self.listMetals, columns=self.listDecades)

        for d in self.listDecades:
            for m in self.listMetals:
                if m in self.OTDyear.index:
                    self.OTDy[d].loc[m] = self.OTDyear[d].loc[m]
                else:
                    self.OTDy[d].loc[m] = 0

    def FutureProd(self):

        '''
        Estimation of future metal production in tons by year, until 2050
        '''

        # Initialisation of a dataframe with metals in lines and years in columns
        Prod = pd.DataFrame(index=self.listMetals, columns=self.listYears)

        # Initialise the datas from 2020 with actual known 2020 metal production
        Prod["2020"] = self.MetalDatas['Production 2020 [t/yr]']
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
        Prod = pd.DataFrame({'2020': self.MetalDatas['Production 2020 [t/yr]'], '2030': ProdFinalbyYear['2030'],
                             '2040': ProdFinalbyYear['2040'], '2050': ProdFinalbyYear['2050']})

        # Correction of future prod with IEA High production case projections for Cobalt, Copper, Nickel
        Prod_IEA = pd.read_excel(self.folder_path + 'Future metal production.xlsx', sheet_name='IEA', index_col=0)

        for m in Prod_IEA.index:
            for d in self.listDecades:
                Prod[d].loc[m] = Prod_IEA[d].loc[m]

        self.Prod = Prod

    def OSD(self):
        '''
        Estimation of the metal demand of the rest of the economy
        Based on literature and Gross Domestic Product (GDP) growth from IAM and SSP-RCP
        '''

        # Creation of the GDP dictionary
        excelGDP = [f for f in os.listdir(self.GDP_folder) if
                    f.endswith('xlsx')]  # Reads the folder with excel files of CF from IAM dataset

        # Creation of a dictionary of dataFrame to stock datas from Excel files
        self.GDP = {}

        for region_name, file in zip(self.listRegions, excelGDP):
            fileGDPbyRegion = os.path.join(self.GDP_folder, file)  # copies datas from files in dataFrames
            self.GDP[region_name] = pd.read_excel(fileGDPbyRegion,
                                                  index_col=0)  # copies the dataFrames in dic, with first col as index

        # New capacity installed between 2010 and 2020, by techno
        new2020Capacity_d = pd.DataFrame(0.0, index=self.listTechnoFlex, columns=['2020'])

        for t in self.listTechnoIAM:
            for r in self.listRegions:
                # Only account if there is an increase in capacity
                if self.s0[r].loc[t, '2020'] - self.s0[r].loc[t, '2010'] > 0:
                    # The new installed capacity in the decade between 2020 and 2010, devided by ten, to have the 2020 yearly value
                    new2020Capacity_d['2020'].loc[t] = sum(
                        self.s0[r]['2020'].loc[t] - self.s0[r]['2010'].loc[t] for r in self.listRegions)

        # Metal demand for the new installed capacity at the year 2020
        PowerMetalDemand_y = pd.DataFrame(0.0, index=self.listMetals, columns=['2020'])

        for m in self.listMetals:
            for t in self.listTechnoFlex:
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

        # Initialisation of a dataframe for other sector demand in metals, with metals in lines and years in columns
        OSDy = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)

        # Initialise the datas from 2020 with actual known metal production at the year 2020 - demand in metal for the transition [t/y]
        for m in self.listMetals:  # Total metal prod in 2020 - Metals for EV, grid and storage - # Metals for power infrastructure installation in 2020
            InitialOTD = (self.Prod['2020'].loc[m] - self.OTDy['2020'].loc[m] / self.MetalDatas['RR (%) Prod'].loc[m]
                          -self.PowerMetalDemand_y['2020'].loc[m] / self.MetalDatas['RR (%) Prod'].loc[m])
            if InitialOTD > 0:
                OSDy["2020"].loc[m] = InitialOTD

            # Initialise a final demand dF to stock future demand in metals, with GDP increase
        # Demand_byYear = {year: Demand[year].copy() for year in listDecades}

        for m in self.listMetals:
            # Loop through list of years, excluding the first year 2020
            for d in range(1, len(self.listDecades)):
                # Create a variable for the demand at the previous year
                previous_demand = OSDy[self.listDecades[d - 1]].loc[m]
                # Calculate the World GDP increase with a sum for every region
                GDP_growth = sum(self.GDP[r][self.listDecades[d]].loc['GDP..PPP'] for r in self.listRegions) / sum(
                    self.GDP[r][self.listDecades[d - 1]].loc['GDP..PPP'] for r in self.listRegions)
                # Calculate Final demand with the previous one and GDP increase
                OSDy[self.listDecades[d]].loc[m] = previous_demand * GDP_growth

        OSD_litt = pd.read_excel(self.folder_path + 'OSD litt.xlsx', index_col=0)

        # Correction with available data for the metal m in the literature of future OSD
        for m in OSD_litt.index:
            for d in self.listDecades:
                OSDy[d].loc[m] = OSD_litt[d].loc[m]

        OSD_litt_byScenario = pd.read_excel(self.folder_path + 'OSD litt.xlsx', sheet_name='OSD by scenario', index_col=0)

        # Correction with projections from IEA, matching scenario
        OSD_ScenarioIEA = OSD_litt_byScenario.loc[OSD_litt_byScenario['Scenario'] == self.scenarioIEA]
        for m in OSD_ScenarioIEA.index:
            for d in self.listDecades:
                OSDy[d].loc[m] = OSD_ScenarioIEA[d].loc[m]

        # Correction with projections for different SSP, matching scenario
        OSD_ScenarioSSP = OSD_litt_byScenario.loc[
            OSD_litt_byScenario['Scenario'] == self.scenario[:4]]  # Match with the corresponding ssp
        for m in OSD_ScenarioSSP.index:
            for d in self.listDecades:
                OSDy[d].loc[m] = OSD_ScenarioSSP[d].loc[m]

        self.OSDy = OSDy

        OSDd = pd.DataFrame(index=self.listMetals, columns=self.listDecades)

        # Initialise the datas from 2020 with actual known metal production at the year 2020, without metal for the transition [t/y]
        OSDd["2020"] = OSDy["2020"]

        for m in self.listMetals:
            # Loop through list of years, excluding the first year 2020
            for y in range(1, len(self.listDecades)):
                # Mean value between the demand at the year y and y-1, devided by 10 to have the value by year
                OSDd[self.listDecades[y]] = (OSDy[self.listDecades[y]] + OSDy[self.listDecades[y - 1]]) * 10 / 2

        self.OSDd = OSDd

    def modelDef(self):
        '''
        Definition of the optimisation model and the variables
        '''

        self.model = ConcreteModel()

        # Declaration of the optimization variable
        self.model.s = Var(self.listRegions, self.listTechno, self.listDecades, initialize=0, domain=NonNegativeReals)
        # Relaxation variables

        # Declaration of the  relaxation variables, one positive, one negative, for the metal reserve constraint
        self.model.Res_relax = Var(self.listMetals_knownRes, initialize=0, within=NonNegativeReals)

        # Declaration of the  relaxation variables, one positive, one negative, for the metal mining constraint
        self.model.Mining_relax = Var(self.listMetals, self.listDecades, initialize=0, within=NonNegativeReals)

        # Declaration of the relaxation variable in case the initial techno mix has a decrease for some technoRen
        self.model.CoherentMix_relax = Var(self.listRegions, self.listTechno_Ren, self.listDecades, initialize=0, within=NonNegativeReals)

    def modelObj(self):
        '''
        Objective of the optimisation model : to minimize the difference with IAM scenario
        '''

        self.model.objective = Objective(
            expr=
            sum(((((sum(self.TechnoMatrix[t].loc[T] * self.model.s[r,t,d] for t in self.listTechno))
                   -self.s0[r].loc[T][d])/self.s0[r].loc[T][d])**2)
                for r in self.listRegions for T in self.listTechnoFlex for d in self.listDecades)
            + sum(self.M * (self.model.Res_relax[m] / self.MetalDatas[self.ResLimit].loc[m]) for m in self.listMetals_knownRes)
            + sum(self.M * (self.model.Mining_relax[m, d] / self.Prod[d].loc[m]) for m in self.listMetals for d in self.listDecades)
            + sum(self.M * (self.model.CoherentMix_relax[r, t_r, d]) for r in self.listRegions for t_r in self.listTechno_Ren for d in self.listDecades)
            , sense=minimize)


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
                        * self.CF_disag[r].loc[T][d] for T in self.listTechnoFlex)
                    >= sum(self.s0[r].loc[T][d] * self.CF_disag[r].loc[T][d] for T in self.listTechnoFlex)
                )

    def CstrCappedFoss(self):
        '''
        Constrain the model to not increase fossil fuel technologies, to limit climate change
        '''
        self.model.constraintCC = ConstraintList()

        for r in self.listRegions:
            for t in self.listTechno_Foss:
                for d in self.listDecades:
                    self.model.constraintCC.add(self.model.s[r, t, d] - self.s0[r].loc[t][d] <= 0)

    def CstrMetalRes(self):
        '''
        Constrain the model to limit cumulated metal demand in 2050 under the ResLimit restriction
        '''

        self.model.constraint_DemandCumulated = ConstraintList()

        for m in self.listMetals_knownRes:
            self.model.constraint_DemandCumulated.add(sum((self.model.s[r,t,self.listDecades[d]]-self.model.s[r,t,self.listDecades[d-1]])
                                          *((self.MetalIntensity_doc[self.listDecades[d]][t].loc[m]+self.MetalIntensity_doc[self.listDecades[d-1]][t].loc[m])/2)
                                          /self.MetalDatas['RR (%) Res'].loc[m] for t in self.listTechno_Ren for r in self.listRegions for d in range(1, len(self.listDecades)))
                                                      + sum(self.OTDd.loc[m][d] / self.MetalDatas['RR (%) Res'].loc[m] for d in self.listDecades)
                                                      + sum(self.OSDd[d].loc[m] for d in self.listDecades)  # Cumulated demand
                                                      - self.model.Res_relax[m]
                                                      <= self.MetalDatas[self.ResLimit].loc[m])

        # model.constraint_DemandCumulated.pprint()

    def CstrMetalMining(self):
        '''
        Constrain the model to limit annual metal demand under the estimated metal production
        or to limit metal demand for energy under an Alpha percentage of the demand for
        the rest of the transition and economy
        '''

        self.model.Constraint_DemandByYear = ConstraintList()

        for d in range(1, len(self.listDecades)):
            for m in self.listMetals:
                self.model.Constraint_DemandByYear.add(sum(self.MetalIntensity_doc[self.listDecades[d]][t].loc[m]
                                                           * (self.model.s[r, t, self.listDecades[d]] -
                                                              self.model.s[r, t, self.listDecades[d - 1]])
                                                           / self.MetalDatas['RR (%) Prod'].loc[m] for t in
                                                           self.listTechno_Ren for r in
                                                           self.listRegions) / 10  # Suppose a linear growth by decade
                                                       <= max(
                    self.Prod[self.listDecades[d]].loc[m]  # Production at the year of the decade
                    - self.OTDy[self.listDecades[d]].loc[m] / self.MetalDatas['RR (%) Prod'].loc[
                        m]  # Metals needed for EV, electric grid, storage, according to IEA
                    - self.OSDy[self.listDecades[d]].loc[m], self.Alpha / 100 * (
                                self.OSDy[self.listDecades[d]].loc[m] + self.OTDy[self.listDecades[d]].loc[m]))
                                                       + self.model.Mining_relax[
                                                           m, self.listDecades[d]])  # Relaxation variable

    def CstrTechnoCoherence(self):
        '''
        Constraint model to keep the initial IAM mix for 2020,
        and to keep the installed capacity during the previous decade
        '''

        # Creation of a list of constraint to add coherence in optimised the technological mix of 2020
        self.model.constraintMixCoherence2020 = ConstraintList()

        for r in self.listRegions:
            for T in self.listTechnoFlex:
                # The initial technological mix of 2020 cannot be changed
                self.model.constraintMixCoherence2020.add(sum(self.TechnoMatrix[t].loc[T] *self.model.s[r,t,'2020'] for t in self.listTechno)
                                                          ==self.s0[r].loc[T]['2020']
                                                          )

        # Creation of a list of constraint to add coherence in the evolution of the technological mix
        self.model.constraintMixCoherence = ConstraintList()

        for r in self.listRegions:
            for t in self.listTechno_Ren:
                for d in range(1, len(self.listDecades)):
                    # If a capacity is installed in d-1, it is still in the cumulated installed capacity in d
                    self.model.constraintMixCoherence.add(
                        self.model.s[r, t, self.listDecades[d]] + self.model.CoherentMix_relax[r, t, self.listDecades[d]] >= self.model.s[r, t, self.listDecades[d - 1]])

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
        self.res = self.model.s.get_values()

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
        self.Obj.to_excel(excel, sheet_name='Obj')
        excel.close()

    def CalculResults (self):
        '''
        Calcul results for initial and optimised power mix, cumulated metal demand and metal demand by year
        '''

        # Stock optimisation results

        # Creation of a dictionary res with the results
        if self.ModelisationType == 'Init':
            self.res_dic = self.s0
        if self.ModelisationType == 'Opti':
            result = self.model.s.get_values()
            res_dic = {}
            for key, value in result.items():
                r, t, d = key
                if r not in res_dic:
                    res_dic[r] = pd.DataFrame()
                # Update the value in the DataFrame at the specified index and column
                res_dic[r].loc[t, d] = value
            self.res_dic = res_dic

        # Calcul cumulated demand
        # Create a list with metals needed for the energy sources
        EnergySectorDemand = pd.Series()
        # Create a list with metals needed for the EV, storage and electric grid
        OtherTransitionDemand = pd.Series()
        # Create a list with metals needed for the other sector of society
        OtherSectorDemand = pd.Series()
        # Add datas for each Series
        for m in self.listMetals:
            EnergySectorDemand[m] = sum(
                (self.res_dic[r].loc[t][self.listDecades[d]] - self.res_dic[r].loc[t][self.listDecades[d - 1]])
                * ((self.MetalIntensity_doc[self.listDecades[d]][t].loc[m] +
                    self.MetalIntensity_doc[self.listDecades[d - 1]][t].loc[m]) / 2)
                / self.MetalDatas['RR (%) Res'].loc[m] for t in self.listTechno_Ren for r in self.listRegions for d
                in range(1, len(self.listDecades)))
            OtherTransitionDemand[m] = sum(
                self.OTDd.loc[m][d] / self.MetalDatas['RR (%) Res'].loc[m] for d in self.listDecades)
            OtherSectorDemand[m] = sum(self.OSDd[d].loc[m] for d in self.listDecades)  # Cumulated demand
            # Concatenate the Series along axis 1 (columns), specifying keys for the resulting DataFrame
        DemandCumulated = pd.concat([OtherSectorDemand, OtherTransitionDemand, EnergySectorDemand], axis=1,
                                    keys=['OtherSectorDemand', 'OtherTransitionDemand', 'EnergySectorDemand'])
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

        # Create a dF with metals needed for the energy sources
        EnergySectorDemandy = pd.DataFrame(0.0, index=self.listMetals, columns=self.listDecades)
        EnergySectorDemandy['2020'] = self.PowerMetalDemand_y['2020']
        # Loop through decades and metals
        for d in range(1, len(self.listDecades)):
            for m in self.listMetals:
                # Add metals need for new installed capacity between the actual decade and the previous one
                EnergySectorDemandy[self.listDecades[d]].loc[m] = sum(
                    self.MetalIntensity_doc[self.listDecades[d]][t].loc[m]
                    * newCapacityD[r][self.listDecades[d]].loc[t] / 10
                    for t in self.listTechno_Ren for r in self.listRegions)
        # Create a three dimensions dataFrame to add all datas of consumption by decade
        DemandByYear_df = {}
        for m in self.listMetals:
            # Index for different datasSets, columns for each decades
            DemandByYear_df[m] = pd.DataFrame(
                index=['OtherSectorDemand', 'OtherTransitionDemand', 'EnergySectorDemand'],
                columns=self.listDecades)
            for i in DemandByYear_df[m].index:
                # Add the values for each dataSets
                DemandByYear_df[m].loc['EnergySectorDemand'] = EnergySectorDemandy.loc[m] / \
                                                                   self.MetalDatas['RR (%) Prod'].loc[m]
                DemandByYear_df[m].loc['OtherSectorDemand'] = self.OSDy.loc[m]
                DemandByYear_df[m].loc['OtherTransitionDemand'] = self.OTDy.loc[m] / \
                                                                  self.MetalDatas['RR (%) Prod'].loc[m]
        self.DemandByYear_df = DemandByYear_df


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

    def Save_OSD_OTD_Prod (self):

        '''
        Save datas for other sector demand and other transition demand by year and cumulated by decade
        Save estimated production by year
        '''

        Res_Folder_BroaderEconomy = self.result_path+ '/BroaderEconomy'
        if not os.path.exists(Res_Folder_BroaderEconomy):
            os.makedirs(Res_Folder_BroaderEconomy)

        # Save OSDd
        # Generate the full path for the Excel file for this scenario and model
        OSDd_Folder = Res_Folder_BroaderEconomy + '/OSDd'
        if not os.path.exists(OSDd_Folder):
            os.makedirs(OSDd_Folder)
        excel_path_OSDd = OSDd_Folder + '/OSDd_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path_OSDd)
        self.OSDd.to_excel(excel, sheet_name='OSDd')
        excel.close()

        # Save OSDy
        # Generate the full path for the Excel file for this scenario and model
        OSDy_Folder = Res_Folder_BroaderEconomy + '/OSDy'
        if not os.path.exists(OSDy_Folder):
            os.makedirs(OSDy_Folder)
        excel_path_OSDy = OSDy_Folder + '/OSDy_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path_OSDy)
        self.OSDy.to_excel(excel, sheet_name='OSDy')
        excel.close()


        # Save OTDy
        # Generate the full path for the Excel file for this scenario and model
        OTDy_Folder = Res_Folder_BroaderEconomy + '/OTDy'
        if not os.path.exists(OTDy_Folder):
            os.makedirs(OTDy_Folder)
        excel_path = OTDy_Folder + '/OTDy_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path)
        self.OTDy.to_excel(excel, sheet_name='OTDy')
        excel.close()

        # Save OTDy
        OTDd_Folder = Res_Folder_BroaderEconomy + '/OTDd'
        if not os.path.exists(OTDd_Folder):
            os.makedirs(OTDd_Folder)
        excel_path = OTDd_Folder + '/OTDd_' + self.model_s0 + '_' + self.scenario + '.xlsx'
        # Create the Excel file and write the data
        excel = pd.ExcelWriter(excel_path)
        self.OTDd.to_excel(excel, sheet_name='OTDd')
        excel.close()

    def SaveEffortBySociety(self):
        '''
        Saves an Excel file, that shows for each decade and each metal, if the mining constraint have been
        respected because total metal demand is under production ('Prod')
        or if the demand of metal for energy is less than Alpha percent of the rest of economy ('Society')
        It would then suppose a requirement of effort from society to decrease consumption
        '''

        ConstraintDemandByYear_file = self.result_path + '/BroaderEconomy/EffortBySociety' + '_' + self.model_s0 + '_' + self.scenario + '.xlsx'

        if not os.path.exists(ConstraintDemandByYear_file):

            ConstraintDemandByYear = pd.DataFrame(index=self.listMetals, columns=self.listDecades[1:])

            for m in self.listMetals:
                for d in self.listDecades[1:]:
                    # If total metal demand is under the estimated production
                    if (self.Prod[d].loc[m] - self.OTDy[d].loc[m] / self.MetalDatas['RR (%) Prod'].loc[m] - self.OSDy[d].loc[m]) >= (
                            self.Alpha / 100 * self.OSDy[d].loc[m]):
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

