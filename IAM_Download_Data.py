

# Import usual libraries
import pandas as pd
import numpy as np
import os
import gzip
import pickle
import matplotlib.pyplot as plt
import math
import xlsxwriter

# Librairies de représentation graphique
import plotly.express as px
import plotly.graph_objs as go
import matplotlib.pyplot as plt

class IAM_Download_Data :

    def __init__(self, file_path_iam_results, folder_path, model_s0, scenario):

        '''
        :param file_path_iam_results: [string] path to the csv file with IAM datas from IIASA database
        :param folder_path: [string] path to the folder where you want to store IAM datas,
        with external datas used in the code
        :param model_s0: Integrated Assessment Model (IAM) that wants to be studied
        (AIM..CGE, GCAM4, IMAGE, MESSAGE-GLOBIOM, REMIND-MAGPIE, WITCH-GLOBIOM)
        :param scenario: scenario of Shared Socioeconomic Pathway and  Representative Concentration Pathway
        (SSP-RCP) that wants to be studied
        '''

        self.file_path_iam_results = file_path_iam_results
        self.folder_path = folder_path
        self.model_s0 = model_s0
        self.scenario = scenario

        self.listRegions = ['ASIA', 'LAM', 'MAF', 'OECD', 'REF']

        # Create dataFrame with csv file of SSP results
        with open(self.file_path_iam_results) as csvfile:
            self.results_read = pd.read_csv(csvfile, sep=';')

    def Download_Cap (self) :

        capacities = {}

        for r in self.listRegions:
            capacity_df = self.results_read[
                self.results_read['MODEL'].str.contains(self.model_s0) & self.results_read['SCENARIO'].str.contains(self.scenario) &
                self.results_read['REGION'].str.contains(r) & self.results_read['VARIABLE'].str.contains("Capacity..Electricity")]
            capacity_df = capacity_df.drop("MODEL", axis=1).drop("SCENARIO", axis=1).drop("REGION", axis=1).drop("UNIT",
                                                                                                                 axis=1)
            capacity_df = capacity_df.set_index("VARIABLE")
            capacities[self.model_s0, self.scenario, r] = capacity_df

            # Create a file to put the excels of the regions in it
            folder_CAP = self.folder_path + 'Power Capacity IAM/Dossier s0 ' + self.model_s0 + '_' + self.scenario

            # Verify if the file does not exist
            if not os.path.exists(folder_CAP):
                os.makedirs(folder_CAP)

                # Save capacities in an excel file of technologies for a particular model, scenario, region
            table = capacities[self.model_s0, self.scenario, r]

            excel = pd.ExcelWriter(
                self.folder_path + 'Power Capacity IAM/Dossier s0 ' + self.model_s0 + '_' + self.scenario + '/ ' + 'cap' + self.model_s0 + '_' + self.scenario + '_' + r + '.xlsx')
            table.to_excel(excel, sheet_name='results')

            excel.close()

        self.capacities =capacities

    def Download_CF(self):

        productions = {}

        for r in self.listRegions:
            production_df = self.results_read[
                self.results_read['MODEL'].str.contains(self.model_s0) & self.results_read['SCENARIO'].str.contains(self.scenario) &
                self.results_read['REGION'].str.contains(r) & self.results_read['VARIABLE'].str.contains(
                    "Secondary Energy..Electricity")]
            production_df = production_df.drop("MODEL", axis=1).drop("SCENARIO", axis=1).drop("REGION", axis=1).drop(
                "UNIT", axis=1)
            production_df = production_df.set_index("VARIABLE")
            productions[self.model_s0, self.scenario, r] = production_df

        capacity_utilization = {}

        # Facteur de charge : energie / (capacité * 0.031536) * 100
        # 0.031536 : 3600*24*365.25 et conversion GW / EJ

        for r in self.listRegions:
            capacity_utilization[self.model_s0, self.scenario, r] = pd.DataFrame(
                columns=self.capacities[self.model_s0, self.scenario, r].columns)
            if (('Secondary Energy..Electricity..' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Total', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                              'Secondary Energy..Electricity..', :] / (
                                                                                          self.capacities[
                                                                                              self.model_s0, self.scenario, r].loc[
                                                                                          'Capacity..Electricity..',
                                                                                          :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Biomass..' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Biomass' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Biomass', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                                'Secondary Energy..Electricity..Biomass..',
                                                                                :] / (self.capacities[
                                                                                          self.model_s0, self.scenario, r].loc[
                                                                                      'Capacity..Electricity..Biomass',
                                                                                      :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Coal..' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Coal' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Coal', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                             'Secondary Energy..Electricity..Coal..',
                                                                             :] / (self.capacities[
                                                                                       self.model_s0, self.scenario, r].loc[
                                                                                   'Capacity..Electricity..Coal',
                                                                                   :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Gas..' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Gas' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Gas', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                            'Secondary Energy..Electricity..Gas..',
                                                                            :] / (self.capacities[self.model_s0, self.scenario, r].loc[
                                                                                  'Capacity..Electricity..Gas',
                                                                                  :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Geothermal' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Geothermal' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Geothermal', :] = productions[
                                                                                       self.model_s0, self.scenario, r].loc[
                                                                                   'Secondary Energy..Electricity..Geothermal',
                                                                                   :] / (self.capacities[
                                                                                             self.model_s0, self.scenario, r].loc[
                                                                                         'Capacity..Electricity..Geothermal',
                                                                                         :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Hydro' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Hydro' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Hydro', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                              'Secondary Energy..Electricity..Hydro',
                                                                              :] / (self.capacities[
                                                                                        self.model_s0, self.scenario, r].loc[
                                                                                    'Capacity..Electricity..Hydro',
                                                                                    :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Nuclear' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Nuclear' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Nuclear', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                                'Secondary Energy..Electricity..Nuclear',
                                                                                :] / (self.capacities[
                                                                                          self.model_s0, self.scenario, r].loc[
                                                                                      'Capacity..Electricity..Nuclear',
                                                                                      :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Solar' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Solar' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Solar', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                              'Secondary Energy..Electricity..Solar',
                                                                              :] / (self.capacities[
                                                                                        self.model_s0, self.scenario, r].loc[
                                                                                    'Capacity..Electricity..Solar',
                                                                                    :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Wind' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Wind' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Wind', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                             'Secondary Energy..Electricity..Wind',
                                                                             :] / (self.capacities[
                                                                                       self.model_s0, self.scenario, r].loc[
                                                                                   'Capacity..Electricity..Wind',
                                                                                   :] * 0.031536) * 100
            if (('Secondary Energy..Electricity..Oil' in productions[self.model_s0, self.scenario, r].index) & (
                    'Capacity..Electricity..Oil' in self.capacities[self.model_s0, self.scenario, r].index)):
                capacity_utilization[self.model_s0, self.scenario, r].loc['Oil', :] = productions[self.model_s0, self.scenario, r].loc[
                                                                            'Secondary Energy..Electricity..Oil', :] / (
                                                                                        self.capacities[
                                                                                            self.model_s0, self.scenario, r].loc[
                                                                                        'Capacity..Electricity..Oil',
                                                                                        :] * 0.031536) * 100

            # Save CF of technologies for a particular model, scenario, region
            table = capacity_utilization[self.model_s0, self.scenario, r]

            # Create a file to put the excels of the regions in it
            folderCF = self.folder_path + 'Capacity Factor IAM/FC' + self.model_s0 + '_' + self.scenario
            if not os.path.exists(folderCF):
                os.makedirs(folderCF)

            # Generate the full path for the Excel file for this region
            excel_path = self.folder_path + 'Capacity Factor IAM/FC' + self.model_s0 + '_' + self.scenario + '/' + 'FC' + self.model_s0 + '_' + self.scenario + '_' + r + '.xlsx'

            # Create the Excel file and write the data
            excel = pd.ExcelWriter(excel_path)
            table.to_excel(excel, sheet_name='results')
            excel.close()

    def Download_GDP(self):
        # Initialize a dictionary to stock GDP datas from IAM database
        GDP = {}

        for r in self.listRegions:
            # Take data from the column GDP, with GDP name
            GDP_df = self.results_read[
                self.results_read['MODEL'].str.contains(self.model_s0) & self.results_read['SCENARIO'].str.contains(self.scenario) & self.results_read[
                    'REGION'].str.contains(r) & self.results_read['VARIABLE'].str.contains("GDP")]
            GDP_df = GDP_df.drop("MODEL", axis=1).drop("SCENARIO", axis=1).drop("REGION", axis=1).drop("UNIT", axis=1)
            GDP_df = GDP_df.set_index("VARIABLE")
            GDP[self.model_s0, self.scenario, r] = GDP_df

            # Save GDP for a particular model, scenario, region
            table = GDP[self.model_s0, self.scenario, r]

            # Create a file to put the excels of the regions in it
            folderGDP = self.folder_path + 'GDP IAM/GDP ' + self.model_s0 + '_' + self.scenario
            if not os.path.exists(folderGDP):
                os.makedirs(folderGDP)

            # Generate the full path for the Excel file for this region
            excel_path = self.folder_path + 'GDP IAM/GDP ' + self.model_s0 + '_' + self.scenario + '/' + 'GDP' + self.model_s0 + '_' + self.scenario + '_' + r + '.xlsx'

            # Create the Excel file and write the data
            excel = pd.ExcelWriter(excel_path)
            table.to_excel(excel, sheet_name='results')
            excel.close()



