## IAM_metal_optimisation 

The python code "IAM_metal_optimisation.py" is a software tool designed, on one hand, to quantify metal supply and demand by sector through 2050, based on energy projections from various IAMs and SSP-RCP scenarios. On the other hand, this code enables the optimization of the minimum variation in the IAM technological mix necessary to limit demand exceedance beyond metal supply constraints.

### How to use the code ? 

1. Clone the repository (or download it)
2. Remove the contents from the 'Datas' folder, so that it is in your "folder path"
3. Download the datas of the IAM you want to use
4. The code is ready to be used !

### How to download datas of the IAM you want

1. You want to study the modelisation from the article
   
The datas needed to generate the results of RCP 2.6, for various SSP, modeled by each IAM markor is available on Zenodo at the following link: 
        — Data: data-zenodo

2. You want to study an other IAM or SSP-RCP scenario

To study an other IAM and SSP-RCP modelisation, use the datas from the IIASA database : https://tntcat.iiasa.ac.at/SspDb/dsd?Action=htmlpage&page=20

After logging in, go to the "IAM Scenarios" section, and download datas for : 
- Regions : choose 5 regions
- Model/Scenarios : choose the scenarios you want to study
- Variable : At least, the GDP, the Secondary Energy (in Energy folder) and the Power Capacity (in Technological                 Indicators Folder > Capacity > Electricity)
The datas are stocked in a csv file for every IAM, every SSP-RCP scenario, every region, for each variable chosen.

To use the data in the code, they have to be downloaded in a folder for each variable. The folder contains a folder for each modelisation (IAM - SSP-RCP), which is filled with 5 Excel files for each region. 
The folder and files name is :
- Capacity Factor IAM > FC models0_scenario > FCmodels0_scenario_region
- GDP IAM > GDP models0_scenario > GDPmodels0_scenario_region
- Power Capacity IAM > Dossier s0 models0_scenario > Cap_models0_scenario_region

### Choose parameters 

The software tool is structured as a class, divided into several methods. To initialize the class, the user must provide seven parameters to the __init__ constructor:
- folder_path: [string] path on the user's computer to the external data used in the code (finishing with /)
- result_path: [string] path on the user's computer for storing the results generated by the code (finishing with /)
- model_s0: [string] IAM that the user wishes to study
- scenario: [string] SSP-RCP scenario combination that the user wants to study, in the format “SSPx-numRCP”
- ModelisationType: [string] ; " Init " to study the metal demand of the initial energy mix ; " Opti " to perform an optimization of the initial energy mix and study the metal demand of the optimized mix
- Alpha: [int] percentage value of the acceptable level of metal demand exceedance in the energy sector compared to OSD and
- OTD for the annual metal production optimization constraint
- ResLimit: [string] limit to apply when optimizing cumulative metal demand, choosing between Reserves or Resources

  All the other methods do not require any parameters. 
