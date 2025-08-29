import streamlit as st
from PIL import Image
import os
import pandas as pd

# Definition of the scope of this study in terms of models, scenarios, regions, decades, metals and techno
studyScope = 'C:/Users/Penel/Documents/Travail/CIRAIG/Maitrise Recherche/8. Optimisation/Scope of the study.xlsx'

with pd.ExcelFile(studyScope) as file:
    # list of the IAM models you can optimise
    listModels = pd.read_excel(file, 'model', index_col=0).squeeze().tolist()
    # list of the SSP - RCP scenarios you can optimise
    listScenarios = pd.read_excel(file, 'scenario', index_col=0).squeeze().tolist()
result_path = 'C:/Users/Penel/Documents/Travail/CIRAIG/Maitrise Recherche/9. Résultats/4- Review/Final_All_Scenario/'

# App title
st.title("Visualize the results for a specific model and scenario")

# Select model and scenario
model = st.selectbox("Choose a model:", listModels)
scenario = st.selectbox("Choose an ssp scenario:", listScenarios)

# File and figure names
image_configs = {
    "Comparing cumulated demand to resources and reserves": {
        "folder": f"{result_path}Resource_images",
        "filename": f"Fig_Resource_{model} - {scenario}.png"
    },
    "Annual demands and mining capacities": {
        "folder": f"{result_path}Mining_images",
        "filename": f"Fig_Mining_{model} - {scenario}.png"
    },

    "Power plant market shares": {
        "folder": f"{result_path}Power_images",
        "filename": f"Fig_PowerComparison_{model} - {scenario}.png"
    },
    "Electric vehicle motor market shares": {
        "folder": f"{result_path}Motor_images",
        "filename": f"Fig_MotorComparison_{model} - {scenario}.png"
    },
    "Electric vehicle battery market shares": {
        "folder": f"{result_path}Battery_images",
        "filename": f"Fig_BatteryComparison_{model} - {scenario}.png"
    }
}

# Loop
for title, cfg in image_configs.items():
    st.subheader(title)
    filepath = os.path.join(cfg["folder"], cfg["filename"])

    if os.path.exists(filepath):
        image = Image.open(filepath)
        st.image(image, caption=f"{model} - {scenario}", use_container_width=True)
    else:
        st.warning(f"No existing figure in {title} for {model} - {scenario}")

