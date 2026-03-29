"""
H/C, O/C, and HCV Predictor
Copyright (C) 2026 Nitish Kapur
GitHub: github.com/nitish-kapur
Licensed under GNU GPLv3

    This script was made as a part of a biofuel research project.

    1.  Defines the elemental compositions (C, H, O, N, S as % mass) for
        each feedstock in the `biomass` dictionary. Update these values to
        match your own feedstocks and experimental data.

    2.  Prompts the user to enter the mass of each feedstock in grams,
        computes the total mass, and converts each input to a mass
        percentage.

    3.  Implements three property calculators for a given blend vector x
        (mass percentages of each feedstock):
            - H/C atomic ratio: molar hydrogen to carbon ratio
            - O/C ratio: mass oxygen to carbon ratio
            - HCV using Dulong's formula (MJ/kg)
            - HCV using Modified Dulong's formula (MJ/kg)

    4.  Prints the following results to the console:
            - Mass percentages and total weight of the blend
            - HCV from both Dulong's and Modified Dulong's formulae
            - H/C atomic ratio
            - O/C ratio
"""

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

"""
        Elemental (CHNS) analysis (% mass) of rice straw (PR 121), LDPE (Verka milk packets), 
        coconut shells (South Indian), and walnut shells (Kashmir, India), was performed using 
        an Elementar Vario EL Cube CHNS Analyser.
"""

biomass = {
    'rice_straw':  {'C': 37.51, 'H': 10.043,  'N': 0.52,  'O': 51.560, 'S': 0.368},
    'LDPE':       {'C': 79.22, 'H': 8.253, 'N': 0.0, 'O': 12.301,  'S': 0.226},
    'coconut':    {'C': 46.53,  'H': 5.829,  'N': 0.30, 'O': 46.991,  'S': 0.350},
    'walnut':     {'C': 45.53,  'H': 5.829,   'N': 0.30,  'O': 47.893, 'S': 0.254},
}

names = ['rice_straw', 'LDPE', 'coconut', 'walnut']
C_atomic_mass = 12
H_atomic_mass = 1
O_atomic_mass = 16


# Function to calculate H/C ratio
def calc_HC_ratio(x):
    C_mass = sum(biomass[n]['C'] * x[i] / 100 for i, n in enumerate(names))
    H_mass = sum(biomass[n]['H'] * x[i] / 100 for i, n in enumerate(names))
    C_moles = C_mass / C_atomic_mass
    H_moles = H_mass / H_atomic_mass
    if C_moles == 0:
        return 0
    return H_moles / C_moles


# Function to calculate O/C ratio
def calc_OC_ratio(x):
    C_mass = sum(biomass[n]['C'] * x[i] / 100 for i, n in enumerate(names))
    O_mass = sum(biomass[n]['O'] * x[i] / 100 for i, n in enumerate(names))
    if C_mass == 0:
        return 0
    return O_mass / C_mass


# Function to calculate HCV using Dulong's Formula
def calc_HCV_Dulong(x):
    total_mass = sum(x)
    C_frac = sum(biomass[n]['C'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    H_frac = sum(biomass[n]['H'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    O_frac = sum(biomass[n]['O'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    S_frac = sum(biomass[n]['S'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    HCV = (33.87 * C_frac) + 122.3 * (H_frac - O_frac / 8) + 9.4 * S_frac
    return HCV


# Function to calculate HCV using Modified Dulong's Formula
def calc_HCV_ModDulong(x):
    total_mass = sum(x)
    C_frac = sum(biomass[n]['C'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    H_frac = sum(biomass[n]['H'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    O_frac = sum(biomass[n]['O'] * x[i] / 100 for i, n in enumerate(names)) / total_mass
    LCV = (38.2 * C_frac) + 84.9 * (H_frac - O_frac / 8) - 0.5 # kJ/g = MJ/kg
    HCV_Mod = (LCV + 0.024 * 2 * H_frac) # 0.024 = Latent heat of vaporization of water (kJ/g)
    # 0.024×2 is used because the molar heat of vaporization of water is approximately 44 kJ/mol, and for each gram of hydrogen, about 2 grams of water are produced
    return HCV_Mod


# Function to get user input in grams and convert to percentages
def get_input_in_grams():
    print("Please enter the mass of each precursor in grams:\n")

    # Get input in grams
    rice_straw_grams = float(input("Rice Straw (grams): "))
    LDPE_grams = float(input("LDPE (grams): "))
    coconut_grams = float(input("Coconut (grams): "))
    walnut_grams = float(input("Walnut (grams): "))

    # Calculate total mass
    total_grams = rice_straw_grams + LDPE_grams + coconut_grams + walnut_grams

    # Convert to percentages
    rice_straw_percent = (rice_straw_grams / total_grams) * 100
    LDPE_percent = (LDPE_grams / total_grams) * 100
    coconut_percent = (coconut_grams / total_grams) * 100
    walnut_percent = (walnut_grams / total_grams) * 100

    # Print total weight
    print(f"\nTotal Weight of Feedstocks: {total_grams:.2f} grams")

    # Print out the calculated percentages
    print("\nMass Percentages of Feedstocks:")
    print(f"Rice Straw: {rice_straw_percent:.2f}%")
    print(f"LDPE: {LDPE_percent:.2f}%")
    print(f"Coconut: {coconut_percent:.2f}%")
    print(f"Walnut: {walnut_percent:.2f}%")

    return [rice_straw_percent, LDPE_percent, coconut_percent, walnut_percent]


# Main function to calculate and display the results
def main():
    mass_percentages = get_input_in_grams()  # Get mass percentages by converting grams to %

    # Calculate HCV, H/C ratio, and O/C ratio
    hcv = calc_HCV_Dulong(mass_percentages)
    hcv_mod = calc_HCV_ModDulong(mass_percentages)
    hc_ratio = calc_HC_ratio(mass_percentages)
    oc_ratio = calc_OC_ratio(mass_percentages)

    # Display the results
    print("\nResults:")
    print(f"HCV as per Dulong's formula: {hcv:.2f} MJ/kg")
    print(f"HCV as per Modified Dulong's Formula: {hcv_mod:.2f} MJ/kg")
    print(f"H/C Ratio: {hc_ratio:.2f}")
    print(f"O/C Ratio: {oc_ratio:.2f}")


# Run the main function
if __name__ == "__main__":
    main()
