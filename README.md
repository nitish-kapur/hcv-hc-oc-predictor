# H/C, O/C, and HCV Predictor

A Python script that takes feedstock masses as input and calculates the H/C ratio, O/C ratio, and Higher Calorific Value (HCV) of a multi-feedstock biomass blend using both Dulong's and Modified Dulong's formulae.

## Author

**Nitish Kapur**<br>
GitHub: [github.com/nitish-kapur](https://github.com/nitish-kapur)

## Feedstocks

The elemental compositions (C, H, O, N, S as % mass) for each feedstock are hard-coded in the `biomass` dictionary at the top of the script. Update these values to match your own feedstocks and experimental data.

## Requirements
```bash
pip install numpy scipy matplotlib
```

## Usage

Run the script directly:
```bash
python hcv-hc-oc-predictor.py
```

The script will prompt the user to enter the mass of each feedstock in grams. It then calculates and prints the blend percentages, H/C ratio, O/C ratio, and HCV to the console.

## Output

- Mass percentages of each feedstock
- Total feedstock weight (grams)
- HCV calculated using Dulong's formula (MJ/kg)
- HCV calculated using Modified Dulong's formula (MJ/kg)
- H/C atomic ratio
- O/C ratio

## Notes

- The number of feedstocks and their names are also hard-coded in the `get_input_in_grams()` function — update the input prompts if you add or remove feedstocks.
- HCV is estimated using both Dulong's formula and a Modified Dulong's formula.
