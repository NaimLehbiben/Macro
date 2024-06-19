# Macroeconomic and Monetary Policy Project: Economic Impacts in the Eurozone Due to Rising Temperatures

## Objectives

This project is part of the Master's program in Economics & Financial Engineering at Paris Dauphine University - PSL. The aim is to analyze the economic impacts of rising temperatures on the Eurozone, focusing on inflation and GDP growth, as well as the resulting monetary policy challenges. 

## Project Structure

Here's an overview of the project's structure:

### Description

- **data/**: Contains climate and macroeconomic datasets.
- **src/**: Contains all the source code for the project.
  - `Macro.m`: The main script for running the analysis.
  - `neweyWest_covmatrix.m`: Script for computing the Newey-West covariance matrix.
  - `plot_IRFs.m`: Script for plotting the Impulse Response Functions (IRFs).
- **static/**: Contains static files such as the research paper, report, and images.
- **README.md**: The readme file.

## Methodology

### Data Collection

The data for this project includes climate and macroeconomic data for 12 Eurozone countries from 1996 to 2021. The climate data, including temperature and precipitation, was obtained from the Climate Change Knowledge Portal. The macroeconomic data includes real GDP per capita and the Harmonized Index of Consumer Prices (HICP).

### Analysis

1. **Temperature Anomaly Calculation**: The temperature anomaly for each country and month is calculated as the deviation from the historical average temperature (1950-1980).

2. **Modeling**: The impact of temperature anomalies on inflation and GDP is modeled using regression equations, incorporating Newey-West standard errors to account for autocorrelation and heteroscedasticity.

3. **Impulse Response Functions (IRFs)**: The IRFs are computed to measure the response of inflation and GDP to temperature shocks over different horizons.

4. **Projections**: Using climate projections (SSP scenarios), the future impacts on inflation and GDP are estimated for the period 2025-2100.

5. **Monetary Policy Implications**: The Taylor rule is applied to determine the appropriate interest rate policy in response to temperature-induced economic changes, and the concept of monetary stress is evaluated.

## Results

The results indicate that temperature anomalies have varying impacts on inflation and GDP across different Eurozone countries. The findings suggest that a one-size-fits-all monetary policy by the ECB may be insufficient to address the diverse effects of climate change within the Eurozone.

## Interface

The project includes scripts for generating visualizations and reports on the analysis. Users can run the main script `Macro.m` to perform the analysis and generate the results.

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Authors

- Antoine Portalier
- Badr Eddine El Hamzaoui
- Naïm Lehbiben
- Vincent Piedbois
