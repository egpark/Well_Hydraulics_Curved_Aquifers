# Well_Hydraulics_Curved_Aquifers
This repository contains MATLAB code for implementing the geodesic distance approach in well hydraulics analysis. The code demonstrates the application of geodesic distances to enhance the accuracy and applicability of conventional analytical solutions for groundwater flow in confined aquifers with complex geometries.

Files
-----
- main_wh_geodesic_aniso.m: The main script file that orchestrates the entire analysis. It sets the domain parameters, well parameters, and other necessary settings. This file calls the various functions and scripts to perform the analysis and generate the results.

- load_data.m: This script file is responsible for loading the pointwise data of aquifer heights and orientations. The data is stored in the variables 'dat_pnt' and 'grad_pnt', respectively. Practitioners can adjust the aquifer height and orientation (strike) information in this file to suit their specific site requirements.

- comp_d_g_aniso.m: This function file contains the implementation of the geodesic distance computation. It calculates the geodesic distances between the pumping well and observation points in the aquifer domain, taking into account the intrinsic geometry of the aquifer.

- GPR_est.m: This function file performs the Gaussian Process Regression (GPR) estimation for interpolating the aquifer topography and gradients. It utilizes the pointwise data loaded from 'load_data.m' and applies the geodesic kernel to capture the intrinsic geometry of the aquifer.

- draw_figures.m: This script file is responsible for visualizing the results of the analysis. It generates various figures, including the aquifer topography, drawdown patterns, and comparisons with numerical simulations (if available).

Usage
-----
1. Ensure that all the necessary files ('main_wh_geodesic_aniso.m', 'load_data.m', 'comp_d_g_aniso.m', 'GPR_est.m', and 'draw_figures.m') are in the same directory.

2. Open 'load_data.m' in MATLAB and adjust the aquifer height and orientation (strike) information according to your specific site requirements. The data is provided in a tabular format, with each row representing a data point and the columns corresponding to the x-coordinate, y-coordinate, height (for 'dat_pnt'), and strike vector components (for 'grad_pnt').

3. Open 'main_wh_geodesic_aniso.m' in the MATLAB editor.

4. Adjust the domain parameters, well parameters, and other settings as needed.

5. Run 'main_wh_geodesic_aniso.m' to execute the analysis. The script will load the data, perform the necessary computations, and generate the results.

6. View the generated figures to analyze the aquifer topography, drawdown patterns, and other relevant information.

Dependencies
------------
The code requires MATLAB to be installed on your system. It has been tested with MATLAB version 1. No additional external dependencies are required.

Note
----
The code provided in this repository is for illustrative purposes and may need to be adapted to suit your specific use case. Make sure to review and modify the code as necessary to align with your research or project requirements.

For more details on the methodology and its applications, please refer to the accompanying research paper: Utilizing Geodesic Distances for Enhanced Aquifer Hydraulic Analysis (Wen et al., WRR, 2024).
