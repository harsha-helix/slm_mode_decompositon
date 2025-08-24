SLM Mode Decomposition Repository
This repository contains the code and notebooks for performing laser mode decomposition using a Spatial Light Modulator (SLM).

#File Descriptions
##mode_decomp.py
This is the core Python script that functions as the main library for this project. It contains all the essential functions required to interface with the hardware, generate phase patterns, and perform the measurements. Key functionalities include:

Connecting to the Holoeye SLM and Basler camera.

Capturing and processing images.

Generating phase-only holograms for Hermite-Gaussian (HG) modes.

Implementing the digital knife-edge technique to find beam properties.

##simulations.ipynb
This Jupyter Notebook is used for running simulations and theoretical analysis related to the mode decomposition process. It likely contains code to model the expected outcomes, visualize theoretical beam profiles, and compare simulated data with experimental results.

##test.ipynb
This Jupyter Notebook serves as the primary interface for running experiments and testing the functions defined in mode_decomp.py. It is used for hands-on control of the SLM and camera, executing the knife-edge scans, displaying holograms, and collecting data in a step-by-step, interactive manner.

## simulations.ipynb
This Jupyter Notebook is used to simulate interaction of the Phase mask on SLM with the input beam to obtain intensity distribution after free space propogation.
