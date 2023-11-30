# Modeling Toolkit for Lithium-Metal Battery Cells

This is the MATLAB toolkit developed by the UCCS team for the Modeling Focus of the Capitalizing on Lithium-Metal Battery Cells (CLiMB) project.  Find a list of the most useful assets below.

- `OCP/` Assets for estimating OCP with MSMR model.
  - `buildOCP.m` Pre-processes OCP data collected by a Gamry Potentiostat
  - `fitOCP.m` Fit MSMR model to laboratory-derived OCP estimates (from `buildOCP.m`)
- `NLEIS/` Assets for estimating linear parameters with linear EIS and reaction symmetry with nonlinear EIS (NLEIS)
  - `fitEIS.m` Regress linear EIS model to spectra collected in the laboratory.
  - `runGPRLinearEIS.m` Produce estimates of solid diffusivity and charge-transfer resistance using Gaussian Process Regression (GPR). (Accepts output of `fitEIS.m` as input.)
  - `fitNLEIS.m` Regress nonlinear EIS model to spectra collected in the laboratory.
  - `simFOMNLEIS_socSeries.m` Simulate medium-signal sinusoidal response of the full-order LMB model in COMSOL and save results to disk. Useful for validating the transfer-function and nonlinear impedance models. **Requires COMSOL LiveLink for MATLAB to be running**.
  - `plotSimFOMNLEIS` Plot results of medium-signal EIS simulation for full-order LMB cell and compare to TF and nonlinear models.
  - `plotLinEISSensitivity.m` Plot variation in impedance predicted by transfer-function model with perturbation in the values of certain parameters.
- `GITT/` Assets for estimating diffusivity from GITT experiments.
  - `runProcessGITT.m` Process raw GITT data collected from Gamry Potentiostat, estimate diffusivity, and plot the results.
  - `compareDs_GITT_EIS.m` Compare diffusivity estimates from GITT to those obtained using linear EIS.
- `ROM/` Assets for generating reduced-order models with the Hybrid Realization Algorithm (HRA)
  - `makeROM.m` Make reduced-order model for LMB battery cell from parameter estimates derived from linear EIS.
  - `compareROM_Lab_GITT.m` Compare prediction of ROM generated using `makeROM.m` to laboratory GITT cycle.
  - `compareROM_Lab_HalfCycleDischarge.m` Compare prediction of ROM generated using `makeROM.m` to laboratory half-sine current cycle.
  - `plotCompareROM_FOM_LAB.m` Plot results of `compareROM_Lab.m`
- `RPT/` Assets for parameterizing aged cells (incomplete).
- `GEN2_TFS/` Functions that compute transfer-function response of internal cell variables.
- `GEN2_UTILITY/` Utility functions.
- `GEN2_XLSX_CELLDEFS/` Excel spreadsheets containing parameter values for simulation cells.
- `GEN2_XLSX_XRACTL/` Excel spreadsheets configuring the HRA.