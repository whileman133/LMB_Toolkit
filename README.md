# Modeling Toolkit for Lithium-Metal Battery Cells

This is the MATLAB toolkit developed for the M.S.E.E. thesis entitled "Modeling and Parameter Estimation Strategies for Rechargeable Lithium-Metal Battery Cells" by W. Hileman. This code base is organized by thesis chapter (`CHxx/` folders) in addition to the following:

- `GEN2_TFS/`: Functions that compute the frequency-response response of internal cell variables by evaluating transfer functions (TFs) derived from the Newman model.
- `GEN2_UTILITY/`: The core utility functions that make up the toolkit.
- `GEN2_XLSX_CELLDEFS/`: Excel spreadsheets containing parameter values for simulation cells.
- `GEN2_XLSX_XRACTL/`: Excel spreadsheets configuring the Hybrid Realization Algorithm (HRA) for model-order reduction.
- `DEMO/`: MATLAB scripts demonstrating usage of the toolkit.

Note: The MATLAB LiveLink is required to run the COMSOL simulations.

More information can be found in the open-access thesis at <https://www.proquest.com/dissertations-theses/physics-based-modeling-parameter-estimation/docview/3050214167/se-2>.