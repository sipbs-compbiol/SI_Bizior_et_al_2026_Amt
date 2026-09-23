# README.md SI_Bizior_et_al_2026_Amt

This repository describes supplementary statistical anakysis for the Bizior _et al._ (2026) manuscript "A bacterial Rhesus transporter retunes a structurally conserved ammonium pore into a reversible nitrogen valve."

The analysis proceeds from the data file `NeRh50_Table_S1_raw_data_used_in_figures.xlsx` provided by Arnaud Javelle (2026-09-02). This comprises a single table, exported into the file `NeRh50_Table_S1_raw_data_used_in_figures.csv` with the columns:

- Series: integer; 1 = variant experiment, 2 = solvent exchange experiment
- Figure number/Panel: string; the figures that use the data
- Variant: factor; the NeRh50 protein variant
- Condition: factor; solvent
- LPR: factor; lipid-protein ratio
- Parameter: factor; the measurement taken
- Unit: string; units for the measurement
- Value: double; measured value



