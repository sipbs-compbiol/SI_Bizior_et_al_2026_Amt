# load_dataframes.R
#
# Provides the load_dataframes() function that returns the four dataframes
# required for the analysis

load_dataframes <- function() {
  # Read in WT data
  dfm_wt_1 <- readr::read_csv("assets/data/01_WT_solvent_exchange_conditions.csv", col_names=TRUE, col_types="fffcdd") |>
    dplyr::mutate(LPR = stringr::str_c("LPR", LPR)) |>
    dplyr::mutate(Variant = "WT")
  dfm_wt_1$LPR <- factor(dfm_wt_1$LPR, levels=c("LPR5", "LPR10", "LPR50"))
  colnames(dfm_wt_1) <- c("Condition", "LPR", "Replicate", "Trace_ID", "Max_Amplitude", "Decay_Rate", "Variant")
  dfm_wt_1 <- dfm_wt_1 |> dplyr::select(Variant, Condition, LPR, Max_Amplitude, Decay_Rate)
  
  dfm_wt_2 <- readr::read_csv("assets/data/02_WT_variant_comparison.csv", col_names=TRUE, col_types="fccddcc") |>
    dplyr::mutate(LPR = stringr::str_c("LPR", LPR)) |>
    dplyr::mutate(Variant = "WT")
  dfm_wt_2$LPR <- factor(dfm_wt_2$LPR, levels=c("LPR5", "LPR10", "LPR50"))
  colnames(dfm_wt_2) <- c("LPR", "Trace_ID", "Measurement", "Max_Amplitude", "Decay_Rate", "Amp_Source", "Decay_Source", "Variant")
  dfm_wt_2 <- dfm_wt_2 |> dplyr::select(Variant, LPR, Max_Amplitude, Decay_Rate)
  
  # Read in D162A data
  dfm_d162 <- readr::read_csv("assets/data/03_D162A_variant_comparison.csv",
                              col_names=TRUE, col_types="ccfcfcfcdc")
  colnames(dfm_d162) <- c("Sheet", "SourceID", "Variant", "Substrate", "LPR", "Replicate", "Parameter", "Unit", "Value", "Status")
  dfm_d162 <- dfm_d162 |> tidyr::pivot_wider(id_cols=c(SourceID, Replicate, Variant, LPR), names_from=Parameter, values_from=Value)
  colnames(dfm_d162) <- c("SourceID", "Replicate", "Variant", "LPR", "Max_Amplitude", "Decay_Rate")
  dfm_d162 <- dfm_d162 |> dplyr::select(Variant, LPR, Max_Amplitude, Decay_Rate)
  
  # Read in H170 variant data
  dfm_h170 <- readr::read_csv("assets/data/04_H170_variant_comparison.csv",
                              col_names=TRUE, col_types="cfcfffcfd") |>
    dplyr::mutate(LPR = stringr::str_c("LPR", LPR))
  dfm_h170$LPR <- factor(dfm_h170$LPR, levels=c("LPR5", "LPR10"))
  colnames(dfm_h170) <- c("Series", "Variant", "Substrate", "Condition", "LPR", "Parameter", "Unit", "ReplicateID", "Value")
  dfm_h170 <- dfm_h170 |>
    dplyr::select(Variant, ReplicateID, Condition, LPR, Parameter, Value) |>
    tidyr::pivot_wider(id_cols=c(Variant, ReplicateID, Condition, LPR),
                       names_from=Parameter, values_from=Value) |>
    dplyr::filter(Variant != "WT")
  colnames(dfm_h170) <- c("Variant", "ReplicateID", "Condition", "LPR", "Max_Amplitude", "Decay_Rate")
  
return(list(dfm_wt_1=dfm_wt_1,
            dfm_wt_2=dfm_wt_2,
            dfm_d162=dfm_d162,
            dfm_h170=dfm_h170))
}