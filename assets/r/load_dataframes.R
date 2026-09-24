# load_dataframes.R
#
# Provides the load_dataframes() function that returns the four dataframes
# required for the analysis

load_dataframes <- function() {
  # Read in combined data
  dfm_combined <- readr::read_csv("assets/data/NeRh50_Table_S1_raw_data_used_in_figures.csv", col_names=TRUE, col_types="fcffffcd") |>
    dplyr::mutate(LPR = stringr::str_c("LPR", LPR)) |>  # add prefix to LPR factor
    dplyr::mutate(Series = stringr::str_c("Series ", Series))  # add prefix to Series factor
  dfm_combined$LPR <- factor(dfm_combined$LPR, levels=c("LPR5", "LPR10", "LPR50"))  # order LPR factor levels
  dfm_combined$Condition <- factor(dfm_combined$Condition, levels=c("H2O", "D2O", "D2O>H2O"))  # order LPR factor levels
  colnames(dfm_combined) <- c("Series", "Figures", "Variant", "Condition", "LPR", "Measurement", "Unit", "Value")  # rename headers
  
  # Supporting figure that is not in the main manuscript
  # For the wild-type comparison between series we need to combine measurements
  # that meet the criteria:
  # WT variant
  # Series 1 or Series 2
  # Condition H2O
  dfm_wt_amp <- dfm_combined |>  # amplitude measurements
    dplyr::filter(Variant == "WT" & LPR %in% c("LPR5", "LPR10") & Condition == "H2O" & Measurement == "Maximum amplitude") |>
    dplyr::select(Series, Variant, Condition, LPR, Value)
  dfm_wt_k<- dfm_combined |>  # rate measurements
    dplyr::filter(Variant == "WT" & LPR %in% c("LPR5", "LPR10") & Condition == "H2O" & Measurement == "k") |>
    dplyr::select(Series, Variant, Condition, LPR, Value)
  
  # Figure 2: NeRh50 generates D2O-resistant, LPR-dependent electrogenic currents in response to NH4+ stimulation
  # For figure 2 we combine measurements that meet criteria:
  # WT Variant
  # Series 2
  dfm_fig2_amp <- dfm_combined |>  # amplitude measurements
    dplyr::filter(Variant == "WT" & Series == "Series 2" & Measurement == "Maximum amplitude") |>
    dplyr::select(Variant, Condition, LPR, Value)
  dfm_fig2_k<- dfm_combined |>  # rate measurements
    dplyr::filter(Variant == "WT" & Series == "Series 2" & Measurement == "k") |>
    dplyr::select(Variant, Condition, LPR, Value)
  
  # Figure 3: D162A uncouples ammonium binding from translocation.
  # For figure 3 we combine measurements that meet criteria:
  # Series 2
  # WT or D162 Variant
  # LPR5 or LPR10
  dfm_fig3_amp <- dfm_combined |>  # amplitude measurements
    dplyr::filter(Variant %in% c("WT", "D162A") & Series == "Series 1" & LPR %in% c("LPR5", "LPR10") & Measurement == "Maximum amplitude") |>
    dplyr::select(Variant, Condition, LPR, Value)
  dfm_fig3_k <- dfm_combined |>  # rate measurements
    dplyr::filter(Variant %in% c("WT", "D162A") & Series == "Series 1" & LPR %in% c("LPR5", "LPR10") &  Measurement == "k") |>
    dplyr::select(Variant, Condition, LPR, Value)
  
  # Figure 5: Electrophysiological characterisation of H170 variants.
  # For figure 5 we combine measurements that meet criteria:
  # WT or H170 Variant
  # Series 1
  # LPR5 or LPR10
  dfm_fig5_amp <- dfm_combined |>  # amplitude measurements
    dplyr::filter(Variant %in% c("WT", "H170A", "H170D", "H170E") & Series == "Series 1" & LPR %in% c("LPR5", "LPR10") &  Measurement == "Maximum amplitude") |>
    dplyr::select(Variant, Condition, LPR, Value)
  dfm_fig5_k <- dfm_combined |>  # rate measurements
    dplyr::filter(Variant %in% c("WT", "H170A", "H170D", "H170E") & Series == "Series 1" & LPR %in% c("LPR5", "LPR10") &  Measurement == "k") |>
    dplyr::select(Variant, Condition, LPR, Value)
  
  return(list(dfm_wt_amp=dfm_wt_amp,
              dfm_wt_k=dfm_wt_k,
              dfm_fig2_amp=dfm_fig2_amp,
              dfm_fig2_k=dfm_fig2_k,
              dfm_fig3_amp=dfm_fig3_amp,
              dfm_fig3_k=dfm_fig3_k,
              dfm_fig5_amp=dfm_fig5_amp,
              dfm_fig5_k=dfm_fig5_k))  
}

load_dataframes_old <- function() {
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