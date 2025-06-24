#' Combine results from different LIANA methods
#'
#' @param ... CSV file paths of each method's output
#' @return A data frame of interactions common to all methods
#' @export
combine_results <- function(cellchat, cellphonedb, connectome, logfc, natmi, sca) {
  cellchat <- read.csv(cellchat) %>% rename_with(~paste0("cellchat_", .), -c(patient_id, source, target, ligand, receptor))
  cellphonedb <- read.csv(cellphonedb) %>% rename_with(~paste0("cellphonedb_", .), -c(patient_id, source, target, ligand, receptor))
  connectome <- read.csv(connectome) %>% rename_with(~paste0("connectome_", .), -c(patient_id, source, target, ligand, receptor))
  logfc <- read.csv(logfc) %>% rename_with(~paste0("logfc_", .), -c(patient_id, source, target, ligand, receptor))
  natmi <- read.csv(natmi) %>% rename_with(~paste0("natmi_", .), -c(patient_id, source, target, ligand, receptor))
  sca <- read.csv(sca) %>% rename_with(~paste0("sca_", .), -c(patient_id, source, target, ligand, receptor))

  common <- cellchat %>%
    inner_join(cellphonedb, by = c("patient_id", "source", "target", "ligand", "receptor")) %>%
    inner_join(connectome, by = c("patient_id", "source", "target", "ligand", "receptor")) %>%
    inner_join(logfc, by = c("patient_id", "source", "target", "ligand", "receptor")) %>%
    inner_join(natmi, by = c("patient_id", "source", "target", "ligand", "receptor")) %>%
    inner_join(sca, by = c("patient_id", "source", "target", "ligand", "receptor"))

  return(common)
}