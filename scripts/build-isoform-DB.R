# Load necessary libraries
library(dplyr)
library(CellChat)
library(biomaRt)
library(SingleCellExperiment)

# ---- 1. Import new ENST - ENST database ----

my_isoform_ppi_df <- read.delim("transcript_specific_interactions.tsv", 
                                header = TRUE,
                                sep = "\t",
                                stringsAsFactors = FALSE)

# ---- 2. Define Helper Functions ----

parse_complex <- function(name) {
  if (name %in% rownames(complex_db)) {
    subunits <- na.omit(unlist(complex_db[name, ]))
    subunits[subunits != ""]
  } else character(0)
}

expand_cellchat_with_complex_handling <- function(cellchat_row) {
  
  ligand_symbol <- cellchat_row$ligand
  receptor_symbol <- cellchat_row$receptor
  
  ligand_is_complex <- ligand_symbol %in% complex_ligands
  receptor_is_complex <- receptor_symbol %in% complex_receptors
  
  # --- HELPER: Get all isoforms (Unfiltered) ---
  get_isoforms_from_symbol <- function(gene_sym) {
    # Check Ligand column
    lig_ensts <- my_isoform_ppi_df %>% 
      filter(ligand_gene == gene_sym) %>% 
      pull(ligand_transcript)
    
    # Check Receptor column
    rec_ensts <- my_isoform_ppi_df %>% 
      filter(receptor_gene == gene_sym) %>% 
      pull(receptor_transcript)
    
    all_ensts <- unique(c(lig_ensts, rec_ensts))
    
    if (length(all_ensts) > 0) {
      return(all_ensts)
    } else {
      # Fallback to symbol if no isoforms found (prevents crash)
      return(gene_sym)
    }
  }

  # --- HELPER: Validate Interaction ---
  # Checks if ANY isoform on the left interacts with ANY isoform on the right
  validate_interaction <- function(ligand_isoforms, receptor_isoforms) {
    # subset PPI to relevant transcripts to speed up check
    relevant_ppi <- my_isoform_ppi_df %>%
      filter(
        (ligand_transcript %in% ligand_isoforms & receptor_transcript %in% receptor_isoforms) |
        (receptor_transcript %in% ligand_isoforms & ligand_transcript %in% receptor_isoforms)
      )
    return(nrow(relevant_ppi) > 0)
  }

  # =========================================================================
  # CASE 1: Simple -> Simple
  # =========================================================================
  if (!ligand_is_complex && !receptor_is_complex) {
    l_isoforms <- get_isoforms_from_symbol(ligand_symbol)
    r_isoforms <- get_isoforms_from_symbol(receptor_symbol)
    
    # Expand Grid
    combos <- expand.grid(L = l_isoforms, R = r_isoforms, stringsAsFactors = FALSE)
    
    # Filter valid interactions
    valid_rows <- do.call(rbind, lapply(1:nrow(combos), function(i) {
      l_curr <- combos$L[i]
      r_curr <- combos$R[i]
      
      if (validate_interaction(l_curr, r_curr)) {
        new_row <- cellchat_row
        new_row$ligand <- l_curr
        new_row$receptor <- r_curr
        new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
        new_row$interaction_name_2 <- paste(l_curr, r_curr, sep = " - ")
        return(new_row)
      } else {
        return(NULL)
      }
    }))
    
    return(list(interactions = valid_rows, complex_mapping = data.frame()))
  }

  # =========================================================================
  # CASE 2: Simple -> Complex
  # =========================================================================
  if (!ligand_is_complex && receptor_is_complex) {
    l_isoforms <- get_isoforms_from_symbol(ligand_symbol)
    
    # 1. Get subunits
    rec_subunits <- parse_complex(receptor_symbol)
    
    # 2. Get all isoforms for each subunit (No filtering yet!)
    rec_subunit_isoforms <- lapply(rec_subunits, get_isoforms_from_symbol)
    names(rec_subunit_isoforms) <- rec_subunits
    
    # 3. Expand all possible receptor complex combinations
    # This assumes independent assembly (we don't check intra-complex binding here)
    rec_combos <- expand.grid(rec_subunit_isoforms, stringsAsFactors = FALSE)
    
    # 4. Iterate and validate against Ligand
    interactions <- list()
    complex_map_list <- list()
    
    if(nrow(rec_combos) > 0) {
      for (i in 1:nrow(rec_combos)) {
        # Construct complex name
        curr_rec_sub_isos <- as.character(rec_combos[i, ])
        complex_name <- paste(curr_rec_sub_isos, collapse = "_")
        
        # Does this specific complex combination interact with ANY ligand isoform?
        # We pass the vector of subunits. If ANY subunit binds the ligand, it counts.
        if (validate_interaction(l_isoforms, curr_rec_sub_isos)) {
          
          # Now find WHICH ligand isoform specifically binds (could be multiple)
          for (l_iso in l_isoforms) {
            if (validate_interaction(l_iso, curr_rec_sub_isos)) {
              new_row <- cellchat_row
              new_row$ligand <- l_iso
              new_row$receptor <- complex_name
              new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
              new_row$interaction_name_2 <- paste(l_iso, complex_name, sep = " - ")
              interactions[[length(interactions) + 1]] <- new_row
            }
          }
          
          # Save Mapping
          map_row <- data.frame(complex_name = complex_name, stringsAsFactors = FALSE)
          for (k in seq_along(curr_rec_sub_isos)) {
            map_row[[paste0("subunit_", k)]] <- curr_rec_sub_isos[k]
          }
          complex_map_list[[length(complex_map_list) + 1]] <- map_row
        }
      }
    }
    
    return(list(
      interactions = do.call(rbind, interactions), 
      complex_mapping = do.call(rbind, complex_map_list)
    ))
  }

  # =========================================================================
  # CASE 3: Complex -> Simple
  # =========================================================================
  if (ligand_is_complex && !receptor_is_complex) {
    r_isoforms <- get_isoforms_from_symbol(receptor_symbol)
    
    lig_subunits <- parse_complex(ligand_symbol)
    lig_subunit_isoforms <- lapply(lig_subunits, get_isoforms_from_symbol)
    lig_combos <- expand.grid(lig_subunit_isoforms, stringsAsFactors = FALSE)
    
    interactions <- list()
    complex_map_list <- list()
    
    if(nrow(lig_combos) > 0) {
      for (i in 1:nrow(lig_combos)) {
        curr_lig_sub_isos <- as.character(lig_combos[i, ])
        complex_name <- paste(curr_lig_sub_isos, collapse = "_")
        
        if (validate_interaction(curr_lig_sub_isos, r_isoforms)) {
          
          for (r_iso in r_isoforms) {
            if (validate_interaction(curr_lig_sub_isos, r_iso)) {
              new_row <- cellchat_row
              new_row$ligand <- complex_name
              new_row$receptor <- r_iso
              new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
              new_row$interaction_name_2 <- paste(complex_name, r_iso, sep = " - ")
              interactions[[length(interactions) + 1]] <- new_row
            }
          }
          
          map_row <- data.frame(complex_name = complex_name, stringsAsFactors = FALSE)
          for (k in seq_along(curr_lig_sub_isos)) {
            map_row[[paste0("subunit_", k)]] <- curr_lig_sub_isos[k]
          }
          complex_map_list[[length(complex_map_list) + 1]] <- map_row
        }
      }
    }
    
    return(list(
      interactions = do.call(rbind, interactions), 
      complex_mapping = do.call(rbind, complex_map_list)
    ))
  }

  # =========================================================================
  # CASE 4: Complex -> Complex
  # =========================================================================
  if (ligand_is_complex && receptor_is_complex) {
    # 1. Expand Ligand
    lig_subunits <- parse_complex(ligand_symbol)
    lig_subunit_isoforms <- lapply(lig_subunits, get_isoforms_from_symbol)
    lig_combos <- expand.grid(lig_subunit_isoforms, stringsAsFactors = FALSE)
    
    # 2. Expand Receptor
    rec_subunits <- parse_complex(receptor_symbol)
    rec_subunit_isoforms <- lapply(rec_subunits, get_isoforms_from_symbol)
    rec_combos <- expand.grid(rec_subunit_isoforms, stringsAsFactors = FALSE)
    
    interactions <- list()
    complex_map_list <- list()
    
    # Loop over Ligand Combos
    if(nrow(lig_combos) > 0 && nrow(rec_combos) > 0) {
      for (i in 1:nrow(lig_combos)) {
        curr_lig_sub_isos <- as.character(lig_combos[i, ])
        l_complex_name <- paste(curr_lig_sub_isos, collapse = "_")
        
        # Loop over Receptor Combos
        for (j in 1:nrow(rec_combos)) {
          curr_rec_sub_isos <- as.character(rec_combos[j, ])
          r_complex_name <- paste(curr_rec_sub_isos, collapse = "_")
          
          # Check if ANY piece of Ligand Complex binds ANY piece of Receptor Complex
          if (validate_interaction(curr_lig_sub_isos, curr_rec_sub_isos)) {
            new_row <- cellchat_row
            new_row$ligand <- l_complex_name
            new_row$receptor <- r_complex_name
            new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
            new_row$interaction_name_2 <- paste(l_complex_name, r_complex_name, sep = " - ")
            interactions[[length(interactions) + 1]] <- new_row
            
            # Save Mapping (Ligand)
            l_map <- data.frame(complex_name = l_complex_name, stringsAsFactors = FALSE)
            for (k in seq_along(curr_lig_sub_isos)) l_map[[paste0("subunit_", k)]] <- curr_lig_sub_isos[k]
            complex_map_list[[length(complex_map_list) + 1]] <- l_map
            
            # Save Mapping (Receptor)
            r_map <- data.frame(complex_name = r_complex_name, stringsAsFactors = FALSE)
            for (k in seq_along(curr_rec_sub_isos)) r_map[[paste0("subunit_", k)]] <- curr_rec_sub_isos[k]
            complex_map_list[[length(complex_map_list) + 1]] <- r_map
          }
        }
      }
    }
    
    return(list(
      interactions = do.call(rbind, interactions), 
      complex_mapping = unique(do.call(rbind, complex_map_list))
    ))
  }
}

# ---- 3. Prepare Base Database ----

db.human <- CellChatDB.human
complex_db <- db.human[["complex"]]
all_complex_names <- rownames(complex_db)
complex_ligands <- unique(db.human$interaction$ligand[db.human$interaction$ligand %in% all_complex_names])
complex_receptors <- unique(db.human$interaction$receptor[db.human$interaction$receptor %in% all_complex_names])

# ---- 4. Filter and Expand ----

genes_in_my_ppi <- unique(c(my_isoform_ppi_df$ligand_gene, my_isoform_ppi_df$receptor_gene))
base_interactions <- db.human$interaction

filtered_interactions <- base_interactions %>%
  filter((ligand %in% genes_in_my_ppi) | (receptor %in% genes_in_my_ppi))

cat("Expanding interactions...\n")
interaction_list <- split(filtered_interactions, seq(nrow(filtered_interactions)))
results_list <- lapply(interaction_list, expand_cellchat_with_complex_handling)

# ---- 5. Assemble New Database ----

isoform_interactions <- dplyr::bind_rows(lapply(results_list, `[[`, "interactions"))
isoform_complexes <- dplyr::bind_rows(lapply(results_list, `[[`, "complex_mapping"))

# Logic for cleaning isoform_complexes_final and creating myIsoformCellChatDB...
# (Omitted repetitive formatting for space)

# ---- 6. Update geneInfo and Cofactors ----

# Create combined map for GeneInfo
ligand_map <- data.frame(Symbol = my_isoform_ppi_df$ligand_transcript, gene.symbol = my_isoform_ppi_df$ligand_gene)
receptor_map <- data.frame(Symbol = my_isoform_ppi_df$receptor_transcript, gene.symbol = my_isoform_ppi_df$receptor_gene)
new_geneInfo <- unique(rbind(ligand_map, receptor_map))
rownames(new_geneInfo) <- new_geneInfo$Symbol

# Update cofactors in GeneInfo
cofactor_genes <- unique(na.omit(unlist(db.human$cofactor)))
new_geneInfo_rows <- data.frame(Symbol = cofactor_genes, gene.symbol = cofactor_genes, row.names = cofactor_genes)
new_geneInfo <- rbind(new_geneInfo, new_geneInfo_rows)

# Save result
saveRDS(isoform_interactions, file = "../data/myIsoformCellChatDB.rds")

# ---- 7. Aggregate Transcript Counts for Cofactors ----

# Requires BioMart connection
ensembl_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "[https://useast.ensembl.org](https://useast.ensembl.org)")

tx_map <- getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'),
                filters = 'hgnc_symbol',
                values = cofactor_genes,
                mart = ensembl_mart)

# Aggregation logic loop (as provided in your qmd)...
# Updates 'sce_modified' with summed transcript counts for gene-level cofactors.
