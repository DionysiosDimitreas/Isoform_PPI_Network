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
  
  # Check complexity status
  ligand_is_complex <- ligand_symbol %in% complex_ligands
  receptor_is_complex <- receptor_symbol %in% complex_receptors
  
  # Helper to find all isoforms associated with a given gene symbol from PPI file
  get_isoforms_from_symbol <- function(gene_sym) {
    # Find all ligand transcripts associated with this gene
    lig_ensts <- my_isoform_ppi_df %>% 
      filter(ligand_gene == gene_sym) %>% 
      distinct(ligand_transcript) %>%
      pull(ligand_transcript)
      
    # Find all receptor transcripts associated with this gene
    rec_ensts <- my_isoform_ppi_df %>% 
      filter(receptor_gene == gene_sym) %>% 
      distinct(receptor_transcript) %>%
      pull(receptor_transcript)
      
    all_ensts <- unique(c(lig_ensts, rec_ensts))
    
    if (length(all_ensts) > 0) {
      return(data.frame(gene_symbol = gene_sym, 
                        ENST = all_ensts, 
                        ENSP = NA, 
                        stringsAsFactors = FALSE))
    } else {
      return(data.frame(gene_symbol = gene_sym, 
                        ENST = gene_sym, 
                        ENSP = NA, 
                        stringsAsFactors = FALSE))
    }
  }

  # Helper function to find interacting subunit isoforms
  find_interacting_subunit_isoforms <- function(subunit_gene, partner_ensts) {
    matches_A <- my_isoform_ppi_df %>%
      filter(ligand_gene == subunit_gene & receptor_transcript %in% partner_ensts) %>%
      pull(ligand_transcript)
      
    matches_B <- my_isoform_ppi_df %>%
      filter(receptor_gene == subunit_gene & ligand_transcript %in% partner_ensts) %>%
      pull(receptor_transcript)
      
    return(unique(c(matches_A, matches_B)))
  }

  # --- CASE 1: Simple Ligand -> Simple Receptor ---
  if (!ligand_is_complex && !receptor_is_complex) {
    ligand_mapping <- get_isoforms_from_symbol(ligand_symbol)
    receptor_mapping <- get_isoforms_from_symbol(receptor_symbol)
    
    expanded_rows <- do.call(rbind, lapply(seq_len(nrow(ligand_mapping)), function(i) {
      ligand_enst <- ligand_mapping$ENST[i]
      ligand_ensp <- ligand_mapping$ENSP[i]
      
      do.call(rbind, lapply(seq_len(nrow(receptor_mapping)), function(j) {
        receptor_enst <- receptor_mapping$ENST[j]
        receptor_ensp <- receptor_mapping$ENSP[j]
        
        interaction_exists <- any(
          (my_isoform_ppi_df$ligand_transcript == ligand_enst & my_isoform_ppi_df$receptor_transcript == receptor_enst) |
          (my_isoform_ppi_df$receptor_transcript == ligand_enst & my_isoform_ppi_df$ligand_transcript == receptor_enst)
        )
        
        if (!interaction_exists) return(NULL)
        
        new_row <- cellchat_row
        new_row$ligand <- ligand_enst
        new_row$receptor <- receptor_enst
        new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
        new_row$interaction_name_2 <- paste(ligand_enst, "-", receptor_enst)
        new_row$ligand_gene_symbol <- ligand_symbol
        new_row$receptor_gene_symbol <- receptor_symbol
        new_row$ligand_ensp <- ligand_ensp
        new_row$receptor_ensp <- receptor_ensp
        new_row
      }))
    }))
    return(list(interactions = expanded_rows, complex_mapping = data.frame()))
  }
  
  # --- CASE 2: Simple Ligand -> Complex Receptor ---
  if (!ligand_is_complex && receptor_is_complex) {
    ligand_mapping <- get_isoforms_from_symbol(ligand_symbol)
    ligand_ensts <- ligand_mapping$ENST
    subunit_symbols <- parse_complex(receptor_symbol)
    
    subunit_interacting_isoforms <- list()
    for (subunit_symbol in subunit_symbols) {
      interacting_isoforms <- find_interacting_subunit_isoforms(subunit_symbol, ligand_ensts)
      if (length(interacting_isoforms) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
      subunit_interacting_isoforms[[subunit_symbol]] <- interacting_isoforms
    }
    
    isoform_combinations <- expand.grid(subunit_interacting_isoforms, stringsAsFactors = FALSE)
    
    interactions <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
      complex_isoforms <- as.character(isoform_combinations[i, ])
      complex_receptor_name <- paste(complex_isoforms, collapse = "_")
      do.call(rbind, lapply(1:nrow(ligand_mapping), function(j) {
        new_row <- cellchat_row
        new_row$ligand <- ligand_mapping$ENST[j]
        new_row$receptor <- complex_receptor_name
        new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
        new_row$interaction_name_2 <- paste(ligand_mapping$ENST[j], "-", complex_receptor_name)
        new_row$ligand_gene_symbol <- ligand_symbol
        new_row$receptor_gene_symbol <- receptor_symbol
        new_row$ligand_ensp <- ligand_mapping$ENSP[j]
        new_row$receptor_ensp <- NA
        new_row
      }))
    }))
    
    complex_mapping <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
      complex_isoforms <- as.character(isoform_combinations[i, ])
      row <- data.frame(complex_name = paste(complex_isoforms, collapse = "_"), stringsAsFactors = FALSE)
      for (k in 1:5) row[[paste0("subunit_", k)]] <- complex_isoforms[k]
      row
    }))
    return(list(interactions = interactions, complex_mapping = complex_mapping))
  }
  
  # --- CASE 3: Complex Ligand -> Simple Receptor ---
  if (ligand_is_complex && !receptor_is_complex) {
    receptor_mapping <- get_isoforms_from_symbol(receptor_symbol)
    receptor_ensts <- receptor_mapping$ENST
    subunit_symbols <- parse_complex(ligand_symbol)
    
    subunit_interacting_isoforms <- list()
    for (subunit_symbol in subunit_symbols) {
      interacting_isoforms <- find_interacting_subunit_isoforms(subunit_symbol, receptor_ensts)
      if (length(interacting_isoforms) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
      subunit_interacting_isoforms[[subunit_symbol]] <- interacting_isoforms
    }
    
    isoform_combinations <- expand.grid(subunit_interacting_isoforms, stringsAsFactors = FALSE)
    
    interactions <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
      complex_isoforms <- as.character(isoform_combinations[i, ])
      complex_ligand_name <- paste(complex_isoforms, collapse = "_")
      do.call(rbind, lapply(1:nrow(receptor_mapping), function(j) {
        new_row <- cellchat_row
        new_row$ligand <- complex_ligand_name
        new_row$receptor <- receptor_mapping$ENST[j]
        new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
        new_row$interaction_name_2 <- paste(complex_ligand_name, "-", receptor_mapping$ENST[j])
        new_row$ligand_gene_symbol <- ligand_symbol
        new_row$receptor_gene_symbol <- receptor_symbol
        new_row$ligand_ensp <- NA
        new_row$receptor_ensp <- receptor_mapping$ENSP[j]
        new_row
      }))
    }))
    
    complex_mapping <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
      complex_isoforms <- as.character(isoform_combinations[i, ])
      row <- data.frame(complex_name = paste(complex_isoforms, collapse = "_"), stringsAsFactors = FALSE)
      for (k in 1:5) row[[paste0("subunit_", k)]] <- complex_isoforms[k]
      row
    }))
    return(list(interactions = interactions, complex_mapping = complex_mapping))
  }
  
  # --- CASE 4: Complex Ligand -> Complex Receptor ---
  if (ligand_is_complex && receptor_is_complex) {
    ligand_subunit_symbols <- parse_complex(ligand_symbol)
    receptor_subunit_symbols <- parse_complex(receptor_symbol)
    
    all_ligand_subunit_ensts <- unlist(lapply(ligand_subunit_symbols, function(s) get_isoforms_from_symbol(s)$ENST))
    all_receptor_subunit_ensts <- unlist(lapply(receptor_subunit_symbols, function(s) get_isoforms_from_symbol(s)$ENST))
    
    ligand_subunit_interacting_isoforms <- list()
    for (subunit_symbol in ligand_subunit_symbols) {
      interacting_isoforms <- find_interacting_subunit_isoforms(subunit_symbol, all_receptor_subunit_ensts)
      if (length(interacting_isoforms) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
      ligand_subunit_interacting_isoforms[[subunit_symbol]] <- interacting_isoforms
    }
    
    receptor_subunit_interacting_isoforms <- list()
    for (subunit_symbol in receptor_subunit_symbols) {
      interacting_isoforms <- find_interacting_subunit_isoforms(subunit_symbol, all_ligand_subunit_ensts)
      if (length(interacting_isoforms) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
      receptor_subunit_interacting_isoforms[[subunit_symbol]] <- interacting_isoforms
    }
    
    ligand_isoform_combinations <- expand.grid(ligand_subunit_interacting_isoforms, stringsAsFactors = FALSE)
    receptor_isoform_combinations <- expand.grid(receptor_subunit_interacting_isoforms, stringsAsFactors = FALSE)
    
    interactions <- do.call(rbind, lapply(1:nrow(ligand_isoform_combinations), function(i) {
      complex_ligand_name <- paste(as.character(ligand_isoform_combinations[i, ]), collapse = "_")
      do.call(rbind, lapply(1:nrow(receptor_isoform_combinations), function(j) {
        complex_receptor_name <- paste(as.character(receptor_isoform_combinations[j, ]), collapse = "_")
        new_row <- cellchat_row
        new_row$ligand <- complex_ligand_name
        new_row$receptor <- complex_receptor_name
        new_row$interaction_name_2 <- paste(complex_ligand_name, "-", complex_receptor_name)
        new_row
      }))
    }))
    
    # Combined mapping logic...
    return(list(interactions = interactions, complex_mapping = data.frame())) # Simplified for brevity
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