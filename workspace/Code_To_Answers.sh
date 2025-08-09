Answers to questions path, please run full script succesfully first to run those answers.

```{r}
# ---- Answers for q1..q5 ----
# Requires objects 'counts', 'dds', and 'res' created above

# q1: How many sequencing lanes were concatenated to form sample T0hrrep1??
# If you later generate a manifest upstream, read it here instead of hardcoding
sample_lanes <- c(T0hrrep1 = 2L, T0hrrep2 = 2L, T0hrCAFrep1 = 2L, T0hrCAFrep2 = 2L, T6hrrep1 = 2L, T6hrrep2 = 2L, T6hrCAFrep1 = 2L, T6hrCAFrep2 = 2L)
ans_q1 <- unname(sample_lanes["ym"])

# q2: What is the library size (total read counts) for T6hrrep1??
libsize <- colSums(counts)
ans_q2 <- unname(libsize["T6hrCAFrep1"])+ unname(libsize["T6hrCAFrep2"])

# q3: How many genes have nonzero counts in sample T0hrCAFrep1??
ans_q3 <- sum(counts[, "T0hrCAFrep1"] > 0)

# q4: How many genes are upregulated 2-fold (log2FC 1) in T0hrCAF vs. T0hrCAF with FDR < 0.01?
res_df$SYMBOL <- rownames(res_df)
up_rank <- res_df %>%
  dplyr::filter(!is.na(log2FoldChange) & log2FoldChange > 0) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange), padj, pvalue)
ans_q5 <- if (nrow(up_rank) >= 1) up_rank$SYMBOL[3] else NA_character_

# q5: Which gene is ranked 3th by log2 fold change (most upregulated) in T0hrCAF vs. T6hrCAF?
res_df$SYMBOL <- rownames(res_df)
up_rank <- res_df %>%
  dplyr::filter(!is.na(log2FoldChange) & log2FoldChange > 0) %>%
  dplyr::arrange(dplyr::desc(log2FoldChange), padj, pvalue)
ans_q5 <- if (nrow(up_rank) >= 2) up_rank$SYMBOL[3] else NA_character_

answers <- tibble::tibble(
  id = c("q1","q2","q3","q4","q5"),
  answer = c(ans_q1, ans_q2, ans_q3, ans_q4, ans_q5)
)

answers
```
