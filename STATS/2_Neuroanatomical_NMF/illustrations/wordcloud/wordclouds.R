################################################################################
# Written by Clément Guichet, PhD Student
# LPNC - CNRS UMR 5105
# 2024

################################################################################

# Word clouds ----

NMF_networks_composition <- rio::import("E:/Research_Projects/MiddleAge_LARA/STATS/2_Neuroanatomical_NMF/output/NMF_composition/NMF_networks_composition.xlsx")

wordcloud_data <- NMF_networks_composition[, c(1, seq(2, 27, 5))]
wordcloud_data$track_name...1 <- gsub("_", " ", wordcloud_data$track_name...1)
wordcloud_data <- wordcloud_data %>%
  t() %>%
  as.data.frame()
rownames(wordcloud_data) <- c("Bundle", "Network 4", "Network 5", "Network 7", "Network 11", "Network 14", "Network 16")


wordcloud_plot <- function(i) {
  library(wordcloud2)
  wordcloud_data_network <<- wordcloud_data[c(1, i + 1), ] %>%
    janitor::row_to_names(1) %>%
    t() %>%
    as.data.frame() %>%
    mutate_at(vars(everything()), funs(as.numeric(.))) %>%
    subset(.[, 1] >= .95) %>%
    rownames_to_column(var = "Bundle") %>%
    dplyr::rename(value = colnames(.)[2]) %>%
    mutate(value = value / 10) %>%
    arrange(desc(value))

  # Make the plot
  wordcloud2(wordcloud_data_network,
    size = 1 / 2, minRotation = -pi / 12, maxRotation = pi / 12, rotateRatio = .2,
    backgroundColor = "white", color = "#FF8243", shuffle = F, gridSize = 20
  )
}

# Net 16
wordcloud_plot(6)
# Net 5
wordcloud_plot(2)
# Net 11
wordcloud_plot(4)
# Net 4
wordcloud_plot(1)
# Net 14
wordcloud_plot(5)
# Net 7
wordcloud_plot(3)
