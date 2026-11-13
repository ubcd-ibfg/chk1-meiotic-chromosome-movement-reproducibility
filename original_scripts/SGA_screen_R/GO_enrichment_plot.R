##SCRIPT FOR GO ENRICHMENT PLOT

library(ggplot2)

data <- data.frame(
  Process = c("Gene expression", "RNA metabolic process",
              "Macromolecular complex assembly", "Macromolecule metabolic process",
              "Macromolecular complex subunit organization",
              "Cellular macromolecule metabolic process"),
  p_value = c(0.0000641, 0.0011759, 0.0240239, 0.0249016, 0.029001, 0.0333114)
)

data$log_p.value <- -log10(data$p_value)

ggplot(data, aes(x = log_p, y = reorder(Process, log_p.value))) +
  geom_point(aes(size = log_p.value, color = log_p.value)) +
  scale_color_gradient(low = "#CAECEA", high = "blue") +
  theme_bw() +
  labs(
    title = "GO Terms",
    x = expression(-log[10](p-value)),
    y = NULL,
    size = "log10(p-value)"
  )
