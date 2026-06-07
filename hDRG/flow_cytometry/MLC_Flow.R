## MLC flow analysis

library(ggplot2)
library(dplyr)

# define custom functions
saveThePlot <- function(path, plot = last_plot(), w = 8, h = 6, u = "in") {
  for(fileformat in c("pdf", "tiff")) {
    ggsave(plot = plot, filename = paste0(path, "_", Sys.Date(), ".", fileformat), 
           device = fileformat, width = w, height = h, units = u)   
  }
}

# set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

flow_p2y12 <- read.csv("mlc_p2y12.csv") |> filter(cell_type != "CD45loCD11b+")
#flow_p2y12$cell_type <- factor(flow_p2y12$cell_type,c("CD45+CD11b+",  "CD45loCD11b+", "SGCs", "CD45+CD11b-", "CD45-"))
flow_p2y12$cell_type <- factor(flow_p2y12$cell_type,c("CD45+CD11b+", "SGCs", "CD45+CD11b-", "CD45-"))

# percent expressing P2Y12

ggplot(flow_p2y12, aes(
  x = cell_type,
  y = pc_p2y12,
  fill = cell_type
)) + geom_col(color="black") + 
  scale_y_continuous(limits=c(0,100)) + 
  scale_fill_manual(values=c("#f6222e", "gray", "#006838", "#21409a")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))

saveThePlot(path=file.path(paste0("MLC_P2Y12")), w=4, h=3.5)

# pie chart for cell types
flow_data <- read.csv("mlc_live.csv") |> filter(cell_type != "Live")

hole_size <- 1.5

ggplot(flow_data, aes(x = hole_size, y = pc_live, fill = cell_type)) +
  geom_col(color = "black") +
  geom_text(aes(label = round(pc_live, 1)),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y", start = 0) +
  xlim(c(0.2, hole_size + 0.5)) +
  scale_fill_manual(values=c("#f6222e", "#31A1B3", "#CCB22B", "gray", "#006838")) + 
  theme_void()

saveThePlot(path=file.path(paste0("Pie_PC_Live")), w=5, h=3.5)
