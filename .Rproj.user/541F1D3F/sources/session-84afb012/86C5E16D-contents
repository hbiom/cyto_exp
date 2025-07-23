# install.packages("BiocManager")
# BiocManager::install("flowCore")
# BiocManager::install("flowViz")
# BiocManager::install("flowWorkspace")
# BiocManager::install("openCyto")
# BiocManager::install("ggcyto")


# Load required libraries
library(flowCore)      # Core functions for reading and processing FCS data
library(openCyto)      # Automated gating tools
library(flowWorkspace) # GatingSet manipulation
library(ggcyto)        # Visualization of flow cytometry data
library(gridExtra)     # Arranging plots
library(ggplot2)       # General plotting tools

# Gating template generator for openCyto
# Save this as e.g. "generate_template.R"

# Set output path
template_path <- "gating_template.csv"

gating_steps <- data.frame(
  alias = c("nonDebris", "singlets", "CD45", "CD3", "Tcell"),
  pop = c("+", "+", "+", "+", "+"),
  parent = c("root", "nonDebris", "singlets", "CD45", "CD3"),
  dims = c("FSC-A", "FSC-A,FSC-W", "APC-A", "FITC-A", "PerCP-Cy5-5-A,APC-Cy7-A"),
  gating_method = c("mindensity", "singletGate", "mindensity", "mindensity", "flowClust.2d"),
  gating_args = c(
    "", 
    "", 
    "", 
    "", 
    "K=4"
  ),
  collapseDataForGating = c("", "", "", "", ""),
  groupBy = c("", "", "", "", ""),
  preprocessing_method = c("", "", "", "", ""),
  preprocessing_args = c("", "", "", "", ""),
  stringsAsFactors = FALSE
)

# Save the updated gating template
write.csv(gating_steps, "gating_template.csv", row.names = FALSE)


cat("✅ Gating template saved to:", normalizePath(template_path), "\n")

gt_df <- read.csv(template_path, stringsAsFactors = FALSE)
print(gt_df$gating_method)


# -------------------------------
# Step 1: Load your FCS file
# -------------------------------
fcs_file <- "~/Desktop/code/cytometry/FCSfiles/yourDATA/Specimen_001_Tube_001_001.fcs"
fr <- read.FCS(fcs_file)

# -------------------------------
# Step 2: Compensation (if present)
# -------------------------------
spill <- keyword(fr)$`SPILL`
if (!is.null(spill)) {
  fr <- compensate(fr, spill)
}

# -------------------------------
# Step 3: Logicle Transformation
# -------------------------------
# View available channels (optional)
print(colnames(fr))

# Update these if needed
channels_to_transform <- c("FSC-A", "SSC-A", "FITC-A", "APC-A", "PerCP-Cy5-5-A", "APC-Cy7-A")
logicle_tf <- estimateLogicle(fr, channels = channels_to_transform)
fr_trans <- transform(fr, logicle_tf)

# -------------------------------
# Step 4: Create GatingSet
# -------------------------------
fs <- flowSet(list(fr_trans))
gs <- GatingSet(fs)

# -------------------------------
# Step 5: Load and apply the gating template
# -------------------------------

gt <- gatingTemplate("gating_template.csv")
gt_gating(gt, gs)

# -------------------------------
# Step 6: Visualize each gate
# -------------------------------
# Create individual ggplots
autoplot(gs, "nonDebris", x = "SSC-A", y = "FSC-A") + ggtitle("Non-debris gate")
autoplot(gs, "singlets", x = "FSC-A", y = "FSC-W") + ggtitle("Singlets gate")
autoplot(gs, "CD45", x = "APC-A") + ggtitle("CD45+ gate")
autoplot(gs, "CD3", x = "FITC-A") + ggtitle("CD3+ gate")
autoplot(gs, "Tcell", x = "PerCP-Cy5-5-A", y = "APC-Cy7-A") + ggtitle("CD4/CD8 gate")

library(ggcyto)

autoplot(gs[[1]])

# Gating strategy
plot(gt)






# Opencyto
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003806
















