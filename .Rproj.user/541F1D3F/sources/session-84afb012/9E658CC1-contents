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

# --------------------
# 1. Load a single FCS file
# --------------------
# Set your file path here:
fcs_file <- "~/Desktop/code/cytometry/FCSfiles/yourDATA/Specimen_001_Tube_001_001.fcs"

# Read the file into a flowFrame
fr <- read.FCS(fcs_file)

# View metadata and channels (optional)
print(keyword(fr))                    # Shows FCS header keywords
print(pData(parameters(fr)))         # Shows channel names and descriptions
print(colnames(exprs(fr)))           # Column names used for plotting/gating

# --------------------
# 2. Compensation (if available in the FCS)
# --------------------
spill <- keyword(fr)$`SPILL`  # or use $`$SPILLOVER` depending on file
if (!is.null(spill)) {
  fr <- compensate(fr, spill)
}

# --------------------
# 3. Transformation (Logicle)
# --------------------
# Define channels to transform (update based on your file)
channels_to_transform <- c("FSC-A", "SSC-A", "FITC-A", "APC-A", "PerCP-Cy5-5-A", "APC-Cy7-A")

# Estimate and apply Logicle transform
logicle_tf <- estimateLogicle(fr, channels = channels_to_transform)
fr_trans <- transform(fr, logicle_tf)

# --------------------
# 4. Gating Steps
# --------------------
# Debris exclusion
gate_debris <- gate_mindensity(fr_trans, channel = "FSC-A")

# Singlet gating (FSC-A vs FSC-W)
gate_singlets <- gate_singlet(fr_trans, channels = c("FSC-A", "FSC-W"))

# CD45+ gating (APC-A)
gate_CD45 <- gate_mindensity(fr_trans, channel = "APC-A")

# CD3+ gating (FITC-A)
gate_CD3 <- gate_mindensity(fr_trans, channel = "FITC-A")

# CD4/CD8 quad gate (PerCP-Cy5-5-A vs APC-Cy7-A)
quad_gate <- gate_quad_tmix(fr_trans, channels = c("PerCP-Cy5-5-A", "APC-Cy7-A"), K = 3)

# --------------------
# 5. Plotting
# --------------------
# Plot each gating step
p1 <- autoplot(fr_trans, x = "FSC-A", y = "SSC-A") + ggcyto::geom_gate(gate_debris) + ggtitle("Non-debris")
p2 <- autoplot(fr_trans, x = "FSC-A", y = "FSC-W") + ggcyto::geom_gate(gate_singlets) + ggtitle("Singlets")
p3 <- autoplot(fr_trans, x = "APC-A") + ggcyto::geom_gate(gate_CD45) + ggtitle("CD45+")
p4 <- autoplot(fr_trans, x = "FITC-A") + ggcyto::geom_gate(gate_CD3) + ggtitle("CD3+")
p5 <- autoplot(fr_trans, x = "PerCP-Cy5-5-A", y = "APC-Cy7-A") + ggcyto::geom_gate(quad_gate) + ggtitle("CD4/CD8 subsets")

# Arrange and save plots
g <- grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
ggsave("single_sample_plots.png", g, width = 12, height = 8)

