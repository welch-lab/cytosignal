# CytoSignal and VelcoCytoSignal

CytoSignal (VeloCytoSignal) is a tool for detecting cell-cell signaling interactions and their dynamics at single-cell resolution from spatial transcriptomic data.

Nearby cells within tissues communicate through ligand and receptor signaling interactions. These interactions play crucial roles in cell differentiation, tissue homeostasis, immune response, and disease, but have been difficult to study in an unbiased fashion. New spatial transcriptomic protocols provide a tremendous opportunity to systematically detect cell-cell signaling. Some computational tools can detect signaling among cell types using dissociated or spatial data, but no method operates at cellular resolution. In addition, no existing approaches can predict the future states of cell signaling from snapshot data.

Here, we address these limitations by developing two computational tools for detecting cell-cell signaling from spatial transcriptomic data. The first tool is CytoSignal, which performs a nonparametric statistical test to identify which cells within a tissue have significant activity for a particular signaling interaction. CytoSignal considers multi-component interactions and separately models interactions mediated by diffusible vs. contact-dependent molecules. Second, we will develop VeloCytoSignal, a method for predicting the rate of change for a signaling interaction–whether the strength of a signaling interaction is increasing or decreasing at each tissue location. This approach combines RNA velocities for ligands and receptors simultaneously to predict future signaling strength using snapshot data from as little as one timepoint.

Both tools will be incorporated into one user-friendly open-source R package. Our work addresses the field's current need of a robust and scalable tool to detect cell-cell signaling interactions and their dynamics at single-cell resolution from spatial transcriptomic data.