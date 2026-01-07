# Metric-Aware Neighborhood Aggregation (MANA)
## Improving compartment identification 
The central idea of this repository is to borrow some of the concepts from [CellCharter](https://github.com/CSOgroup/cellcharter) and extend them with a metric-aware, continuous neighborhood aggregation strategy. Instead of treating all neighboring cells equally within discrete hop-based rings, MANA weights the contribution of each neighboring cell by its physical distance to the cell of interest, such that closer cells contribute more strongly than more distant ones.

This results in a spatial context representation that is continuous, distance-aware, and smoothly varying across tissue space, making it better suited for capturing spatial gradients, microenvironmental transitions, and lesion-associated tissue states in high-resolution spatial transcriptomics data (e.g. Xenium, Visium, etc.).

Conceptually, MANA replaces discrete, ring-based neighborhood pooling with a kernel-weighted, diffusion-like aggregation of cellular representations.