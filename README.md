# Improving nichefinding 
The central idea of this repository is to borrow some of the concepts from [CellCharter](https://github.com/CSOgroup/cellcharter) and build on that to weight the contribution of neighouring cells by the proximity of cells to the cell in question. 

# Some ideas
Choosing the decay scale (what to tune): If our coordinates are in microns, start with sigma (or λ) around 1–2 cell diameters (e.g. 20–40 µm for many tissues, but depends on platform + segmentation). A nice heuristic: set sigma to the **median of your nonzero neighbor distances**.