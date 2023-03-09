
# Morphological feature extractions:
## Dendritic feature extraction with MorphoPy:

https://github.com/berenslab/MorphoPy

#### Morphometric statistics

`MorphoPy` offers a default selection of 28 single-valued morphometric statistics, namely:
- number of branch points
- width (x-extent), depth (y-extent), height (z-extent)
- number of tips
- number of neurites extending from the soma directly (stems)
- the total path length (in microns)
- average and maximal radius thickness (with the soma excluded)
- total surface and volume
- maximal path distance to the soma
- maximal branch order
- maximal, min and median path angle
- average soma exit angle
- maximal path length of a segment
- median intermediate and median terminal segment length
- log of max, min and median tortuosity across all edges (= path length/euclidean length)
- max, min and average branch angle
- maximal branching degree (with soma excluded)
- _weighted proportional sum of absolute deviations_ as a measure of tree asymmetry (for more details see https://www.sciencedirect.com/science/article/pii/0165027086901196)

![Morphometric statistics that can be queried.](https://user-images.githubusercontent.com/520137/80974473-0f4d2380-8e21-11ea-8ce2-acb8153cece4.png)

*Fig. 2: Explanatory schematic of the morphometric statistics that can be computed on all nodes. Left: distance measures, Right: angles.*

## Soma feature extraction and sholl analysis with Neurolucida Explorer:
### Soma feature extraction:
https://www.mbfbioscience.com/help/neurolucida_explorer/Content/Analyze/BranchedStructure/neuronSumm.htm

![image](https://user-images.githubusercontent.com/42681557/224106730-f3319249-09d6-4b4d-b5ed-82727c8caec8.png)

- Perimeter

* Length of the contour representing the cell body.

- Area

* The 2-dimensional cross-sectional area contained within the boundary of the cell body.

- Feret maximum and Feret minimum
- Aspect ratio
- Compactness
- Convexity
- Form factor
- Roundness
- Solidity


https://www.mbfbioscience.com/help/neurolucida_explorer/Content/Analyze/Sholl.htm
### 
