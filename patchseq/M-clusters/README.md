
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

### stems length calculations:

## Soma feature extraction and sholl analysis with Neurolucida Explorer:
### Soma feature extraction:
https://www.mbfbioscience.com/help/neurolucida_explorer/Content/Analyze/BranchedStructure/neuronSumm.htm

![image](https://user-images.githubusercontent.com/42681557/224106730-f3319249-09d6-4b4d-b5ed-82727c8caec8.png)

- Perimeter

*Length of the contour representing the cell body.*

- Area

*The 2-dimensional cross-sectional area contained within the boundary of the cell body.*

- Feret maximum and Feret minimum

![image](https://user-images.githubusercontent.com/42681557/224108051-35341cb7-18c3-4a49-b50f-2200f2bd05ff.png)

*Feret maximum and Feret minimum refer to the largest and smallest dimensions of a contour, respectively, as if a caliper were used for measurement. The two measurements are independent of one another and not necessarily at right angles to each other.Feret maximum and Feret minimum are reported for each contour drawn.*

- Aspect ratio

*Aspect ratio = Feret Maximum ÷ Feret Minimum*

- Compactness

![image](https://user-images.githubusercontent.com/42681557/224108350-87191a80-7411-46b0-a597-25c55bf3a427.png)

*The range of values is 0 to 1. A circle is the most compact shape (compactness = 1).A square has a compactness of 0.8.*

- Convexity

*Convexity = [Convex Perimeter] / [Perimeter]*

*A completely convex object does not have indentations, and has a convexity value of 1 (e.g., circles, ellipses, and squares). Concave objects have convexity values less than 1. Contours with low convexity have a large boundary between inside and outside areas.*

- Form factor

![image](https://user-images.githubusercontent.com/42681557/224108768-310379eb-e42b-4ccd-ae82-0f8142fd0eec.png)

*The form factor differs from the compactness by considering the complexity of the perimeter of the object. For example, a circle with a smooth perimeter has a compactness of 1 and a form factor of 1. If the smooth perimeter is replaced with a finely jagged edge (like a cell covered in microvilli), the compactness is still near 1, but the form factor is much smaller since the perimeter is lengthened considerably.*


- Roundness

*Roundness = [Compactness]^2 Use to differentiate objects that have small compactness values.*

- Solidity

*Solidity = [Area] / [Convex Area]*

### Sholl analysis:

https://www.mbfbioscience.com/help/neurolucida_explorer/Content/Analyze/Sholl.htm

