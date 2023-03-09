
## Morphological feature extractions:
https://github.com/berenslab/MorphoPy
Module description
Important: MorphoPy requires the soma to be one single point. If more than 3 somatic points are present in the reconstruction file they will be automatically collapsed to the centroid of their convex hull on file loading. If the soma is described by 2 to 3 points they will be automatically collapsed to their mean (also see utils.get_standardized_swc).

MorphoPy currently only supports neurites that connect back to the soma. This means, axons that emerge from dendritic structures can not be handled.

A neuron is represented as a directed acyclic graph with node attributes id, x-, y-, z- position, radius and type_id (soma: 1, axon: 2, dendrite: 3, apical dendrite: 4), and with edge attributes path_length and euclidean_dist. Positions, radius and length mesaures are assumed to be given in microns.

Node and edge attributes

Fig. 1: Node and edge attributes associated with each neuron graph.

All data is stored in the tidy data format.

Please also refer to our tutorial and the documentation.

Density maps
Density maps are marginal histograms over the neural mass. MorphoPy allows you to create density maps of different projections through the function compute_density_maps(). Per default it computes x, y, z, xy, xz and yz density maps from the point cloud of the original reconstruction. The point cloud is constructed through resampling along all neurites with a default distance of 1 micron. The resulting point cloud is then binned into bins of 20 microns and smoothed using Gaussian smoothing with std of 1.

However, you can customize all these parameters by passing a config file to the function (see above).

Morphometric statistics
MorphoPy offers a default selection of 28 single-valued morphometric statistics, namely:

number of branch points
width (x-extent), depth (y-extent), height (z-extent)
number of tips
number of neurites extending from the soma directly (stems)
the total path length (in microns)
average and maximal radius thickness (with the soma excluded)
total surface and volume
maximal path distance to the soma
maximal branch order
maximal, min and median path angle
average soma exit angle
maximal path length of a segment
median intermediate and median terminal segment length
log of max, min and median tortuosity across all edges (= path length/euclidean length)
max, min and average branch angle
maximal branching degree (with soma excluded)
weighted proportional sum of absolute deviations as a measure of tree asymmetry (for more details see https://www.sciencedirect.com/science/article/pii/0165027086901196)
Morphometric statistics that can be queried.

Fig. 2: Explanatory schematic of the morphometric statistics that can be computed on all nodes. Left: distance measures, Right: angles.
