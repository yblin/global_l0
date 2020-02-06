# Global L0 algorithm for regularity-constrained plane fitting

## Introduction
We introduce a method based on energy minimization to reconstruct the planes consistent with our constraint model. The proposed algorithm is efficient, easily to understand, and simple to implement. The experimental results show that our algorithm successfully fits planes under high percentages of noise and outliers. This is superior to other state-of-the-art regularity-constrained plane fitting methods in terms of speed and robustness.

The details can be found in the following [ISPRS 2020 paper](http://www.sciencedirect.com/science/article/pii/S0924271620300095)

## Citing our work
If you find our works useful in your research, please consider citing:

	Lin Y, Li J, Wang C, Chen Z, Wang Z, and Li J. Fast regularity-constrained plane fitting. Isprs Journal of Photogrammetry & Remote Sensing, vol. 161, pages 208-217, 2020.

### BibTex
	@article{LIN2020208,
		title = "Fast regularity-constrained plane fitting",
		journal = "ISPRS Journal of Photogrammetry and Remote Sensing",
		volume = "161",
		pages = "208 - 217",
		year = "2020",
		issn = "0924-2716",
		doi = "https://doi.org/10.1016/j.isprsjprs.2020.01.009",
		url = "http://www.sciencedirect.com/science/article/pii/S0924271620300095",
		author = "Yangbin Lin and Jialian Li and Cheng Wang and Zhonggui Chen and Zongyue Wang and Jonathan Li",
	}

## Install & complie

Please directly copy the code into your workspace and compile it with any complier that supports C++11. In addition, we also use the output functions provided by the fmt library, you can refer it by the link: [fmt: A modern formatting library](https://fmt.dev/)

## Sample usage:
First, we construct a extractor, specifying the value of **k_neighbors**, **min_support_points**, **n_constraints** and **outlier_penalty** in turn. You can refer the detailed usage of the last two parameters in out paper.
```c++
cl::point_cloud::GlobalL0Extractor extractor(15, 50, 3, 1.0);
```
Then, you can extract the planes through calling on the  **ExtractPlanes** method.
```c++
extractor.ExtractPlanes(kd_tree, normals, &planes, &labels);
```
Where, 'kd_tree' is the input 3D point cloud organized by the KD-Tree. It can be read the points from XYZ file and construct the KD-Tree by calling: 
```c++
cl::geometry::io::ReadXYZPoints(filename.c_str(), &points);
KDTree<RPoint3D> kd_tree;
kd_tree.SwapPoints(&points);
```
Given the KD-tree, we can use the PCA method to calculate the normal vector for each point. For example:
```c++
Array<RVector3D> normals;
geometry::PCAEstimateNormals(kd_tree, 30, &normals);
```
It indicates we use 30 nearest neighbors of the current point to estimate its normal vector.

At last, the extracted planes and the labels for all points are stored in the 'planes' and 'labels' variables. Then we can compute the root-mean-square (RMS) of the distance from inliers to the fitting planes mentioned by our paper. The usage is as follows:

```c++
// Note! We need to ensure that the 'points' variable is not null.
kd_tree.SwapPoints(&points);
kd_tree.clear();

// Print metrics.
Metric(points, labels, planes);
```
​	
Please see main.cpp for more details.

The mentioned file 'foam_box.xyz' is also provided in the project.

## Comparison
We compare our method to the Efficient RANSAC and Pearl algorithm. You can find the implementation by the following links:

Efficient RANSAC implementation: [Efficient RANSAC for Point-Cloud Shape Detection](https://cg.cs.uni-bonn.de/en/publications/paper-details/schnabel-2007-efficient/)

PEARL implementation: [Computer Vision at Waterloo - Code](https://vision.cs.uwaterloo.ca/code/)

## Sample results. 

The first column is the orignal point cloud. The second column is Efficient RANSAC, the third is Pearl and the last is our method. The calculated normals have been mapped on the two-dimensinal space and visualized in the second row.

<img src="https://github.com/yblin/global_l0/blob/master/result.jpg" width="1000">

## Contact

Please feel free to leave suggestions or comments to Dr. Lin (yblin@jmu.edu.cn), or Prof. Wang (cwang@xmu.edu.cn)

