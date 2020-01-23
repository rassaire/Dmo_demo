# DMO-GPM
Scalable dynamic Gaussian process models (DMO-GPM)


# Description

This project aims to offer easy access to the modeling framework for multiple 3D rigid biological structures of interest from medical images (MRI, CT, etc.). It is a framework that allows for the easy creation of a dynamic multi-object models, which can be trained to detect and segment objects of interest from 3D/2D medical images. For more details please consult the paper (http://arxiv.org/abs/2001.07904)

The project code is implemented in Scalismo software(https://unibas-gravis.github.io/scalismo-microsite) and is being cleaned up for public research-related use. Accompanying data are will be provided to run the preset examples and make sure that the system is functioning on your system.
The gif below shows an illustrative summary of the main functionality, in order for the user to understand the main processing cycle. The image shows the variation of the first five principal components of the DMO-GPM of the shoulder joint (scapula and humerus), sampled from -3 to + 3 standard deviations. The first block shows the shape-pose variation, the second the marginalised pose variation, and the third the marginalised shape variation. 



<p align="center">
<img src="DMO-animation.gif" width="60%" hight="50%">
</p>

We hope this project will serve well in making the state-of-the-art dynamic shape modelling more accessible in the field of medical imaging.



