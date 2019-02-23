# OJLA
Official code for BMVC 2018 paper 'Incremental Tube Construction for Human Action Detection'
If you use this code, please cite this [Incremental Tube Construction for Human Action Detection](https://arxiv.org/pdf/1704.01358.pdf):


	@inproceedings{Behl2018IncrementalTC,
	  title={Incremental Tube Construction for Human Action Detection},
	  author={Harkirat S. Behl and Michael Sapienza and Gurkirt Singh and Suman Saha and Fabio Cuzzolin and Philip H. S. Torr},
	  booktitle={BMVC},
	  year={2018}
	}


## Pre Requisites

### STEP-1 Installing required libraries


The code for the tube construction algorithm is built in Matlab. Most of the functions used usually come with most Matlab distributions. Users might need to install Bioinformatics Toolbox, if it is not already present in the Matlab distribution. This is needed for the [graphconncomp](https://uk.mathworks.com/help/bioinfo/ref/graphconncomp.html) function.

Caffe and Matcaffe
putting sym link


### STEP-2 Downloading and preparing the dataset

**UCF-101**
Download the  dataset from [link]().
After downloading the dataset, extract it within the  directory.

**For your own dataset**

### STEP-3 Downloading the trained Model
The trained model for UCF-101 can be downloaded from [link]().

After downloading the frame-level detections, please put it in the folder


## TESTING
Run the file 'mian.py'

	python main.py

To visualize the results, make the flag 'DISPLAY = 1' or 'DISPLAY = 2' in file 'main.py'
