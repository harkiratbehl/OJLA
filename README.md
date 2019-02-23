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


### STEP-2 Downloading and preparing the dataset

**UCF-101:**
Download the images and ssd detections for UCF-101 from [link](https://unioxfordnexus-my.sharepoint.com/:u:/g/personal/engs1635_ox_ac_uk/EYWs-k3zy1pHtOhQ3U0BDhIBCzvQc3wXvC4271dn4D0tsA?e=XVw8eh). After downloading the dataset, extract it within the repository as it is.

## TESTING
Run the file 'main.py'

	python main.py

To visualize the results, make the flag 'DISPLAY = 1' or 'DISPLAY = 2' in file 'main.py'

## TO-DO
[to-do list](https://gist.github.com/harkiratbehl/cb499c4ec19b2ab3e81c7576b3990839)

## Contact
If you have any queries at all, or need help with understanding/using the code or any part of it, please do not hesitate to contact me (Harkirat) at harkiratbehl1@gmail.com .