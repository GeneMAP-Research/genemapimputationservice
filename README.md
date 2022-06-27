# Phasing and Imputation
Haplotype estimation (phasing) and imputation.

Involves:

- Phasing: With reference and without reference
- Custom panel creation
- Imputation

Phasing tools implemented/to be implemented:
- EAGLE2
- SHAPEIT4
- BEAGLE5

Imputation tools implemented/to be implemented:
- Minimac4
- IMPUTE2
- BEAGLE5

Custom panel creation uses
- Minimac3

Uses ```bcftools``` throughtout and ```bgzip``` occasionally

***Requires editting config scripts found in the configs directory before running***

The following tutorial will describe clearly how to prepare the workflow and run it on a new cluster

1) Choose a location for all containers that are required: This location is called ```containers_dir``` 
and it is found in the **containersDir.config** script located in configs directory

* In my case, my ```containers_dir = '/mnt/lustre/groups/CBBI1243/KEVIN/containers/'``` 

Edit this to your preferred location. Note that the total size of all the containers needed can be greater than 10GB

2) Now pull the containers using singularity on your cluster as follows:
* Load singularity on your cluster
* Navigate to your chosen containers_dir
* Run the following commands
```

```
