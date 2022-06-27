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


## Clone the workflow to your computer or cluster
git clone 


### Requires editting config scripts found in the configs directory before running

The following tutorial will describe clearly how to prepare the workflow and run it on a new cluster

1) Choose a location for all containers that are required
This location is called ```containers_dir``` 
and it is found in the ***containersDir.config*** script located in configs directory

  * In my case ```containers_dir = '/mnt/lustre/groups/CBBI1243/KEVIN/containers/'``` 

Edit this to your preferred location. 

> Note that the total size of all the containers needed can be upto 5GB

2) Now pull the containers using singularity on your cluster as follows
  * Load singularity on your cluster: ```module load ...```
  * Navigate to your chosen containers_dir
  * Run the following commands
  ```
  singularity pull vcftools_latest.sif docker://biocontainers/vcftools:v0.1.16-1-deb_cv1
  singularity pull bcftools_latest.sif docker://sickleinafrica/bcftools:1.11
  singularity pull beagle_latest.sif docker://sickleinafrica/beagle:latest
  singularity pull eagle2_latest.sif docker://sickleinafrica/eagle2:latest
  singularity pull plink2_latest.sif docker://sickleinafrica/plink2:latest
  singularity pull minimac4_latest.sif docker://biocontainers/minimac4:v1.0.0-2-deb_cv1
  singularity pull minimac3_latest.sif docker://rhughwhite/minimac3:2.0.1
  ```
