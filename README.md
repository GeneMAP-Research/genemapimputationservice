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

Uses ```bcftools``` throughtout and ```bgzip``` and ```plink2``` occasionally.

# Installing and using the workflow
Clone the workflow to your computer or cluster.
```
$ git clone https://github.com/esohkevin/siaImputationService.git
```
---------------------

There are three main configuration scripts that need to be editted before running the workflow

1) **containers.config**: Edit the following line by replacing "${HOME}/singularity" with your own path
```
cacheDir = "${HOME}/singularity"
```
This allows you to choose a location where the containers that are required will be stored.

It is best to pick a location with enough storage capacity as the total size of all the
containers could be upto 5GB

2) > pbspro.config
The following tutorial will describe clearly how to prepare the workflow and run it on a new cluster

