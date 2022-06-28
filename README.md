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

1. > containers.config

The script is located in the configs directory. Edit the following line by replacing "${HOME}/singularity" with your own path
```
cacheDir = "${HOME}/singularity"
```
This allows you to choose a location where the containers that are required will be stored.
It is best to pick a location with enough storage capacity as the total size of all the
containers could be upto 5GB

2. > pbspro.config

The script is located in the configs directory. Edit this script to suit your cluster options.
```
executor {
    name = 'pbspro'
    queue = 'normal'                                            // replace 'normal' with your cluster queue option
    queueSize = 10                                              // increase or decrease the queue size according to your privileges 
}

.
.
.

process {
    beforeScript = 'module load chpc/singularity/3.5.3'         // load singularity according to your cluster directives
    cpus = 24
    time = 3.h
    clusterOptions = '-P XXXXXXX'                              // replace 'XXXXXXX' with the project name you are assigned to on your cluster
}
``` 

The above two scripts are only edited once to setup the workflow on a new cluster.

The next script is the job configuration script which will be updated each time a different job is to be run

3. > nextflow.config 

The following tutorial will describe clearly how to use the nextflow.config file to run jobs

A. PHASING WITHOUT REFERENCE (NO IMPUTATION)
  - Make the following changes
    ```
     
    ```

