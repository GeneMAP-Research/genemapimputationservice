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

3. > reference-hg*.config

There are two config files for reference sequences/panels: reference-hg19.config and reference-hg38.config

The files are located in the configs directory. Edit them accordingly to point the workflow to the appropriate reference sequence/panels location

For instance, for reference-hg19.config
```
params {
    panel_dir = ''			// directory where imputation reference panels are stored or will be stored
    ref_dir = ''			// directory containing fasta reference sequence
    fastaRef = ''			// name of the fasta reference sequence
    buildVersion = 'hg19'		// DO NOT EDIT THIS LINE!!!
}

```

The panel directory should contain subdirectories with different reference data types. See the example file tree below
```
referencePanels
├── hapmap
│   ├── genetic_map_chr10_combined_b37.txt
│   ├── genetic_map_chr11_combined_b37.txt
│   ├── genetic_map_chr12_combined_b37.txt
│   ├── genetic_map_chr13_combined_b37.txt
│   ├── genetic_map_chr14_combined_b37.txt
│   ├── genetic_map_chr15_combined_b37.txt
│   ├── genetic_map_chr16_combined_b37.txt
│   ├── genetic_map_chr17_combined_b37.txt
│   ├── genetic_map_chr18_combined_b37.txt
│   ├── genetic_map_chr19_combined_b37.txt
│   ├── genetic_map_chr1_combined_b37.txt
│   ├── genetic_map_chr20_combined_b37.txt
│   ├── genetic_map_chr21_combined_b37.txt
│   ├── genetic_map_chr22_combined_b37.txt
│   ├── genetic_map_chr23_combined_b37.txt
│   ├── genetic_map_chr24_combined_b37.txt
│   ├── genetic_map_chr2_combined_b37.txt
│   ├── genetic_map_chr3_combined_b37.txt
│   ├── genetic_map_chr4_combined_b37.txt
│   ├── genetic_map_chr5_combined_b37.txt
│   ├── genetic_map_chr6_combined_b37.txt
│   ├── genetic_map_chr7_combined_b37.txt
│   ├── genetic_map_chr8_combined_b37.txt
│   ├── genetic_map_chr9_combined_b37.txt
│   ├── genetic_map_chrX_combined_b37.txt
│   ├── genetic_map_chrX_nonPAR_combined_b37.txt
│   ├── genetic_map_chrX_PAR1_combined_b37.txt
│   └── genetic_map_chrX_PAR2_combined_b37.txt
├── kgp
│   ├── chr10.1kg.phase3.v5a.vcf.gz
│   ├── chr10.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr11.1kg.phase3.v5a.vcf.gz
│   ├── chr11.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr1.1kg.phase3.v5a.vcf.gz
│   ├── chr1.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr12.1kg.phase3.v5a.vcf.gz
│   ├── chr12.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr13.1kg.phase3.v5a.vcf.gz
│   ├── chr13.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr14.1kg.phase3.v5a.vcf.gz
│   ├── chr14.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr15.1kg.phase3.v5a.vcf.gz
│   ├── chr15.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr16.1kg.phase3.v5a.vcf.gz
│   ├── chr16.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr17.1kg.phase3.v5a.vcf.gz
│   ├── chr17.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr18.1kg.phase3.v5a.vcf.gz
│   ├── chr18.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr19.1kg.phase3.v5a.vcf.gz
│   ├── chr19.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr20.1kg.phase3.v5a.vcf.gz
│   ├── chr20.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr21.1kg.phase3.v5a.vcf.gz
│   ├── chr21.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr2.1kg.phase3.v5a.vcf.gz
│   ├── chr2.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr22.1kg.phase3.v5a.vcf.gz
│   ├── chr22.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr3.1kg.phase3.v5a.vcf.gz
│   ├── chr3.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr4.1kg.phase3.v5a.vcf.gz
│   ├── chr4.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr5.1kg.phase3.v5a.vcf.gz
│   ├── chr5.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr6.1kg.phase3.v5a.vcf.gz
│   ├── chr6.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr7.1kg.phase3.v5a.vcf.gz
│   ├── chr7.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr8.1kg.phase3.v5a.vcf.gz
│   ├── chr8.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chr9.1kg.phase3.v5a.vcf.gz
│   ├── chr9.1kg.phase3.v5a.vcf.gz.tbi
│   ├── chrX.1kg.phase3.v5a.vcf.gz
│   └── chrX.1kg.phase3.v5a.vcf.gz.tbi
├── plinkmap
│   ├── plink.chr10.GRCh37.map
│   ├── plink.chr11.GRCh37.map
│   ├── plink.chr12.GRCh37.map
│   ├── plink.chr13.GRCh37.map
│   ├── plink.chr14.GRCh37.map
│   ├── plink.chr15.GRCh37.map
│   ├── plink.chr16.GRCh37.map
│   ├── plink.chr17.GRCh37.map
│   ├── plink.chr18.GRCh37.map
│   ├── plink.chr19.GRCh37.map
│   ├── plink.chr1.GRCh37.map
│   ├── plink.chr20.GRCh37.map
│   ├── plink.chr21.GRCh37.map
│   ├── plink.chr22.GRCh37.map
│   ├── plink.chr2.GRCh37.map
│   ├── plink.chr3.GRCh37.map
│   ├── plink.chr4.GRCh37.map
│   ├── plink.chr5.GRCh37.map
│   ├── plink.chr6.GRCh37.map
│   ├── plink.chr7.GRCh37.map
│   ├── plink.chr8.GRCh37.map
│   ├── plink.chr9.GRCh37.map
│   └── plink.chrX.GRCh37.map
├── shapeit
│   ├── chr10.b37.gmap.gz
│   ├── chr10.b38.gmap.gz
│   ├── chr11.b37.gmap.gz
│   ├── chr11.b38.gmap.gz
│   ├── chr12.b37.gmap.gz
│   ├── chr12.b38.gmap.gz
│   ├── chr13.b37.gmap.gz
│   ├── chr13.b38.gmap.gz
│   ├── chr14.b37.gmap.gz
│   ├── chr14.b38.gmap.gz
│   ├── chr15.b37.gmap.gz
│   ├── chr15.b38.gmap.gz
│   ├── chr16.b37.gmap.gz
│   ├── chr16.b38.gmap.gz
│   ├── chr17.b37.gmap.gz
│   ├── chr17.b38.gmap.gz
│   ├── chr18.b37.gmap.gz
│   ├── chr18.b38.gmap.gz
│   ├── chr19.b37.gmap.gz
│   ├── chr19.b38.gmap.gz
│   ├── chr1.b37.gmap.gz
│   ├── chr1.b38.gmap.gz
│   ├── chr20.b37.gmap.gz
│   ├── chr20.b38.gmap.gz
│   ├── chr21.b37.gmap.gz
│   ├── chr21.b38.gmap.gz
│   ├── chr22.b37.gmap.gz
│   ├── chr22.b38.gmap.gz
│   ├── chr2.b37.gmap.gz
│   ├── chr2.b38.gmap.gz
│   ├── chr3.b37.gmap.gz
│   ├── chr3.b38.gmap.gz
│   ├── chr4.b37.gmap.gz
│   ├── chr4.b38.gmap.gz
│   ├── chr5.b37.gmap.gz
│   ├── chr5.b38.gmap.gz
│   ├── chr6.b37.gmap.gz
│   ├── chr6.b38.gmap.gz
│   ├── chr7.b37.gmap.gz
│   ├── chr7.b38.gmap.gz
│   ├── chr8.b37.gmap.gz
│   ├── chr8.b38.gmap.gz
│   ├── chr9.b37.gmap.gz
│   ├── chr9.b38.gmap.gz
│   ├── chrX.b37.gmap.gz
│   ├── chrX.b38.gmap.gz
│   ├── chrX_par1.b37.gmap.gz
│   ├── chrX_par1.b38.gmap.gz
│   ├── chrX_par2.b37.gmap.gz
│   └── chrX_par2.b38.gmap.gz
└── tables
    ├── genetic_map_1cMperMb.txt
    ├── genetic_map_hg17_withX.txt.gz
    ├── genetic_map_hg18_withX.txt.gz
    ├── genetic_map_hg19_withX.txt.gz
    └── genetic_map_hg38_withX.txt.gz
```


The above three scripts are only edited once to set up the workflow on a new cluster.

The following script is the job configuration script which will be updated each time a different job is to be run

4. > nextflow.config 

Edit the following according to you input paramters.
```
params {
    phase = true
    with_ref = false
    impute = false
    phase_tool = 'shapeit4'			// options: shapeit4, eagle2, beagle5. [default: shapeit4] 
    impute_tool = 'minimac4'			// options: minimac4, impute2, beagle5. [default: minimac4]

    vcf = ''		    			// insert the VCF file to process here, with the full (absolute) path.
    out_dir = ''			    	// full path where the results will be saved
    out_prefix = ''				// base name of all output files

    RefereceSamplesToExclude = ''

    burninIter = 2000                   	// BEAGLE5 (this is only applicable to BEALGE5)
    mainIter = 1000                     	// BEAGLE5 (this is only applicable to BEALGE5)
    kpbwt = 20000                       	// BEAGLE5 | EAGLE2 (this is applicable to both EAGLE2 and BEAGLE5)
    pbwt = 8                            	// SHAPEIT4 (this is only applicable to SHAPEIT4)
    autosome = false

    // CREATING IMPUTATION PANEL
    input_dir = ''                      	// directory containing phased VCF files (single chromosomes) for which imputation panel is to be created
}
```

#### The following tutorial will describe how to use the nextflow.config file to run jobs

--------------------
A. PHASING WITHOUT REFERENCE (NO IMPUTATION)
The parameter scope should look like this
```
phase = true
with_ref = false
impute = false
phase_tool = 'shapeit4'
impute_tool = 'minimac4'
```
Then run the workflow as follows
```
nextflow run phaseAndImputeGenotypes.nf -w /path/to/work_directory/ -profile pbspro,hg19
```
- -w: this is the path where all temporary files will be stored. It should prefarably be a location with enough storage capacity
- -profile: the profile selector is added to make it easy to run the workflow on any cluster/system and using reference panels in different genome builds. 
Therefore, if your data set is in GRCh37 coordinate (hg19 or b37), then select hg19, etc.

------------------
B. IMPUTE GENOTYPES ONLY

This requires pre-phased data in single chromosome VCF files, bgzipped and in a single directory. The directory must also contain the tabix indexed for each VCF file.

The directory of the pre-phased data is supplied using ```input_dir``` in the nextflow.config file.

Example
```
chr1.vcf.gz
chr1.vcf.gz.tbi
.
.
.
chrX.vcf.gz
chrX.vcf.gz.tbi
```
NB: The directory must only contain the VCF files to be process (and their indexes) and must not contain any other files ending in '.vcf.gz'

---------------------
C. PHASING WITHOUT REFERENCE AND THEN IMPUTE GENOTYPES

Make the following change in the nextflow.config script ```impute = true``` then run the workflow
```
nextflow run phaseAndImputeGenotypes.nf -w /path/to/work_directory/ -profile pbspro,hg19
```


