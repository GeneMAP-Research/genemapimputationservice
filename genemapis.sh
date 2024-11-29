#!/usr/bin/env sh

#--- genemapimputation workflow wrapper ---#

function usage() {
   echo """
   =================================================================
   GeneMAP-NGS ~ a wrapper for the nextflow-based genemapis workflow
   =================================================================

   Usage: genemapis <workflow> <profile> [options] ...

           workflows:
           ---------
               test: Run test to see if workfow installed correctly.
              align: Check strand and reference allele overlap.
	      phase: Phase (and Impute) genotypes.
             impute: Only impute missing genotypes.
           xxxxxxxx: .
           xxxxxxxx: .
          xxxxxxxxx: .


           profiles: <executor>,<container>,<reference>
           ---------
          executors: local, slurm
         containers: singularity, aptainer, docker
          reference: hg19, hg38, t2t


           examples:
           ---------
         genemapis test slurm,singularity,hg38 [options]
         genemapis align local,singularity,hg19 --vcf my-vcf.vcf.gz 
   """
}


#####################################################################################################################
function checkprofile() {
   #profile was passed as only argument 
   #so it takes position $1 here
   if [[ "$1" == "" ]]; then
      echo "ERROR: please specify a profile to use!";
      usage;
      exit 1;
   elif [[ $1 == -* ]]; then
      echo "ERROR: please specify a valid profile!";
      usage;
      exit 1;
   else
      local profile="$1"
   fi
}

function check_required_params() {
   for params_vals in $@; do
      #get each param and its value as an array
      param_val=( $(echo ${params_vals} | sed 's/,/ /g') )
      
      #slice the array to its consituent params and values
      param=${param_val[0]}
      val=${param_val[1]}

      #now check each param and its value
      if [[ $val == -* ]] || [[ $val == NULL ]]; then
         echo "ERROR: Invalid paramter value for option '--${param}'";
         break;
         exit 1;
      fi
   done
}

function check_optional_params() {
   for params_vals in $@; do
      #get each param and its value as an array
      param_val=( $(echo ${params_vals} | sed 's/,/ /g') )

      #slice the array to its consituent params and values
      param=${param_val[0]}
      val=${param_val[1]}

      #now check each param and its value
      if [[ $val == -* ]]; then
         echo "ERROR: Invalid paramter value for option '--${param}'";
         break;
         exit 1;
      fi
   done
}

function setglobalparams() {
#- create the project nextflow config file
echo """includeConfig \"\${projectDir}/nextflow.config\"
includeConfig \"\${projectDir}/configs/profile-selector.config\"
includeConfig \"\${projectDir}/configs/resource-selector.config\"
"""
}

##################################################### USAGE #########################################################
function alignusage() {
   echo -e "\nUsage: genemapis align <profile> [options] ..."
   echo """
           options:
           --------
           --autosome           : (optional) specify this flag to process autosomes only.
           --vcf                : (required) VCF file. Must specify full path.
           --output_prefix      : (optional) output prefix [default: myout].
           --output_dir         : (required) path to save output files.
           --exclude_sample     : (optional) reference samples to exclude. Single column, one sample per line.
           --threads            : (optional) number of computer cpus to use [default: 8].
           --njobs              : (optional) number of jobs to submit at once [default: 4]
           --help               : print this help message.

   """
}

function phaseusage() {
   echo -e "\nUsage: genemapis phase <profile> [options] ..."
   echo """
           options:
           --------
	   --with_ref           : (optional) specify this flag to phase with reference
	   --phase_tool         : (optional) shapeit4, eagle2, beagle5 [default: shapeit4]
	   --impute             : (optional) to phase and impute in one run, add this flag
	   --impute_tool        : (optional) minimac4, impute2, beagle5 [default: minimac4] 
           --vcf                : (required) VCF file. Must specify full path
	   --autosome           : (optional) specify this flag to process autosomes only
           --output_prefix      : (optional) output prefix [default: myout]
           --output_dir         : (required) path to save output files
           --burnit             : (optional) number of BEAGLE burn in iterations [default: 2000]
           --mainit             : (optional) number of BEAGLE main iterations [default: 1000]
           --kpbwt              : (optional) number of BEAGLE and EAGLE conditioning haplotypes [default: 50000]
	   --pbwt               : (optional) number of SHAPEIT PBWT iterations [default: 8]
           --threads            : (optional) number of computer cpus to use [default: 8]
           --njobs              : (optional) number of jobs to submit at once [default: 4]
           --help               : print this help message

   """
}

function imputeusage() {
   echo -e "\nUsage: genemapis impute <profile> [options] ..."
   echo """
           options:
           --------
           --impute_tool        : (optional) minimac4, impute2, beagle5 [default: minimac4] 
	   --panel              : (optional) available panels: custom, kgp, h3a [default: kgp]
           --vcf                : (required) VCF file. Must specify full path
           --autosome           : (optional) specify this flag to process autosomes only
           --output_prefix      : (optional) output prefix [default: myout]
           --output_dir         : (required) path to save output files
           --burnit             : (optional) number of BEAGLE burn in iterations [default: 2000]
           --mainit             : (optional) number of BEAGLE main iterations [default: 1000]
           --kpbwt              : (optional) number of BEAGLE and EAGLE conditioning haplotypes [default: 50000]
           --pbwt               : (optional) number of SHAPEIT PBWT iterations [default: 8]
           --threads            : (optional) number of computer cpus to use [default: 8]
           --njobs              : (optional) number of jobs to submit at once [default: 4]
           --help               : print this help message

   """
}


############################################# CONFIGURATION FILES ####################################################
function testconfig() {
#check and remove test config file if it exists
[ -e test.config ] && rm test.config

# $indir $bpm $csv $cluster $fasta $bam $out $outdir $thrds
echo """includeConfig \"\${projectDir}/nextflow.config\"
includeConfig \"\${projectDir}/configs/profile-selector.config\"
includeConfig \"\${projectDir}/configs/test.config\"
""" >> test.config
}


function alignconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${3}-align.config ] && rm ${3}-align.config

#alignconfig $autosome $vcf $output_prefix $output_dir $exclude_sample $threads $njobs
echo """`setglobalparams`

params {
  //====================================
  // genemapis align workflow parameters 
  //====================================
  autosome        = $1                          // (optional) include this argument if you wish to only process autosomes [default: false]
  vcf             = '$2'                        // (required) VCF file. Must specify full path.
  output_prefix   = '$3'                        // (optional) output prefix [default: myout]
  output_dir      = '$4'                        // (required) path to save output files.
  exclude_sample  = '$5'                        // (optional) reference samples to exclude. Single column, one sample per line.
  threads         = $6                          // number of computer cpus to use [default: 8]
  njobs           = $7                          // (optional) number of jobs to submit at once [default: 4]
}
""" >> ${3}-align.config
}


function phaseconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${7}-phase.config ] && rm ${7}-phase.config

echo """
params {
  //====================================
  // genemapis phase workflow parameters
  //====================================

  phase           = true
  with_ref        = ${1} 
  phase_tool      = '${2}'
  impute          = ${3}
  impute_tool     = '${4}'
  vcf             = '${5}'
  autosome        = ${6}
  output_prefix   = '${7}'
  output_dir      = '${8}'
  burnit          = ${9}
  mainit          = ${10}
  kpbwt           = ${11}
  pbwt            = ${12}
  threads         = ${13}
  njobs           = ${14}

  /***************************************************************************************
  ~ with_ref: (optional) specify this flag to phase with reference
    available ref: kgp
  ~ phase_tool: (optional) shapeit4, eagle2, beagle5 [default: shapeit4]
  ~ impute: (optional) specify this flag to impute
  ~ impute_tool: (optional) minimac4, impute2, beagle5 [default: minimac4]
  ~ panel: (optional) available imputation panels: custom, kgp, h3a [default: kgp]
  ~ vcf: (required) VCF file. Must specify full path
  ~ autosome: (optional) specify this flag to process autosomes only
  ~ output_prefix: (optional) output prefix [default: myout]
  ~ output_dir: (required) path to save output files
  ~ burnit: (optional) number of BEAGLE burn in iterations [default: 2000]
  ~ mainit: (optional) number of BEAGLE main iterations [default: 1000]
  ~ kpbwt: (optional) number of BEAGLE and EAGLE conditioning haplotypes [default: 50000]
  ~ pbwt: (optional) number of SHAPEIT PBWT iterations [default: 8]
  ~ threads: (optional) number of computer cpus to use [default: 8]
  ~ njobs: number of simultaneous jobs to submit [default: 4]
  ****************************************************************************************/
}

`setglobalparams`
""" >> ${7}-phase.config
}

function imputeconfig() { #params passed as arguments
#check and remove config file if it exists
[ -e ${5}-impute.config ] && rm ${5}-impute.config

echo """
params {
  //=====================================
  // genemapis impute workflow parameters
  //=====================================

  impute          = true
  phase           = false
  impute_tool     = '${1}'
  panel           = '${2}'
  vcf             = '${3}'
  autosome        = ${4}
  output_prefix   = '${5}'
  output_dir      = '${6}'
  burnit          = ${7}
  mainit          = ${8}
  kpbwt           = ${9}
  pbwt            = ${10}
  threads         = ${11}
  njobs           = ${12}

  /***************************************************************************************
  ~ impute_tool: (optional) minimac4, impute2, beagle5 [default: minimac4]
  ~ panel: (optional) available imputation panels: custom, kgp, h3a [default: kgp]
  ~ vcf: (required) VCF file. Must specify full path
  ~ autosome: (optional) specify this flag to process autosomes only
  ~ output_prefix: (optional) output prefix [default: myout]
  ~ output_dir: (required) path to save output files
  ~ burnit: (optional) number of BEAGLE burn in iterations [default: 2000]
  ~ mainit: (optional) number of BEAGLE main iterations [default: 1000]
  ~ kpbwt: (optional) number of BEAGLE and EAGLE conditioning haplotypes [default: 50000]
  ~ pbwt: (optional) number of SHAPEIT PBWT iterations [default: 8]
  ~ threads: (optional) number of computer cpus to use [default: 8]
  ****************************************************************************************/
}

`setglobalparams`
""" >> ${5}-impute.config
}


if [ $# -lt 1 ]; then
   usage; exit 1;
else
   case $1 in
      test)
         profile='local,singularity,hg19'
         testconfig
      ;;
      align)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            alignusage;
            exit 1;
         fi

         prog=`getopt -a --long "help,autosome,vcf:,output_prefix:,output_dir:,exclude_sample:,threads:,njobs:" -n "${0##*/}" -- "$@"`;
         
         # defaults
         autosome=false        #// (optional) include this argument if you wish to only process autosomes [default: false]
         vcf=NULL             #// (required) VCF file. Must specify full path.
         output_prefix=myout  #// (optional) output prefix [default: myout]
         output_dir=NULL      #// (required) path to save output files.
         exclude_sample=NULL  #// (optional) reference samples to exclude. Single column, one sample per line.
         threads=8            #// number of computer cpus to use [default: 8]
         njobs=4              #// (optional) number of jobs to submit at once [default: 4]
         
         eval set -- "$prog"

         while true; do
            case $1 in
               --autosome) autosome=true; shift;;
               --vcf) vcf="$2"; shift 2;;
               --output_prefix) output_prefix="$2"; shift 2;;
               --output_dir) output_dir="$2"; shift 2;;
               --exclude_sample) exclude_sample="$2"; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; alignusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; alignusage; 1>&2; exit 1;;
            esac
         done

         #- check required options
         check_required_params vcf,$vcf
         check_output_dir $output_dir
         check_optional_params output_prefix,$output_prefix exclude_sample,$exclude_sample threads,$threads njobs,$njobs

         alignconfig $autosome $vcf $output_prefix $output_dir $exclude_sample $threads $njobs
         #echo `nextflow -c ${out}-qc.config run qualitycontrol.nf -profile $profile -w ${outdir}/work/`
      ;;
      phase)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            phaseusage;
            exit 1;
         fi

         prog=`getopt -a --long "with_ref,phase_tool:,impute,impute_tool:,vcf:,autosome,output_prefix:,output_dir:,burnit:,mainit:,kpbwt:,pbwt:,threads:,njobs:" -n "${0##*/}" -- "$@"`; 


         #- defaults

         with_ref=false
         phase_tool=shapeit4
         impute=false
         impute_tool=minimac4
         vcf=NULL
         autosome=false
         output_prefix=my-phase
         output_dir=NULL
         burnit=2000
         mainit=1000
         kpbwt=50000
         pbwt=8
         threads=8
	 
         eval set -- "$prog"

         while true; do
            case $1 in
	       --with_ref) with_ref=true; shift;;
	       --phase_tool) phase_tool=$2; shift 2;;
	       --impute) impute=true; shift;;
	       --impute_tool) impute_tool=$2; shift 2;;
	       --vcf) vcf=$2; shift 2;;
	       --autosome) autosome=true; shift;;
	       --output_prefix) output_prefix=$2; shift 2;;
	       --output_dir) output_dir=$2; shift 2;;
	       --burnit) burnit=$2; shift 2;;
	       --mainit) mainit=$2; shift 2;;
	       --kpbwt) kpbwt=$2; shift 2;;
	       --pbwt) pbwt=$2; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; phaseusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; phaseusage; 1>&2; exit 1;;
            esac
            continue; shift;
         done

         #- check required options
         #check_ftype $ftype
         check_required_params \
	    vcf,$vcf \
	    output_dir,$output_dir && \
         check_optional_params \
	    phase_tool,$phase_tool \
	    impute_tool,$impute_tool \
	    output_prefix,$output_prefix \
	    burnit,$burnit \
	    mainit,$mainit \
	    kpbwt,$kpbwt \
	    pbwt,$pbwt \
	    threads,$threads \
	    njpbs,$njobs && \
         phaseconfig \
	    $with_ref \
	    $phase_tool \
	    $impute \
	    $impute_tool \
	    $vcf \
	    $autosome \
	    $output_prefix \
	    $output_dir \
	    $burnit \
	    $mainit \
	    $kpbwt \
	    $pbwt \
	    $threads \
	    $njobs
      ;;
      impute)
         #pass profile as argument
         checkprofile $2;
         profile=$2;
         shift;
         if [ $# -lt 2 ]; then
            imputeusage;
            exit 1;
         fi

         prog=`getopt -a --long "impute_tool:,panel:,vcf:,autosome,output_prefix:,output_dir:,burnit:,mainit:,kpbwt:,pbwt:,threads:,njobs:" -n "${0##*/}" -- "$@"`; 


         #- defaults

         impute_tool=minimac4
	 panel=kgp
         vcf=NULL
         autosome=false
         output_prefix=my-impute
         output_dir=NULL
         burnit=2000
         mainit=1000
         kpbwt=50000
         pbwt=8
         threads=8
	 
         eval set -- "$prog"

         while true; do
            case $1 in
	       --impute_tool) impute_tool=$2; shift 2;;
	       --panel) panel=$2; shift 2;;
	       --vcf) vcf=$2; shift 2;;
	       --autosome) autosome=true; shift;;
	       --output_prefix) output_prefix=$2; shift 2;;
	       --output_dir) output_dir=$2; shift 2;;
	       --burnit) burnit=$2; shift 2;;
	       --mainit) mainit=$2; shift 2;;
	       --kpbwt) kpbwt=$2; shift 2;;
	       --pbwt) pbwt=$2; shift 2;;
               --threads) threads="$2"; shift 2;;
               --njobs) njobs="$2"; shift 2;;
               --help) shift; imputeusage; 1>&2; exit 1;;
               --) shift; break;;
               *) shift; imputeusage; 1>&2; exit 1;;
            esac
            continue; shift;
         done

         #- check required options
         #check_ftype $ftype
         check_required_params \
	    vcf,$vcf \
	    output_dir,$output_dir && \
         check_optional_params \
	    impute_tool,$impute_tool \
	    panel,$panel \
	    output_prefix,$output_prefix \
	    burnit,$burnit \
	    mainit,$mainit \
	    kpbwt,$kpbwt \
	    pbwt,$pbwt \
	    threads,$threads \
	    njpbs,$njobs && \
         imputeconfig \
	    $impute_tool \
	    $panel \
	    $vcf \
	    $autosome \
	    $output_prefix \
	    $output_dir \
	    $burnit \
	    $mainit \
	    $kpbwt \
	    $pbwt \
	    $threads \
	    $njobs
      ;;      
      *) shift; usage; exit 1;;
   esac
fi


