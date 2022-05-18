#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

    chromosome = getChromosomes()
    geneticMap = getGeneticMapFiles()
    //[chrom, thousandGenomesReference, kgp_index] = getThousandGenomesReference()
    thousandGenomesReference = getThousandGenomesReference().view()


    vcf = getVcf()

    chromosome
        .combine(vcf)
        .set { splitvcf_input }
    
    perChromosomeVcfFiles = splitVcfByChrom( splitvcf_input )   
    perChromosomeVcfFiles
	.map { chrom, vcfFile -> tuple( "${chrom}", vcfFile ) }
	.join( thousandGenomesReference )
	.join( geneticMap )
	.set { checkstrand_input }

    alignedVcfs = alignGenotypesToReference( checkstrand_input )
    alignedVcfs
        .map { chrom, vcfFile, logFile -> tuple( "${chrom}", vcfFile, logFile ) }
        .join( thousandGenomesReference )
        .join( geneticMap )
        .set { phasing_input }

/*
    phasedVcfFiles = phaseGenotypes( phasing_input )
    getVcfIndex(phasedVcfFiles).view()

    prePhasingQualityReports = getCheckStrandReports()
*/
	
}

def getVcf() {
    return channel.fromPath( params.inputDir + params.vcf )
}

def getChromosomes() {
    return channel.of(1..22, 'X')
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.referenceDir + "chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr.replaceFirst(/chr/,""), ref_file.first(), ref_file.last())
	    	  }
}

/*
def getThousandGenomesReference() {
    return channel.fromFilePairs( params.referenceDir + 'chr*1kg.phase3.v5a.vcf.gz', size: 1 )
	   	  .map { group_key, ref_file -> 
			 tuple( group_key.replaceFirst(/chr/,""), ref_file.first())
	    	  }
}
*/

def getGeneticMapFiles() {
    return channel.fromFilePairs( params.referenceDir + 'plink.chr*.GRCh37.map', size: 1 )
		  .map { group_key, map_file ->
			 tuple( group_key.replaceFirst(/^plink\.chr/,""), map_file.first() )
		  }
}

/*
process getPalindromicSnvs() {
    tag 'FILTER PALINDROMES'
    label 'Rbase'
    label 'extractPalindromes'
    input:
        tuple val(plink_base), path(plink_fileset)
    output:
        //publishDir path: "${params.outputDir}"
        path "${params.cohortName}.palindromic.snvs.txt"
    script:
        template "filterPalindromicVariants.r"
}
*/

process getVcfIndex() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        path "*.tbi"
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process splitVcfByChrom() {
    tag "processing chr${chrom}"
    label 'plink2'
    label 'splitVCF'
    cache 'lenient'
    input:
        tuple val(chrom), path(input_vcf)
    output:
        publishDir path: "${params.outputDir}"
        tuple val(chrom), path("chr${chrom}.vcf.gz")
    script:
        """
        plink2 \
            --vcf ${input_vcf} \
            --chr ${chrom} \
            --maf 0.01 \
            --max-alleles 2 \
            --vcf-half-call missing \
            --export vcf-4.2 bgz id-paste='iid' \
            --threads ${task.cpus} \
            --out "chr${chrom}"
        """
}

process alignGenotypesToReference() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'alignGenotypes'
    input:
	tuple val(chrom), path(input_vcf), path(ref_vcf), path(ref_index), path(geneticMap)
    output:
	publishDir path: "${params.outputDir}", mode: 'copy'
        tuple val("${chrom}"), path("chr${chrom}-aligned.vcf.gz"), path("chr${chrom}-aligned.log")
    script:
	"""
	conform-gt \
	    gt="${input_vcf}" \
	    ref="${ref_vcf}" \
	    chrom=${chrom} \
	    match=POS \
            excludesamples="${params.RefereceSamplesToExclude}" \
	    out="chr${chrom}-aligned"
	"""
}

process phaseGenotypes() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'phaseGenotypes'
    input:
	tuple val(chrom), path(input_vcf), path(vcf_log), path(ref_vcf), path(ref_index), path(geneticMap)
    output:
	publishDir path: "${params.outputDir}", mode: 'copy'
	tuple \
            val(chrom), \
            path("chr${chrom}-phased.vcf.gz")
    script:
	"""
	beagle \
	    gt="${input_vcf}" \
	    ref="${ref_vcf}" \
	    map="${geneticMap}" \
	    chrom=${chrom} \
	    burnin="${params.burninIter}" \
	    iterations="${params.mainIter}" \
	    ne=20000 \
	    impute=${params.impute} \
	    nthreads=${task.cpus} \
	    out="chr${chrom}-phased"
	"""
}

process getVcfIntersect() {
    tag "processing chr${chrom}"
    label 'bcftools'
    label 'mediumMemory'
    cache 'lenient'
    input:
        tuple val(chrom), path(input_vcf), path(vcf_index), path(ref_vcf), path(ref_index)
    output:
        publishDir path: "${params.outputDir}chr${chrom}/", mode: 'copy'
        path "000{2,3}.vcf.gz"
    script:
        """
        mkdir -p "${params.outputDir}/chr${chrom}"
        
        bcftools \
            isec \
            -r ${chrom} \
            --threads ${task.cpus} \
            -p . \
            -Oz \
            -o "chr${chrom}" \
            ${input_vcf} \
            ${ref_vcf}
        """
}

process concatenateIntersectVcf() {
    tag "CONCAT"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path intersects
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path ""
    script:
        """
        bcftools \
            concat \
            --threads ${task.cpus} \
            -Oz \
            -o "${params.outputPrefix}-kgp-intersect.vcf.ga" \
            ${intersects}
        """

}

