#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

    chromosome = getChromosomes()
    geneticMap = getGeneticMapFiles()
    thousandGenomesReference = getThousandGenomesReference()
  
    vcf = getVcf()
    vcf.join(thousandGenomesReference)
       .map{ chr, vcf, ref_file, ref_index -> tuple(chr.replaceFirst(/chr/,""), vcf, ref_file, ref_index) }
       .set{ intersect_input }

    intersects = getVcfIntersect( intersect_input )
    mergedVcfs = mergeIntersectVcf(intersects).collect().flatten()
    merged_indexed =  getVcfIndex(mergedVcfs).collect().view()
    concatenateVcfs(merged_indexed).view()
}

def getVcf() {
    return channel.fromFilePairs( params.inputDir + "*.vcf.{gz,gz.tbi}" )
                  .ifEmpty { error "\nERROR: Please check that the vcf file and its index '.tbi' exist...\n" }
                  .map{ chr, vcf -> tuple( chr.replaceFirst(/-phased/,""), vcf ) }
}

def getChromosomes() {
    return channel.of(1..22, 'X')
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.referenceDir + "chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr, ref_file.first(), ref_file.last())
	    	  }
}

def getGeneticMapFiles() {
    return channel.fromFilePairs( params.referenceDir + 'plink.chr*.GRCh37.map', size: 1 )
		  .map { group_key, map_file ->
			 tuple( group_key.replaceFirst(/^plink\.chr/,""), map_file.first() )
		  }
}

process getVcfIndex() {
    tag "BCFTOOLS INDEX: ${input_vcf}"
    label 'bcftools'
    label 'index_vcf'
    input:
        path input_vcf
    output:
        tuple \
            path("${input_vcf}"), \
            path("${input_vcf}.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads ${task.cpus} \
            ${input_vcf}
        """
}

process getVcfIntersect() {
    tag "processing chr${chrom}"
    label 'bcftools'
    label 'mediumMemory'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(vcf), \
            path(ref_vcf), \
            path(ref_index)
    output:
        publishDir path: "${params.outputDir}chr${chrom}/", mode: 'copy'
        tuple \
            val(chrom), \
            path("000{2,3}.vcf.gz")
    script:
        input_vcf = vcf[0]
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

process mergeIntersectVcf() {
    tag "Merging cohort and reference data"
    label 'bcftools'
    label 'mediumMemory'
    input:
        tuple \
            val(chrom), \
            path(intersects)
    output:
        publishDir path: "${params.outputDir}/chr${chrom}", mode: 'copy'
        path "chr${chrom}-cohort-ref.vcf.gz"
    script:
        """
        for i in $intersects; do bcftools index -ft --threads 24 \$i; done
        bcftools \
            merge \
            --force-samples \
            --threads ${task.cpus} \
            -Oz \
            -o "chr${chrom}-cohort-ref.vcf.gz" \
            ${intersects[0]} \
            ${intersects[1]}
        """

}

process concatenateVcfs() {
    tag "Concatenating cohort-ref VCFs"
    label 'bcftools'
    label 'mediumMemory'
    input:
        path vcfs
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${params.outputPrefix}-cohort-ref.vcf.gz"
    script:
        """
        bcftools \
            concat \
            -a -D \
            --threads ${task.cpus} \
            -Oz \
            -o "${params.outputPrefix}-cohort-ref.vcf.gz" \
            chr{1..22}-cohort-ref.vcf.gz \
            chrX-cohort-ref.vcf.gz
        """

}

