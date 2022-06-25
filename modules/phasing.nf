def getVcf() {
    return channel.fromPath( params.vcf_dir + params.vcf )
}

def getChromosomes() {
    return channel.of(1..22, 'X')
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.panel_dir + "chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr.replaceFirst(/chr/,""), ref_file.first(), ref_file.last())
	    	  }
}

def getCostumReferencePanel() {
    return channel.fromFilePairs( params.panel_dir + "wgs_chr*_phased.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
                  .map { chr, ref_file ->
                         tuple( chr.replaceFirst(/wgs_chr/,""), ref_file.first(), ref_file.last() )
                  }
}

def getEagleHapmapGeneticMap() {
    return channel.fromPath( params.panel_dir + "/tables/genetic_map_hg19_withX.txt.gz" )
}

def getShapeitGeneticMap() {
    return channel.fromPath( params.panel_dir + "/shapeit/chr*.b37.gmap.gz" )
                  .map { map -> tuple(map.simpleName.replaceAll(/chr/,""), map) }
}

def getHapmapGeneticMap() {
    return channel.fromFilePairs( params.panel_dir + "genetic_map_chr*_combined_b37.txt", size: 1 )
                  .map { chr, map_file ->
                         tuple( chr.replaceFirst(/genetic_map_chr/,""), map_file.first() )
                  }
}

def getPlinkGeneticMap() {
    return channel.fromFilePairs( params.panel_dir + 'plink.chr*.GRCh37.map', size: 1 )
                  .map { chr, map_file ->
                         tuple( chr.replaceFirst(/^plink\.chr/,""), map_file.first() )
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
        tuple \
            val(chrom), \
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

process splitVcfByChrom() {
    tag "processing chr${chrom}"
    label 'plink2'
    label 'splitVCF'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        tuple \
            val(chrom), \
            path("chr${chrom}.vcf.gz")
    script:
        """
        plink2 \
            --vcf ${input_vcf} \
            --chr ${chrom} \
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
	publishDir path: "${params.out_dir}", mode: 'copy'
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

process beaglephase() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'phaseGenotypes'
    input:
	tuple val(chrom), path(input_vcf), path(vcf_index), path(ref_vcf), path(ref_index), path(geneticMap)
    output:
	publishDir path: "${params.out_dir}", mode: 'copy'
	tuple \
            val(chrom), \
            path("chr${chrom}.*.vcf.gz")
    script:
	"""
        if [[ ${params.impute} == "true" ]]; then
           out_suffix=imputed;
        else
           out_suffix=phased
        fi
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
	    out="chr${chrom}.\${out_suffix}"
	"""
}

process eaglePhaseWithoutRef() {
    tag "processing chr${chrom}"
    label 'eagle'
    label 'phaseGenotypes'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(geneticMap)
    output:
        publishDir path: "${params.out_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path("${params.out_prefix}_chr${chrom}_phased.vcf.gz")
    script:
        """
        eagle \
          --vcf=${input_vcf} \
          --geneticMapFile=${geneticMap} \
          --chrom=${chrom} \
          --numThreads=${task.cpus} \
          --Kpbwt=${params.kpbwt} \
          --vcfOutFormat=z \
          --outPrefix=${params.out_prefix}_chr${chrom}_phased 2>&1 | \
          tee ${params.out_prefix}_chr${chrom}_phased.log
        """
}

process eaglePhaseWithRef() {
    tag "processing chr${chrom}"
    label 'eagle'
    label 'phaseGenotypes'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(ref_vcf), \
            path(ref_index), \
            path(geneticMap)
    output:
        publishDir path: "${params.out_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path("${params.out_prefix}_chr${chrom}_phased.vcf.gz")
    script:
        """
        eagle \
          --vcfTarget=${input_vcf} \
          --vcfRef=${ref_vcf} \
          --geneticMapFile=${geneticMap} \
          --chrom=${chrom} \
          --numThreads=${task.cpus} \
          --Kpbwt=${params.kpbwt} \
          --vcfOutFormat=z \
          --outPrefix=${params.out_prefix}_chr${chrom}_phased 2>&1 | \
          tee ${params.out_prefix}_chr${chrom}_phased.log
        """
}

process shapeitPhaseWithoutRef() {
    tag "processing chr${chrom}"
    label 'shapeit'
    label 'phase_shapeit'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(geneticMap)
    output:
        publishDir path: "${params.out_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path("${params.out_prefix}_chr${chrom}_phased_noref.vcf.gz")
    script:
        """
        shapeit4 \
          --input ${input_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --pbwt-depth ${params.pbwt} \
          --output ${params.out_prefix}_chr${chrom}_phased_noref 2>&1 | \
          tee ${params.out_prefix}_chr${chrom}_phased_noref.log
        """
}

process shapeitPhaseWithRef() {
    tag "processing chr${chrom}"
    label 'shapeit'
    label 'phase_shapeit'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index), \
            path(ref_vcf), \
            path(ref_index), \
            path(geneticMap)
    output:
        publishDir path: "${params.out_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path("${params.out_prefix}_chr${chrom}_phased.vcf.gz")
    script:
        """
        shapeit4 \
          --input ${input_vcf} \
          --reference ${ref_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --pbwt-depth ${params.pbwt} \
          --output ${params.out_prefix}_chr${chrom}_phased 2>&1 | \
          tee ${params.out_prefix}_chr${chrom}_phased.log
        """
}

process createLegendFile() {
    tag "processing chr${chrom}"
    label 'vcftools'
    label 'mediumMemory'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index)
    output:
        publishDir path: "${params.out_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path(input_vcf), \
            path("${input_vcf.simpleName}.legend.gz")
    script:
        """
        echo "id position a0 a1 all.aaf" > header
        vcftools \
            --gzvcf ${input_vcf} \
            --freq \
            --out ${chrom}
        sed 's/:/\\t/g' ${chrom}.frq | \
            sed 1d | \
            awk '{print \$1":"\$2" "\$2" "\$5" "\$7" "\$8}' \
            > ${input_vcf.simpleName}.legend
        cat header ${input_vcf.simpleName}.legend | \
            bgzip \
            > ${input_vcf.simpleName}.legend.gz
        """
}

process prepareChrXPanel() {
    tag "processing chr${chrom}"
    label 'eagle'
    label 'phaseGenotypes'
    input:
        tuple \
            val(chrom), \
            path(input_vcf)
    output:
        publishDir path: "${params.out_dir}", mode: 'copy'
        tuple \
            val(chrom), \
            path("${input_vcf.simpleName}_PAR1.vcf.gz"), \
            path("${input_vcf.simpleName}_PAR1,nonPAR.vcf.gz.tbi"), \
            path("${input_vcf.simpleName}_PAR2.vcf.gz"), \
            path("${input_vcf.simpleName}_PAR2,nonPAR.vcf.gz.tbi"), \
            path("${input_vcf.simpleName}_nonPAR.vcf.gz"), \
            path("${input_vcf.simpleName}_nonPAR.vcf.gz.tbi")
    script:
        """
        bcftools \
            index \
            -ft \
            --threads \
            ${task.cpus} \
            ${input_vcf}
        bcftools \
            view \
            -r X:60001-2699520 \
            ${input_vcf} \
            -Oz | \
            tee ${input_vcf.simpleName}_PAR1.vcf.gz | \
            bcftools \
                index \
                -ft \
                --threads ${task.cpus} \
                --output ${input_vcf.simpleName}_PAR1.vcf.gz.tbi
        bcftools \
            view \
            -r X:154931044-155260560 \
            ${input_vcf} \
            -Oz | \
            tee ${input_vcf.simpleName}_PAR2.vcf.gz | \
            bcftools \
                index \
                -ft \
                --threads ${task.cpus} \
                --output ${input_vcf.simpleName}_PAR2.vcf.gz.tbi
        bcftools \
            view \
            -r X:2699521-154931043 \
            ${input_vcf} \
            -Oz | \
            tee ${input_vcf.simpleName}_nonPAR.vcf.gz | \
            bcftools \
                index \
                -ft \
                --threads ${task.cpus} \
                --output ${input_vcf.simpleName}_nonPAR.vcf.gz.tbi            
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

