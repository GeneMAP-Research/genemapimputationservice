def getVcf() {
    return channel.fromPath( params.vcf )
}

def listChromosomes() {
    return channel.of(1..22, 'X')
}

def getChromosomes() {
    if(params.autosome == true) {
       channel.of(1..22)
    } else {
       channel.of(1..22, 'X')
    }
}

/*
def getVcfFileset() {
    getVcf().set { vcf }
    getChromosomes().set { chromosome }
    chromosome.combine(vcf).set { split_vcf_input }
    splitVcfByChrom(split_vcf_input).set { per_chr_vcf }
    getVcfIndex(per_chr_vcf).set { vcf_fileset }
    vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }
    return vcfFileset
}
*/

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.panel_dir + "/kgp/chr*1kg.phase3.v5a.vcf.{gz,gz.tbi}", size: 2 )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the reference file and its index '.tbi' exist...\n" }
	   	  .map { chr, ref_file -> 
			 tuple( chr.replaceFirst(/chr/,""), ref_file.first(), ref_file.last())
	    	  }
}

def getMinimacReference() {
    return channel.fromFilePairs( params.panel_dir + "/m3vcfs/*.{m3vcf.gz,rec,erate}", size: 3 )
                  .ifEmpty { error: println "\nAn error occurred! Please check that the reference files exist...\n" }
                  .map { chr, ref_fileset ->
                         tuple( chr.replaceFirst(/chr/,""), ref_fileset[0], ref_fileset[1], ref_fileset[2])
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
                  .ifEmpty { error: println "\nAn error occurred! Please check that the map files exist...\n" }
                  .map { map -> tuple(map.simpleName.replaceAll(/chr/,""), map) }
}

def getHapmapGeneticMap() {
    return channel.fromFilePairs( params.panel_dir + "/hapmap/genetic_map_chr*_combined_b37.txt", size: 1 )
                  .map { chr, map_file ->
                         tuple( chr.replaceFirst(/genetic_map_chr/,""), map_file.first() )
                  }
}

def getPlinkGeneticMap() {
    return channel.fromFilePairs( params.panel_dir + '/plinkmap/plink.chr*.GRCh37.map', size: 1)
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

process checkStrand() {
    tag "processing chr${chrom}"
    label 'beagle'
    label 'alignGenotypes'
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
            val("${chrom}"), \
            path("chr${chrom}-aligned.vcf.gz"), \
            path("chr${chrom}-aligned.log")
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
            path("chr${chrom}.${params.out_prefix}.eagle2.noref.vcf.gz")
    script:
        """
        eagle \
          --vcf=${input_vcf} \
          --geneticMapFile=${geneticMap} \
          --chrom=${chrom} \
          --numThreads=${task.cpus} \
          --Kpbwt=${params.kpbwt} \
          --vcfOutFormat=z \
          --outPrefix=chr${chrom}.${params.out_prefix}.eagle2.noref 2>&1 | \
          tee chr${chrom}.${params.out_prefix}.eagle2.noref.log
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
            path("chr${chrom}.${params.out_prefix}.eagle2.vcf.gz")
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
          --outPrefix=chr${chrom}.${params.out_prefix}.eagle2 2>&1 | \
          tee chr${chrom}.${params.out_prefix}.eagle2.log
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
            path("chr${chrom}.${params.out_prefix}.shapeit4.noref.vcf.gz")
    script:
        """
        shapeit4 \
          --input ${input_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --pbwt-depth ${params.pbwt} \
          --log chr${chrom}.${params.out_prefix}.shapeit4.noref.log \
          --output chr${chrom}.${params.out_prefix}.shapeit4.noref.vcf.gz
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
            path("chr${chrom}.${params.out_prefix}.shapeit4.vcf.gz")
    script:
        """
        shapeit4 \
          --input ${input_vcf} \
          --reference ${ref_vcf} \
          --map ${geneticMap} \
          --region ${chrom} \
          --thread ${task.cpus} \
          --pbwt-depth ${params.pbwt} \
          --log chr${chrom}.${params.out_prefix}.shapeit4.log \
          --output chr${chrom}.${params.out_prefix}.shapeit4.vcf.gz
        """
}

process getm3vcf() {
    tag "processing chr${chrom}"
    label 'minimac3'
    label 'phaseGenotypes'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(input_vcf), \
            path(vcf_index) 
    output:
        publishDir path: "${params.panel_dir}/m3vcfs/", mode: 'copy'
        tuple \
            val(chrom), \
            path("chr${chrom}*")
    script:
        """
        Minimac3 \
          --refHaps ${input_vcf} \
          --processReference \
          --prefix chr${chrom} \
          --cpus ${task.cpus} \
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
    label 'bcftools'
    label 'mediumMemory'
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

process imputeVariantsWithMinimac4() {
    tag "processing chr${chrom}"
    label 'minimac4'
    label 'mem_minimac4'
    cache 'lenient'
    input:
        tuple \
            val(chrom), \
            path(target_vcf), \
            path(target_vcf_index), \
            path(erate), \
            path(ref_vcf), \
            path(rec)
    output:
        publishDir path: "${params.out_dir}/imputed/", mode: 'copy'
        tuple \
            val(chrom), \
            path("chr${chrom}.dose.vcf.gz"), \
            path("chr${chrom}.info.gz"), \
            path("chr${chrom}.logfile")
    script:
        """
        minimac4 \
          --refHaps ${ref_vcf} \
          --haps ${target_vcf} \
          --log chr${chrom} \
          --prefix chr${chrom} \
          --cpus ${task.cpus} \
          --allTypedSites \
          --ChunkLengthMb 5 \
          --ChunkOverlapMb 0.0025 \
          --noPhoneHome

        bgzip -f chr${chrom}.info
        """
}
