def getPhasedVcf() {
    return channel.fromFilePairs(params.input_dir + "*.{vcf.gz,vcf.gz.tbi}", checkIfEmpty: true)
                  .map { vcf_name, vcf_fileset -> tuple(vcf_name, vcf_fileset.first(), vcf_fileset.last()) }
                  
}

def getValidationExitStatus(exit_code) {
    if( exit_code == "success" ) {
       println"\nExit status:" + exit_code + "\n"
       println "\nSuccess! Valid VCF file\n"
    } else {
       println"\nExit status: " + exit_code + "\n"
       error: println "\nERROR: Invalid VCF file! Please check that only single chromosome VCF files are provided.\n"
    }
}

process validateVcf() {
    tag "processing file ${vcf_input}"
    label 'bcftools'
    label 'mediumMemory'
    cache 'lenient'
    input:
        tuple \
            val(vcf_name), \
            path(vcf_input), \
            path(vcf_index)
    output:
        tuple \
            path("*"), \
            path(vcf_input), \
            path(vcf_index)
    script:
        //exit_code = 'success'
        """
        number_of_chrs=\$(bcftools query -f '%CHROM\n' ${vcf_input} | sort | uniq | wc -l)
        if [ \${number_of_chrs} -eq 1 ]; then
            chr_number=\$(bcftools query -f '%CHROM\n' ${vcf_input} | sort | uniq)
            echo \${chr_number} > \${chr_number}
            exit_code='success'
        else
            exit_code='failure'
            exit 1
        fi
        echo \${exit_code}
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
          --rsid
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

