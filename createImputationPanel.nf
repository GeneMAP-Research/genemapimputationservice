#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getVcf;
    splitVcfByChrom;
    beaglephase;
    eaglePhaseWithoutRef;
    getEagleHapmapGeneticMap;
    getVcfIndex;
} from "${projectDir}/modules/phasing_and_imputation.nf"

include {
    getPhasedVcf;
    validateVcf;
    getValidationExitStatus;
    createLegendFile;
    getm3vcf;
    prepareChrXPanel;
} from "${projectDir}/modules/custom_panel.nf"

workflow {
    chromosome = getChromosomes()
    vcf = getPhasedVcf()
    validateVcf(vcf).map { chr, vcf_file -> tuple(chr.baseName, vcf_file) }.set { chrom_vcf }
    getVcfIndex(chrom_vcf).view().set { vcf_fileset }

/*
    hapmapGeneticMap = getEagleHapmapGeneticMap()
    vcf_fileset.combine(hapmapGeneticMap).set { eagle_no_ref_phase_input }
    phased_vcf = eaglePhaseWithoutRef(eagle_no_ref_phase_input)
    createLegendFile(phased_vcf).view()
*/

    m3vcf = getm3vcf(vcf_fileset).view()

}
