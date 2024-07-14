#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getEagleHapmapGeneticMap;
    getShapeitGeneticMap;
    getThousandGenomesReference;
    getVcf;
    splitVcfByChrom;
    sortIndexVcf;
    checkStrand;
    concatVcfs;
    beaglephase;
    eaglePhaseWithoutRef;
    eaglePhaseWithRef;
    shapeitPhaseWithRef;
    shapeitPhaseWithoutRef;
    getVcfIndex;
    removeDupVarsAndIndexVcf;
    getVcfIndex as getPhasedVcfIndex;
    getCostumReferencePanel;
    getm3vcf;
    getMinimacReference;
    imputeVariantsWithMinimac4;
} from "${projectDir}/modules/phasing_and_imputation.nf"

include {
    getPhasedVcf;
    validateVcf;
} from "${projectDir}/modules/custom_panel.nf"

workflow {

    if(params.autosome == true) {
       getChromosomes()
          .filter(Number)
          .set { chromosome }
    } else {
       getChromosomes()
          .set { chromosome }
    }

    getVcf()
       .set { vcf }

    chromosome
       .combine(vcf)
       .set { split_vcf_input }

    getThousandGenomesReference()
       .set { thousandGenomesReference }

    getHapmapGeneticMap()
       .set { geneticMap }

    splitVcfByChrom(split_vcf_input)
       .set { per_chr_vcf }

    removeDupVarsAndIndexVcf(per_chr_vcf)
       .set { vcf_fileset }

    vcf_fileset
       .map { chr, vcf, index -> tuple("${chr}", vcf, index) }
       .set { vcfFileset }

    vcfFileset
	.map { chrom, vcfFile, vcfIndex -> tuple( "${chrom}", vcfFile, vcfIndex ) }
	.join( thousandGenomesReference )
	.join( geneticMap )
	.set { checkstrand_input }

    alignedVcfs = checkStrand( checkstrand_input ).map { chr, vcf, log -> tuple(chr, vcf) }
    aligned_indexed_vcfs = getVcfIndex(alignedVcfs).map { chr, vcf, index -> tuple(vcf, index) }.collect()
    concatVcfs(aligned_indexed_vcfs)

//    getPlinkGeneticMap()
//       .set { plinkGeneticMap }
//    getEagleHapmapGeneticMap()
//       .set { eagleGeneticMap }
//
//    refPanel = getCostumReferencePanel()
//    vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }
//    refPanel.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { ref_panel }
//    vcfFileset.join(ref_panel).set { phase_input }
//    phase_input.join(plinkGeneticMap).set { impute_input }
//
//    beaglephase(impute_input)
}
