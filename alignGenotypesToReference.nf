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
    checkStrand;
    beaglephase;
    eaglePhaseWithoutRef;
    eaglePhaseWithRef;
    shapeitPhaseWithRef;
    shapeitPhaseWithoutRef;
    getVcfIndex;
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
    getThousandGenomesReference().view()
       .set { thousandGenomesReference }
    getHapmapGeneticMap()
       .set { geneticMap }
    splitVcfByChrom(split_vcf_input)
       .set { per_chr_vcf }
    getVcfIndex(per_chr_vcf)
       .set { vcf_fileset }
    vcf_fileset
       .map { chr, vcf, index -> tuple("${chr}", vcf, index) }
       .set { vcfFileset }

    vcfFileset
	.map { chrom, vcfFile, vcfIndex -> tuple( "${chrom}", vcfFile, vcfIndex ) }
	.join( thousandGenomesReference )
	.join( geneticMap )
	.set { checkstrand_input }

    alignedVcfs = checkStrand( checkstrand_input )

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
