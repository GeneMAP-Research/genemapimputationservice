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
    alignGenotypesToReference;
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

    if(params.phase == true) {
       println "\nPHASE MODE ON!\n"

       chromosome = getChromosomes()
       vcf = getVcf()
       chromosome.combine(vcf).set { split_vcf_input }
       per_chr_vcf = splitVcfByChrom(split_vcf_input)
       vcf_fileset = getVcfIndex(per_chr_vcf)
       vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }

       if(params.phase_tool == 'eagle2') {
            if(params.with_ref == true) {
                geneticmap = getEagleHapmapGeneticMap()
                refpanel = getThousandGenomesReference()
                refpanel.combine(geneticmap).set { panel_map }
                vcfFileset.join(panel_map).set { phase_input }
                phased = eaglePhaseWithRef(phase_input)
            } else {
                geneticmap = getEagleHapmapGeneticMap()
                vcfFileset.combine(geneticmap).set { phase_input }
                phased = eaglePhaseWithoutRef(phase_input)
            }
       } else if(params.phase_tool == 'shapeit4') {
            if(params.with_ref == true) {
                geneticmap = getShapeitGeneticMap()
                refpanel = getThousandGenomesReference()
                refpanel.join(geneticmap).set { panel_map }
                vcfFileset.join(panel_map).set { phase_input }
                phased = shapeitPhaseWithRef(phase_input).view()
            } else {
                geneticmap = getShapeitGeneticMap()
                vcfFileset.join(geneticmap).set { phase_input }
                shapeitPhaseWithoutRef(phase_input).view().set { phased }
          }
       } else if(params.phase_tool == 'beagle5') {
          geneticmap = getPlinkGeneticMap()
       }

    if(params.impute == true) {
        println "\nMODE: IMPUTE\n"
        if(params.impute_tool == 'minimac4') {
            getPhasedVcfIndex(phased).set { vcf_fileset }
            getMinimacReference().set{ minimac_ref_panel }
            vcf_fileset.join(minimac_ref_panel).set{ minimac_input }
            imputeVariantsWithMinimac4(minimac_input).view()
        }
    }

    } else if(params.impute == true) {
        if(params.impute_tool == 'minimac4') {
            println "\nMODE: IMPUTE ONLY\n"
            vcf = getPhasedVcf()
            validateVcf(vcf).map { chr, vcf_file, vcf_index -> tuple(chr.baseName, vcf_file, vcf_index) }.set { vcf_fileset }
            getMinimacReference().set{ minimac_ref_panel }
            vcf_fileset.join(minimac_ref_panel).set{ minimac_input }
            imputeVariantsWithMinimac4(minimac_input)
        }
    } else if(params.impute == false) {
       error: println "\nWORKFLOW STOPPED: Please select a run mode - 'phase' and/or 'impute' -\n"
    } 



/*
    chromosome = getChromosomes()

    vcf = getVcf()
    chromosome.combine(vcf).set { split_vcf_input }
    per_chr_vcf = splitVcfByChrom(split_vcf_input)
    vcf_fileset = getVcfIndex(per_chr_vcf)
    plinkGeneticMap = getPlinkGeneticMap()
    eagleGeneticMap = getEagleHapmapGeneticMap().view()

    refPanel = getCostumReferencePanel()
    vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }
    refPanel.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { ref_panel }
    vcfFileset.join(ref_panel).set { phase_input }
    phase_input.join(plinkGeneticMap).set { impute_input }

    beaglephase(impute_input)


    thousandGenomesReference = getThousandGenomesReference().view()

    //vcf = getVcf()

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

    phasedVcfFiles = phaseGenotypes( phasing_input )
    getVcfIndex(phasedVcfFiles).view()

    prePhasingQualityReports = getCheckStrandReports()
*/
}
