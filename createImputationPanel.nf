#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getReferencePanel;
    getVcf;
    splitVcfByChrom;
    alignGenotypesToReference;
    beaglephase;
    eaglePhaseWithoutRef;
    getEagleHapmapGeneticMap;
    getVcfIndex;
    createLegendFile;
} from "${projectDir}/modules/phasing.nf"

workflow {
    chromosome = getChromosomes()
    vcf = getVcf()
    chromosome.combine(vcf).set { split_vcf_input }
    per_chr_vcf = splitVcfByChrom(split_vcf_input)
    vcf_fileset = getVcfIndex(per_chr_vcf)

    hapmapGeneticMap = getEagleHapmapGeneticMap()
    vcf_fileset.combine(hapmapGeneticMap).set { eagle_no_ref_phase_input }
    phased_vcf = eaglePhaseWithoutRef(eagle_no_ref_phase_input)
    createLegendFile(phased_vcf).view()

/*
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

    phasedVcfFiles = phaseGenotypes( phasing_input )
    getVcfIndex(phasedVcfFiles).view()

    prePhasingQualityReports = getCheckStrandReports()
*/
	
}
