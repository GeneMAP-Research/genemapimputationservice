#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
    getChromosomes;
    listChromosomes;
    getHapmapGeneticMap;
    getPlinkGeneticMap;
    getEagleHapmapGeneticMap;
    getShapeitGeneticMap;
    getKgpPanel;
    getVcf;
    //getVcfFileset;
    splitVcfByChrom;
    beaglephase;
    eaglePhaseWithoutRef;
    eaglePhaseWithRef;
    shapeitPhaseWithRef;
    shapeitPhaseWithoutRef;
    getVcfIndex;
    getVcfIndex as getPhasedVcfIndex;
    getCustomVcfPanel;
    getm3vcf;
    getCustomM3Panel;
    imputeVariantsWithMinimac4;
} from "${projectDir}/modules/phasing_and_imputation.nf"

include {
    getPhasedVcf;
    validateVcf;
} from "${projectDir}/modules/custom_panel.nf"

workflow {

    if(params.phase == true && params.impute == false) {
       println "\nMODE: PHASE ONLY\n"
   
       getChromosomes()
           .set { chromosome }
       
       vcf = getVcf()
       chromosome.combine(vcf).set { split_vcf_input }
       per_chr_vcf = splitVcfByChrom(split_vcf_input)
       vcf_fileset = getVcfIndex(per_chr_vcf)
       vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }

       if(params.phase_tool == 'eagle2') {
            if(params.with_ref == true) {
                geneticmap = getEagleHapmapGeneticMap()
                refpanel = getKgpPanel()
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
                refpanel = getKgpPanel()
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
    } else if(params.phase == false && params.impute == true) {
        println "\nMODE: IMPUTE ONLY\n"
        
        getVcf().set { vcf }
        getChromosomes().set { chromosome }
        chromosome.combine(vcf).set { split_vcf_input }
        splitVcfByChrom(split_vcf_input).set { per_chr_vcf }
        getVcfIndex(per_chr_vcf).set { vcf_fileset }
        vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }

        if(params.impute_tool == 'minimac4') {
            vcf = getPhasedVcf()
            validateVcf(vcf).map { chr, vcf_file, vcf_index -> tuple(chr.baseName, vcf_file, vcf_index) }.set { vcf_fileset }
            getCustomM3Panel().set{ minimac_ref_panel }
            vcf_fileset.join(minimac_ref_panel).set{ minimac_input }
            imputeVariantsWithMinimac4(minimac_input)
        }
    } else if(params.phase == true && params.impute == true) {
        println "\nMODE: PHASE AND IMPUTE\n"

        getVcf().set { vcf }
        getChromosomes().set { chromosome }
        chromosome.combine(vcf).set { split_vcf_input }
        splitVcfByChrom(split_vcf_input).set { per_chr_vcf }
        getVcfIndex(per_chr_vcf).set { vcf_fileset }
        vcf_fileset.map { chr, vcf, index -> tuple("${chr}", vcf, index) }.set { vcfFileset }

        if(params.phase_tool == 'eagle2' && params.impute_tool == 'minimac4') {
             if(params.with_ref == true) {
                 getEagleHapmapGeneticMap().set { geneticmap }
                 getKgpPanel().set { refpanel }
                 refpanel.combine(geneticmap).set { panel_map }
                 vcfFileset.join(panel_map).set { phase_input }
                 eaglePhaseWithRef(phase_input).set { phased }
             } else {
                 getEagleHapmapGeneticMap().set { geneticmap}
                 vcfFileset.combine(geneticmap).set { phase_input }
                 eaglePhaseWithoutRef(phase_input).set {phased }
             }

             getPhasedVcfIndex(phased).set { phased_vcf_fileset }
             getCustomM3Panel().set{ minimac_ref_panel }
             phased_vcf_fileset.join(minimac_ref_panel).set{ minimac_input }
             imputeVariantsWithMinimac4(minimac_input)

        } else if(params.phase_tool == 'shapeit4' && params.impute_tool == 'minimac4') {
             if(params.with_ref == true) {
                 geneticmap = getShapeitGeneticMap()
                 refpanel = getKgpPanel()
                 refpanel.join(geneticmap).set { panel_map }
                 vcfFileset.join(panel_map).set { phase_input }
                 phased = shapeitPhaseWithRef(phase_input).view()
             } else {
                 geneticmap = getShapeitGeneticMap()
                 vcfFileset.join(geneticmap).set { phase_input }
                 shapeitPhaseWithoutRef(phase_input).view().set { phased }
             }

             getPhasedVcfIndex(phased).set { phased_vcf_fileset }
             getCustomM3Panel().set{ minimac_ref_panel }
             phased_vcf_fileset.join(minimac_ref_panel).set{ minimac_input }
             imputeVariantsWithMinimac4(minimac_input)

        } else if(params.phase_tool == 'beagle5') {
           geneticmap = getPlinkGeneticMap()
        }
    } else {
        error: println "\nWORKFLOW STOPPED: Please select a run mode - 'phase' and/or 'impute' -\n"
    }
} 
