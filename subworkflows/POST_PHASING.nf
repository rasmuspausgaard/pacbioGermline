#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { 
        kivvi_d4z4;
        kivvi05_d4z4;
        pbCPGtools;
        paraphase;
        paraphase35;
        paraphase4;
        starphase;
        methBat;
        methBatNEW_profile_single;
        methBatNEW_pileup;
        multiQC;
        multiQC_ALL;
        mosdepthROI;
        cramino;
        nanoStat;
        whatsHap_stats;
        svTopo;
        svTopo_filtered;
        mitorsaw;
        svdb_SawFish;
        sawFish2_jointCall_all;
        svdb_sawFish2_jointCall_all;
        sawFish2_jointCall_caseID;
        svdb_sawFish2_jointCall_caseID;
        } from "../modules/dnaModules.nf" 

process POST_PHASING_DONE {
    tag "post_phasing_done"
    label "low"

    input:
    val trigger

    output:
    path ".post_phasing.done", emit: done

    script:
    """
    touch .post_phasing.done
    """
}

workflow POST_PHASING {

    take:
    phasedAll
    sawfish_supporting_reads
    mosdepth
    nanoStat

    main:
        mirror_items_ch = Channel.empty()
       // pbCPGtools(phasedAll)
       // methBat(pbCPGtools.out)
        methBatNEW_pileup(phasedAll)
        methBatNEW_profile_single(methBatNEW_pileup.out.met5mC)
        cramino(phasedAll)
        mitorsaw(phasedAll)
        whatsHap_stats(phasedAll)
        //paraphase(phasedAll)
        //paraphase35(phasedAll)
        paraphase4(phasedAll)
        kivvi_d4z4(phasedAll)
        //kivvi05_d4z4(phasedAll)
        starphase(phasedAll)
        svTopo(phasedAll)
        svdb_SawFish(phasedAll)

        mirror_items_ch = mirror_items_ch.mix(methBatNEW_pileup.out[0].map { meta, met_files, bedgraph_files -> tuple(meta, 'specialAnalysis/methylation/5mC_bedgraphs', bedgraph_files) })
        mirror_items_ch = mirror_items_ch.mix(methBatNEW_pileup.out.met5mC.map { meta, bed, tbi -> tuple(meta, 'specialAnalysis/methylation/5mC_pileup', [bed, tbi]) })
        mirror_items_ch = mirror_items_ch.mix(methBatNEW_profile_single.out.map { meta, profile -> tuple(meta, 'specialAnalysis/methylation/5mC_profile', profile) })
        mirror_items_ch = mirror_items_ch.mix(cramino.out.map { meta, txt -> tuple(meta, 'QC/cramino', txt) })
        mirror_items_ch = mirror_items_ch.mix(mitorsaw.out.map { meta, files -> tuple(meta, 'specialAnalysis/mitochondrialVariants', files) })
        mirror_items_ch = mirror_items_ch.mix(whatsHap_stats.out.multiqc.map { meta, tsv -> tuple(meta, 'QC/whatsHap', tsv) })
        mirror_items_ch = mirror_items_ch.mix(paraphase4.out[0].map { meta, files -> tuple(meta, 'specialAnalysis/paraphase4', files) })
        mirror_items_ch = mirror_items_ch.mix(paraphase4.out[1].map { meta, files -> tuple(meta, 'specialAnalysis/paraphase4', files) })
        mirror_items_ch = mirror_items_ch.mix(kivvi_d4z4.out.map { meta, outdir -> tuple(meta, 'repeatExpansions/Kivvi_D4Z4_contraction', outdir) })
        mirror_items_ch = mirror_items_ch.mix(starphase.out.map { meta, files -> tuple(meta, 'specialAnalysis/starphase', files) })
        mirror_items_ch = mirror_items_ch.mix(svTopo.out.map { meta, outdir -> tuple(meta, 'structuralVariants/SVtopo', outdir) })
        mirror_items_ch = mirror_items_ch.mix(svdb_SawFish.out[0].map { meta, files -> tuple(meta, 'structuralVariants', files) })

        /*
            hiPhase_OUT.hiphase_bam
            .join(svdb_SawFish.out.sawfishAF10)
            .join(sawfish_supporting_reads)
            | map {meta,bam,bai,sv10_vcf,sv10_idx,sv_jsonReads -> 
            tuple(meta,[bam:bam,bai:bai,sawfish10_vcf:sv10_vcf,sawfish10_idx:sv10_idx,sawfish_reads:sv_jsonReads])}
            |set {phasedSawfishAF10}   
        */

        phasedAll
        .join(svdb_SawFish.out.sawfishAF10)
        .join(sawfish_supporting_reads)
        | map { meta, data, sv10_vcf, sv10_idx, sv_jsonReads ->
            tuple(meta, [
                bam:            data.bam,
                bai:            data.bai,
                sawfish10_vcf:  sv10_vcf,
                sawfish10_idx:  sv10_idx,
                sawfish_reads:  sv_jsonReads
            ])
        }
        | set { phasedSawfishAF10 }


        svTopo_filtered(phasedSawfishAF10)
        mirror_items_ch = mirror_items_ch.mix(svTopo_filtered.out.map { meta, outdir -> tuple(meta, 'structuralVariants/SVtopo_filtered', outdir) })

        if (!params.skipQC) {
            Channel.empty()
            .mix(mosdepth)
            .mix(nanoStat)
            .mix(whatsHap_stats.out.multiqc)
            .map { meta, qcfile ->
                tuple(params.multiqcKey(meta), meta, qcfile)
            }
            .groupTuple(by: 0)
            .map { key, metas, qcfiles ->

                // pick one representative meta for publishDir + naming
                def meta0 = metas.find { it.relation == 'index' } ?: metas[0]

                tuple(meta0, qcfiles)
            }
            .set { multiqc_inputs_ch }
            multiQC(multiqc_inputs_ch)
            mirror_items_ch = mirror_items_ch.mix(multiQC.out.map { meta, html -> tuple(meta, 'QC', html) })
        }

        post_phasing_done_inputs_ch = Channel.empty()
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(methBatNEW_profile_single.out.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(cramino.out.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(mitorsaw.out.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(whatsHap_stats.out.multiqc.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(paraphase4.out[0].map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(paraphase4.out[1].map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(kivvi_d4z4.out.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(starphase.out.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(svTopo.out.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(svdb_SawFish.out.sawfishAF10.map { 1 })
        post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(svTopo_filtered.out.map { 1 })

        if (!params.skipQC) {
            post_phasing_done_inputs_ch = post_phasing_done_inputs_ch.mix(multiQC.out.map { 1 })
        }

        POST_PHASING_DONE(post_phasing_done_inputs_ch.collect())

    emit:
    done = POST_PHASING_DONE.out.done
    mirror_items = mirror_items_ch
}
