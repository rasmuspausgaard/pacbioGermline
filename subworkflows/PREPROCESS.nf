#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include {pbmm2_align;
        create_fofn;
        inputFiles_symlinks_ubam;
        pbmm2_align_mergedData;
        extractHifi;
        } from "../modules/dnaModules.nf" 

workflow PREPROCESS {

    take:
    finalUbamInput     
   
    main:

    inputFiles_symlinks_ubam(finalUbamInput)
    create_fofn(finalUbamInput)
    pbmm2_align_mergedData(create_fofn.out)

    mirror_items_ch = channel.empty()
    mirror_items_ch = mirror_items_ch.mix(inputFiles_symlinks_ubam.out.map { meta, data -> tuple(meta, 'documents/inputSymlinks', data) })
    mirror_items_ch = mirror_items_ch.mix(create_fofn.out.map { meta, fofn -> tuple(meta, 'documents', fofn) })

    if (!params.failedReads && !params.allReads && !params.hifiReads) {
        extractHifi(pbmm2_align_mergedData.out.bamAll)
        extractHifi.out.alignedHifi
            .join(pbmm2_align_mergedData.out.bamAll)
            | map {meta,bamHifi,baiHifi,bamAll,baiAll ->
            tuple(meta, [mainBamFile:bamHifi, mainBaiFile:baiHifi, bamAll:bamAll, baiAll:baiAll])}
            | set {alignedFinal_ch}
    }
    if (params.allReads || params.hifiReads || params.failedReads) {
        pbmm2_align_mergedData.out.bamAll
            | map {meta,bamAll,baiAll ->
            tuple(meta,[mainBamFile:bamAll,mainBaiFile:baiAll])}
            | set {alignedFinal_ch}
    }

    emit:
    alignedFinal = alignedFinal_ch
    mirror_items = mirror_items_ch

}
