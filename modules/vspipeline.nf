#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
 VSpipeline VarSeq module
================================================================================

Input channel:
  tuple val(meta), val(data)

Expected meta fields:
  meta.rekv      order / rekvisition id
  meta.id        NPN / sample id
  meta.testlist  analysis testlist

Expected data fields:
  data.bam or data.mainBamFile       phased BAM used for VarSeq sample_fields_file
  data.dv_vcf                        phased DeepVariant VCF
  data.sawfish10_vcf or sawfish_vcf  SV/BND VCF for update_bnd_import

The main workflow filters the input channel so this process only receives samples
where meta.testlist normalizes to params.vspipeline_testlist.
*/

params.vspipeline_testlist                  = params.vspipeline_testlist ?: 'SL_NGC_HJERTESYGDOM'
params.vs_sif                               = params.vs_sif ?: '/lnx01_data2/shared/testdata/RunData/vspipeline_latest.sif'
params.vspipeline_template_path             = params.vspipeline_template_path ?: (params.template_path ?: '/appdata/VarSeq/User Data/ProjectTemplates/LRS-WGS Arvelig hjertesygdom v260512.vsproject-template')
params.vspipeline_project_base              = params.vspipeline_project_base ?: '/data/shared/VarSeq/projects/WGS_NGC/ArveligeHjertesygdomme/00_Pacbio'
params.vspipeline_project_base_in_container = params.vspipeline_project_base_in_container ?: '/appdata/projects/WGS_NGC/ArveligeHjertesygdomme/00_Pacbio'

process VSpipeline {

    errorStrategy 'ignore'

    tag { "${meta.rekv}_${meta.id}" }

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path('*_varseq.done'), emit: done

    script:
    def order_id = meta.rekv.toString().padLeft(10, '0')
    def npn_id   = meta.id.toString()

    def template_project_stem = new File(params.vspipeline_template_path)
        .name
        .replaceAll(/\.vsproject-template$/, '')
        .replaceAll(/\s+/, '_')

    def varseq_project_folder = "${order_id}_${npn_id}_${template_project_stem}"
    def varseq_project_name   = "${order_id}_${npn_id}_${template_project_stem}"

    def vcf_path      = data.dv_vcf.toString()
    def bnd_vcf_path  = (data.sawfish10_vcf ?: data.sawfish_vcf).toString()
    def bam_path      = (data.bam ?: data.mainBamFile).toString()
    def sample_file   = "${varseq_project_folder}_sample_fields_file.csv"
    def done_file     = "${varseq_project_folder}_varseq.done"

    """
    set -euo pipefail

    export XDG_RUNTIME_DIR="\$PWD/.xdg-runtime"
    mkdir -p "\$XDG_RUNTIME_DIR"
    chmod 700 "\$XDG_RUNTIME_DIR"

    cat > "${sample_file}" << 'EOF'
Sample,BAM Path
${npn_id},${bam_path}
EOF

    mkdir -p "${params.vspipeline_project_base}/${varseq_project_folder}"

    full_path_in="${params.vspipeline_project_base_in_container}/${varseq_project_folder}/${varseq_project_name}"

    echo "=== VSpipeline input ==="
    echo "sample=${npn_id}"
    echo "testlist=${meta.testlist}"
    echo "vcf_path=${vcf_path}"
    echo "bnd_vcf_path=${bnd_vcf_path}"
    echo "bam_path=${bam_path}"
    echo "sample_file=${sample_file}"
    echo "varseq_project_folder=${varseq_project_folder}"
    echo "varseq_project_name=${varseq_project_name}"
    echo "full_path_in=\${full_path_in}"

    test -f "${vcf_path}"
    test -f "${bnd_vcf_path}"
    test -f "${bam_path}"
    test -f "${sample_file}"

    echo "=== CONTAINER DEBUG ==="
    singularity exec \
      --env XDG_RUNTIME_DIR="\$XDG_RUNTIME_DIR" \
      --bind /data/shared/VarSeq/:/appdata \
      --bind /lnx01_data2:/lnx01_data2 \
      --bind /lnx01_data3:/lnx01_data3 \
      ${params.vs_sif} \
      bash -lc 'test -f '"${vcf_path}"' && echo CONTAINER_VCF_FOUND || echo CONTAINER_VCF_MISSING; \
                test -f '"${bnd_vcf_path}"' && echo CONTAINER_BND_VCF_FOUND || echo CONTAINER_BND_VCF_MISSING; \
                test -f '"${bam_path}"' && echo CONTAINER_BAM_FOUND || echo CONTAINER_BAM_MISSING'

    singularity run \
      --env XDG_RUNTIME_DIR="\$XDG_RUNTIME_DIR" \
      --bind /data/shared/VarSeq/:/appdata \
      --bind /data/shared/VarSeq/:/data/shared/VarSeq \
      --bind /lnx01_data2:/lnx01_data2 \
      --bind /lnx01_data3:/lnx01_data3 \
      ${params.vs_sif} \
      -c login ${params.user_email} ${params.user_login} \
      -c license_activate ${params.license_key} \
      -c project_create "\${full_path_in}" "${params.vspipeline_template_path}" \
      -c import files="${vcf_path}" sample_fields_file="\$PWD/${sample_file}" \
      -c update_bnd_import files="${bnd_vcf_path}" table_id="BND" \
      -c download_required_sources \
      -c task_wait \
      -c get_data_list table \
      -c get_task_list

    touch "${done_file}"
    """
}
