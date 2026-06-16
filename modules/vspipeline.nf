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

  For BND/SV import:
    data.sawfish10_vcf or data.sawfish_vcf

  For CNV import:
    data.sawfish_cnv_vcf or data.cnv_vcf
    fallback: data.sawfish10_vcf or data.sawfish_vcf
*/

params.vs_sif = params.vs_sif ?: '/lnx01_data2/shared/testdata/RunData/vspipeline_latest.sif'

/*
================================================================================
 Testlist-specific VarSeq configuration

Note:
  project_base_host is the host path.
  project_base_container is the corresponding path inside the VarSeq container.

Because /data/shared/VarSeq is mounted as /appdata, this mapping is:

  host:      /data/shared/VarSeq/projects/...
  container: /appdata/projects/...
================================================================================
*/

params.vspipeline_configs = params.vspipeline_configs ?: [

    SL_NGC_HJERTESYGDOM: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/LRS-WGS Arvelig hjertesygdom v260512.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/ArveligeHjertesygdomme/00_Pacbio',
        project_base_container: '/appdata/projects/WGS_NGC/ArveligeHjertesygdomme/00_Pacbio',
        bnd_table_id: 'BND'
    ],

    SL_NGC_UNGE_VOKSNE: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV NPN RUN KRFT LRS-WGS v260529.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/UngeVoksne.ArveligKrft/PacBio/2026',
        project_base_container: '/appdata/projects/WGS_NGC/UngeVoksne.ArveligKrft/PacBio/2026',
        bnd_table_id: 'BND'
    ],

    SL_NGC_ARVELIG_KRFT: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV NPN RUN KRFT LRS-WGS v260529.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/UngeVoksne.ArveligKrft/PacBio/2026',
        project_base_container: '/appdata/projects/WGS_NGC/UngeVoksne.ArveligKrft/PacBio/2026',
        bnd_table_id: 'BND'
    ],

    SL_NGC_NEUROGENETIK: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260529.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/Neurogenetiske sygdomme/LRS/2026',
        project_base_container: '/appdata/projects/WGS_NGC/Neurogenetiske sygdomme/LRS/2026',
        bnd_table_id: 'BND'
    ],

    SL_LWG_GENOM: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260529.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS/PacBio/2026',
        project_base_container: '/appdata/projects/WGS/PacBio/2026',
        bnd_table_id: 'BND'
    ],

    SL_NGC_SJAELDNE: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260530.2.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/Sjaeldne.sygdomme/Pacbio/2026',
        project_base_container: '/appdata/projects/WGS_NGC/Sjaeldne.sygdomme/Pacbio/2026',
        bnd_table_id: 'Breakends'
    ],

    SL_NGC_NYRESVIGT: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260530.2.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/Kronisk nyresvigt/PacBio_LRS/2026',
        project_base_container: '/appdata/projects/WGS_NGC/Kronisk nyresvigt/PacBio_LRS/2026',
        bnd_table_id: 'Breakends'
    ],

    SL_NGC_ENDOKRINOLOG: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260530.2.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/Endokrinologi/2026',
        project_base_container: '/appdata/projects/WGS_NGC/Endokrinologi/2026',
        bnd_table_id: 'Breakends'
    ],

    SL_NGC_OFTALMOLOGI: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260530.2.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/Oftalmologi/PacBio/2026',
        project_base_container: '/appdata/projects/WGS_NGC/Oftalmologi/PacBio/2026',
        bnd_table_id: 'Breakends'
    ],

    SL_NGC_HUDSYGDOM: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN_LRS-WGS v260530.2.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS_NGC/Dermatology/2026',
        project_base_container: '/appdata/projects/WGS_NGC/Dermatology/2026',
        bnd_table_id: 'Breakends'
    ],

    SL_LWG_CNV: [
        template_path: '/appdata/VarSeq/User Data/ProjectTemplates/REKV_NPN LRS CNV v260324.vsproject-template',
        project_base_host: '/data/shared/VarSeq/projects/WGS CNV/2026',
        project_base_container: '/appdata/projects/WGS CNV/2026',
        cnv_table_id: 'Sawfish CNV'
    ]
]


process VSpipeline {

    errorStrategy 'ignore'

    tag { "${meta.rekv}_${meta.id}_${meta.testlist}" }

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path('*_varseq.done'), emit: done

    script:
    def order_id = meta.rekv.toString().padLeft(10, '0')
    def npn_id   = meta.id.toString()

    def testlist_key = (meta.testlist ?: '')
        .toString()
        .trim()
        .replaceAll('-', '_')

    def cfg = params.vspipeline_configs[testlist_key]

    if (!cfg) {
        throw new IllegalArgumentException(
            "No VSpipeline config found for testlist='${meta.testlist}' normalized as '${testlist_key}'"
        )
    }

    def template_path          = cfg.template_path
    def project_base_host      = cfg.project_base_host
    def project_base_container = cfg.project_base_container

    /*
    ============================================================================
     Decide whether this testlist should use BND/SV import or CNV import
    ============================================================================
    */

    def is_cnv_import = cfg.containsKey('cnv_table_id')

    def variant_import_command = is_cnv_import
        ? 'update_cnv_import'
        : 'update_bnd_import'

    def variant_table_id = is_cnv_import
        ? cfg.cnv_table_id
        : (cfg.bnd_table_id ?: 'BND')

    /*
    ============================================================================
     Build VarSeq project name
    ============================================================================
    */

    def template_project_stem = new File(template_path)
        .name
        .replaceAll(/\.vsproject-template$/, '')
        .replaceAll(/\s+/, '_')

    def varseq_project_folder = "${order_id}_${npn_id}_${template_project_stem}"
    def varseq_project_name   = "${order_id}_${npn_id}_${template_project_stem}"

    /*
    ============================================================================
     Input files
    ============================================================================
    */

    def vcf_path = data.dv_vcf.toString()

    def variant_vcf_candidate = is_cnv_import
        ? (data.sawfish_cnv_vcf ?: data.cnv_vcf ?: data.sawfish10_vcf ?: data.sawfish_vcf)
        : (data.sawfish10_vcf ?: data.sawfish_vcf)

    if (!variant_vcf_candidate) {
        throw new IllegalArgumentException(
            "Missing variant VCF for testlist='${meta.testlist}'. " +
            "For CNV import expected one of: data.sawfish_cnv_vcf, data.cnv_vcf, data.sawfish10_vcf, data.sawfish_vcf. " +
            "For BND import expected one of: data.sawfish10_vcf, data.sawfish_vcf."
        )
    }

    def variant_vcf_path = variant_vcf_candidate.toString()

    def bam_candidate = data.bam ?: data.mainBamFile

    if (!bam_candidate) {
        throw new IllegalArgumentException(
            "Missing BAM for testlist='${meta.testlist}'. Expected data.bam or data.mainBamFile."
        )
    }

    def bam_path    = bam_candidate.toString()
    def sample_file = "${varseq_project_folder}_sample_fields_file.csv"
    def done_file   = "${varseq_project_folder}_varseq.done"

    """
    set -euo pipefail

    export XDG_RUNTIME_DIR="\$PWD/.xdg-runtime"
    mkdir -p "\$XDG_RUNTIME_DIR"
    chmod 700 "\$XDG_RUNTIME_DIR"

    cat > "${sample_file}" << EOF
Sample,BAM Path
${npn_id},${bam_path}
EOF

    mkdir -p "${project_base_host}/${varseq_project_folder}"

    full_path_in="${project_base_container}/${varseq_project_folder}/${varseq_project_name}"

    echo "=== VSpipeline VarSeq DEBUG ==="
    echo "order_id=${order_id}"
    echo "npn_id=${npn_id}"
    echo "testlist=${meta.testlist}"
    echo "testlist_key=${testlist_key}"
    echo "template_path=${template_path}"
    echo "project_base_host=${project_base_host}"
    echo "project_base_container=${project_base_container}"
    echo "varseq_project_folder=${varseq_project_folder}"
    echo "varseq_project_name=${varseq_project_name}"
    echo "full_path_in=\${full_path_in}"
    echo "vcf_path=${vcf_path}"
    echo "variant_vcf_path=${variant_vcf_path}"
    echo "bam_path=${bam_path}"
    echo "sample_file=\$PWD/${sample_file}"
    echo "variant_import_command=${variant_import_command}"
    echo "variant_table_id=${variant_table_id}"

    echo "=== HOST FILE CHECK ==="
    test -f "${vcf_path}" && echo "HOST DV VCF FOUND" || echo "HOST DV VCF MISSING"
    test -f "${variant_vcf_path}" && echo "HOST VARIANT VCF FOUND" || echo "HOST VARIANT VCF MISSING"
    test -f "${bam_path}" && echo "HOST BAM FOUND" || echo "HOST BAM MISSING"
    test -f "\$PWD/${sample_file}" && echo "HOST SAMPLE FILE FOUND" || echo "HOST SAMPLE FILE MISSING"

    singularity run \\
      --env XDG_RUNTIME_DIR="\$XDG_RUNTIME_DIR" \\
      --bind /data/shared/VarSeq/:/appdata \\
      --bind /data/shared/VarSeq/:/data/shared/VarSeq \\
      --bind /lnx01_data2:/lnx01_data2 \\
      --bind /lnx01_data3:/lnx01_data3 \\
      --bind /fast:/fast \\
      ${params.vs_sif} \\
      -c login ${params.user_email} ${params.user_login} \\
      -c license_activate ${params.license_key} \\
      -c project_create "\${full_path_in}" "${template_path}" \\
      -c import files="${vcf_path}" sample_fields_file="\$PWD/${sample_file}" \\
      -c ${variant_import_command} files="${variant_vcf_path}" table_id="${variant_table_id}" \\
      -c download_required_sources \\
      -c task_wait \\
      -c get_data_list table \\
      -c get_task_list

    touch "${done_file}"
    """
}
