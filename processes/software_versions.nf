// Gather software versions
params.publish_dir = "pipeline_info"
params.pipeline_version = "2.1.0"
params.nextflow_version = "20.07.1"

process software_versions {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: { it == "software_versions.csv" ? it : null }

    input:
    path version_files

    output:
    path 'software_versions_mqc.yaml', emit: report
    path 'software_versions.csv'

    script:
    """
    echo $params.pipeline_version > v_pipeline.txt
    echo $params.nextflow_version > v_nextflow.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}