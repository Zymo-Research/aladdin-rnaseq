// Collect deduplication results for multiqc report

process parse_dedup_stat {

    input:
    path deduplog_files

    output:
    path "*mqc.json", emit: dedup_stats

    script:
    """
    parse_umitoolsdedup_for_multiqc.py $deduplog_files
    """
}