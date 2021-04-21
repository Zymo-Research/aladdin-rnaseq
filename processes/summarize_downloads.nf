// Collect and parse information about files for download on the aladdin platform
params.publish_dir = 'download_data'

process summarize_downloads {
    label 'no_cache'
    publishDir "${params.publish_dir}", mode: 'copy'

    input:
    path locations
    path md5sums
    path design

    output:
    path 'files_to_download.json'

    script:
    md5_option = md5sums.size() > 0 ? "-m $md5sums" : ''
    """
    summarize_downloads.py $locations $md5_option -d $design
    """
}