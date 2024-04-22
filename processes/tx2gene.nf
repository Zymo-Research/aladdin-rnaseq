// Prepare a transcript-to-gene mapping file from GTF file
params.fc_extra_attributes = 'gene_name'

process tx2gene {
    input:
    path gtf

    output:
    path 'tx2gene.tsv', emit: mapping

    script:
    extra_attributes = params.fc_extra_attributes ? "--extra_attributes ${params.fc_extra_attributes}" : ""
    """
    tx2gene.py $gtf $extra_attributes
    """
}