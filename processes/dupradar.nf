// Run dupRadar
params.publish_dir = "dupRadar"
params.skip_qc = false
params.skip_dupradar = false
params.strandedness = 0

process dupradar {
    label 'mid_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_duprateExpDens.pdf") > 0) "scatter_plots/$filename"
            else if (filename.indexOf("_duprateExpBoxplot.pdf") > 0) "box_plots/$filename"
            else if (filename.indexOf("_expressionHist.pdf") > 0) "histograms/$filename"
            else if (filename.indexOf("_dupMatrix.txt") > 0) "gene_data/$filename"
            else if (filename.indexOf("_intercept_slope.txt") > 0) "intercepts_slopes/$filename"
            else null
        }

    when:
    !params.skip_qc && !params.skip_dupradar

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf

    output:
    path "*.pdf"
    path "*_dupMatrix.txt"
    path "*_intercept_slope.txt"
    path "*_mqc.txt", emit: report
    path "v_dupRadar.txt", emit: version

    script:
    paired = meta.single_end ? 'single' : 'paired'
    """
    dupRadar.r $bam $gtf ${params.strandedness} $paired ${task.cpus}
    Rscript -e "write(x=as.character(packageVersion('dupRadar')), file='v_dupRadar.txt')"
    """
}