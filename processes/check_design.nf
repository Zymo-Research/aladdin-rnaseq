// Sanity check design file
process check_design {
    tag "$design"

    input:
    path design
    path comparisons

    output:
    path "checked_${design}", emit: checked_design
    path "DESeq2_${design}", emit: deseq2_design

    script:
    comparison_file = params.comparisons ? "-c $comparisons" : ''
    """
    check_design.py $comparison_file $design
    """
}