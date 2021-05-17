// Run RSeQC for QC
params.publish_dir = "RSeQC"
params.skip_qc = false
params.skip_rseqc = false
params.csi_index = false
params.bacteria = false

process rseqc {
    label 'mid_memory'
    tag "${meta.name}"
    publishDir "${params.publish_dir}" , mode: 'copy',
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("RPKM_saturation.eRPKM.xls") > 0)         "RPKM_saturation/rpkm/$filename"
            else if (filename.indexOf("RPKM_saturation.rawCount.xls") > 0)      "RPKM_saturation/counts/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.pdf") > 0)    "RPKM_saturation/$filename"
            else if (filename.indexOf("RPKM_saturation.saturation.r") > 0)      "RPKM_saturation/rscripts/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junction_annotation_log.txt") > 0)       "junction_annotation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else null
        }

    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    tuple val(meta), path(bam), path(bai)
    path bed12

    output:
    path "${meta.name}*", emit: report
    path "v_rseqc.txt", emit: version

    script:
    junc_anno = params.bacteria ? '' : "junction_annotation.py -i $bam -o ${meta.name}.rseqc -r $bed12 2> ${meta.name}.rseqc.junction_annotation_log.txt"
    junc_satu = params.bacteria ? '' : "junction_saturation.py -i $bam -o ${meta.name}.rseqc -r $bed12"
    inner_dis = (params.csi_index || meta.single_end) ? '' : "inner_distance.py -i $bam -o ${meta.name}.rseqc -r $bed12"
    read_dis = params.csi_index ? '' : "read_distribution.py -i $bam -r $bed12 > ${meta.name}.rseqc.read_distribution.txt"
    read_dup = (bam.size() > 6000000000) ? '' : "read_duplication.py -i $bam -o ${meta.name}.rseqc.read_duplication"
    """
    infer_experiment.py -i $bam -r $bed12 -s 2000000 > ${meta.name}.rseqc.infer_experiment.txt
    $junc_anno
    bam_stat.py -i $bam > ${meta.name}.rseqc.bam_stat.txt
    $junc_satu
    $inner_dis
    $read_dis
    $read_dup
    read_duplication.py --version &> v_rseqc.txt
    """
}