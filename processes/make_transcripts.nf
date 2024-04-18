// Make a transcripts FASTA file from genome FASTA file and GTF file

process make_transcripts {

    container 'quay.io/biocontainers/rsem:1.3.3--pl526ha52163a_0'

    input:
    path genome
    path gtf

    output:
    path "reference.transcripts.fa", emit: transcripts

    script:
    """
    rsem-prepare-reference --gtf $gtf $genome reference
    """
}