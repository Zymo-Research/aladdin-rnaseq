// Extract various settings associated with a genome

import groovy.json.JsonSlurper

def parse_genome(genome, genomes_path, public_genomes_only) {

    json_slurper = new JsonSlurper()

    // load the genomes from the genome JSON file
    try {
        genomes_file = file(genomes_path)
        genomes = json_slurper.parseText(genomes_file.text)
    } catch (Exception e) {
         exit 1, "Trouble reading file at ${genomes_path}: $e"
    }
    // Check if genome exists in the config file
    if (!(genome in genomes.keySet())) {
        exit 1, "The provided genome '${genome}' is not found in the genomes settings file ${genomes_path}."
    }
    // Check if selected genome is public and if public_genomes_only is set
    if (public_genomes_only && !(genomes[genome]['public'])) {
        exit 1, "'${genome}' is a private genome resource. Try using one of the public genomes in ${genomes_path} or use your own genomes setting file."
    }
    // Set undefined items to false
    genome_settings = genomes[genome]
    keys = ['star','gtf','bed12','gprofiler','ensembl_web','rRNA_gtf','bacteria','csi_index']
    for (key in keys) {
        if (!(key in genome_settings.keySet())) {
            genome_settings[key] = false
        }
    }
    return genome_settings
}