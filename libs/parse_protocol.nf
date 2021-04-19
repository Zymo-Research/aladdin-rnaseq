// Extract various settings associated with a protocol

import groovy.json.JsonSlurper

def parse_protocol(protocol, protocols_path) {

    json_slurper = new JsonSlurper()

    // load the protocols from the protocol JSON file
    try {
        protocols_file = file(protocols_path)
        protocols = json_slurper.parseText(protocols_file.text)
    } catch (Exception e) {
         exit 1, "Trouble reading file at ${protocols_path}: $e"
    }
    if (!(protocol in protocols.keySet())) {
        exit 1, "The provided protocol '${protocol}' is not found in the protocols settings file ${protocols_path}."
    }
    return protocols[protocol]
}