// Function to parse design file and set appropriate channels
def parse_design(LinkedHashMap row, ignore_r1) {
    def meta = [:]
    meta.name         = row.sample
    meta.single_end   = ignore_r1 ?: !row.read_2

    def array = []
    if (meta.single_end) {
        if (ignore_r1) {
            // When ignoring Read 1, Read 2 is passed on as "Read 1" for single-end data processing.
            array = [ meta, [ file(row.read_2, checkIfExists: true) ] ] 
        } else {
            array = [ meta, [ file(row.read_1, checkIfExists: true) ] ]
        }
    } else {
        array = [ meta, [ file(row.read_1, checkIfExists: true), file(row.read_2, checkIfExists: true) ] ]
    }
    return array
}