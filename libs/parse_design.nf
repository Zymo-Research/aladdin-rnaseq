// Function to parse design file and set appropriate channels
def parse_design(LinkedHashMap row) {
    def meta = [:]
    meta.name         = row.sample
    meta.single_end   = !row.read_2

    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(row.read_1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.read_1, checkIfExists: true), file(row.read_2, checkIfExists: true) ] ]
    }
    return array
}