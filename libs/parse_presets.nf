def parse_presets(preset) {
    // Protocol presets
    presets = [ "zymo_ribofree" : [ "clip_r1":0, "clip_r2":10, "three_prime_clip_r1":0, "three_prime_clip_r2":0,  
                                    "adapter":"NNNNNNNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "adapter2":"AGATCGGAAGAGCGTCGTGTAGGGAAAGA", "strandedness":2], 
                "zymo_3mrna"    : ["clip_r1":10, "clip_r2":0, "three_prime_clip_r1":0, "three_prime_clip_r2":0,  
                                    "adapter":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "strandedness":1],
                "zymo_3mrna_nodedup" : [ "clip_r1":10, "clip_r2":0, "three_prime_clip_r1":0, "three_prime_clip_r2":0,  
                                    "adapter":"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "strandedness":1],
                "illumina"      : [ "clip_r1":0, "clip_r2":0, "three_prime_clip_r1":0, "three_prime_clip_r2":0, 
                                    "adapter":false, "adapter2":false, "strandedness":0],
                "pico"          : [ "clip_r1":3, "clip_r2":0, "three_prime_clip_r1":0, "three_prime_clip_r2":3, 
                                    "adapter":false, "adapter2":false, "strandedness":1]
    ]
    if (!(preset in presets.keySet())) {
        exit 1, "The provided protocol '${preset}' is not supported. The supported protocols are ${presets.keySet().join(', ')}"
    }
    return presets[preset]
}