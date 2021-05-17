// Check STAR log and filter those below a cutoff of alignment rate
def check_star_log(star_log, cutoff, log) {
    def percent_aligned = 0
    star_log.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = star_log.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() < cutoff.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}