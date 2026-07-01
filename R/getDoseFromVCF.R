#' Read genotype dosages from an indexed VCF/BCF file
#'
#' `getDoseFromVCF()` extracts genotype dosages for selected variants from an
#' indexed VCF or BCF file using `bcftools query`. Variants are specified as
#' character strings in the form `CHR:POS:REF:ALT`. The function first extracts
#' records overlapping the requested genomic positions and then keeps only records
#' whose chromosome, position, reference allele, and alternate allele exactly match
#' the requested variants.
#'
#' The function currently reads dosage values from the `FORMAT/DS` field and
#' returns a numeric dosage matrix with variants in rows and samples in columns.
#'
#' @param vcf Character scalar. Path to a bgzip-compressed and indexed VCF file
#'   (`.vcf.gz`) or an indexed BCF file (`.bcf`).
#' @param variants Character vector. Variants to extract, specified as
#'   `CHR:POS:REF:ALT`, for example `"1:12345:A:G"` or `"chr1:12345:A:G"`.
#' @param samples Optional character vector. Sample IDs to extract. If `NULL`,
#'   all samples in the VCF/BCF are extracted.
#' @param BCFTOOLS Character scalar. Path to the `bcftools` executable, or the
#'   command name if `bcftools` is available in the system `PATH`.
#' @param tmpdir Character scalar. Directory used to write temporary files such
#'   as the BED file and optional sample list.
#' @param nline Integer. Number of `bcftools query` output lines read from the
#'   pipe at a time. Smaller values reduce peak memory usage, while larger values
#'   may improve speed.
#' @param progress Logical. If `TRUE`, print progress information while reading
#'   and parsing records.
#'
#' @return A numeric dosage matrix with variants in rows and samples in columns.
#'   Row names are the requested `CHR:POS:REF:ALT` variant IDs, and column names
#'   are sample IDs. Variants not found in the VCF/BCF are returned as rows filled
#'   with `NA`.
#'
#' @details
#' This function requires `bcftools` and an index for the input VCF/BCF file.
#' For `.vcf.gz` files, a tabix index such as `.tbi` or `.csi` is expected.
#' For `.bcf` files, a `.csi` index is expected.
#'
#' The extraction is position-based at the `bcftools` step. Exact allele matching
#' is performed afterwards in R using `CHR:POS:REF:ALT`.
#'
#' @examples
#' \dontrun{
#' variants <- c("chr1:12345:A:G", "chr1:67890:C:T")
#' Gall <- getDoseFromVCF(
#'   vcf = "imputed.vcf.gz",
#'   variants = variants,
#'   BCFTOOLS = "bcftools",
#'   nline = 20,
#'   progress = TRUE
#' )
#'
#' dim(Gall)
#' Gall[1:2, 1:5]
#' }
#'
#' @export
getDoseFromVCF = function(vcf, variants, samples=NULL, BCFTOOLS="bcftools", tmpdir=tempdir(), nline=20, progress=TRUE){
    variants = as.character(variants)
    sp = strsplit(variants, ":", fixed=TRUE)
    vdat = data.frame(CHR=sapply(sp, `[`, 1), POS=as.integer(sapply(sp, `[`, 2)), REF=sapply(sp, `[`, 3), ALT=sapply(sp, `[`, 4), stringsAsFactors=FALSE)

    reg = unique(vdat[,c("CHR","POS")])
    bed = tempfile(tmpdir=tmpdir, fileext=".bed")
    write.table(data.frame(reg$CHR, reg$POS-1L, reg$POS), bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    if(is.null(samples)){
        sample_ids = system2(BCFTOOLS, c("query","-l",vcf), stdout=TRUE)
        sarg = ""
    }else{
        sample_ids = as.character(samples)
        sfile = tempfile(tmpdir=tmpdir)
        writeLines(sample_ids, sfile)
        sarg = paste("-S", shQuote(sfile))
    }

    dose = matrix(NA_real_, nrow=length(variants), ncol=length(sample_ids), dimnames=list(variants, sample_ids))
    idx = seq_along(variants)
    names(idx) = variants

    fmt = "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n"
    errfile = tempfile(tmpdir=tmpdir)
    cmd = paste(shQuote(BCFTOOLS), "query -R", shQuote(bed), sarg, "-f", shQuote(fmt), shQuote(vcf), "2>", shQuote(errfile))

    bar = function(x, n, width=30){
        p = if(n <= 0) 1 else min(1, x/n)
        k = floor(width*p)
        paste0("[", paste(rep("=", k), collapse=""), paste(rep(" ", width-k), collapse=""), "]")
    }
    draw = function(nread, nparse, nhit){
        if(!progress) return(invisible())
        cat(sprintf("\rhit %s %d/%d | read lines %d | parsed lines %d",
                    bar(nhit, length(variants)), nhit, length(variants), nread, nparse))
        flush.console()
    }

    con = pipe(cmd, open="r")
    closed = FALSE
    on.exit({
        if(!closed) try(close(con), silent=TRUE)
        unlink(c(bed, errfile, if(exists("sfile")) sfile))
    }, add=TRUE)

    nread = 0L
    nparse = 0L
    nhit = 0L
    seen = rep(FALSE, length(variants))
    draw(0,0,0)

    repeat{
        lines = readLines(con, n=nline, warn=FALSE)
        if(length(lines)==0) break

        nread = nread + length(lines)
        z = strsplit(lines, "\t", fixed=TRUE)

        for(a in z){
            vid = paste(a[1], a[2], a[3], a[4], sep=":")
            i = unname(idx[vid])
            if(!is.na(i)){
                x = a[-(1:4)]
                x[x=="."] = NA_character_
                dose[i,] = as.numeric(x)
                if(!seen[i]){
                    seen[i] = TRUE
                    nhit = nhit + 1L
                }
            }
            nparse = nparse + 1L
        }

        draw(nread, nparse, nhit)
    }

    status = close(con)
    closed = TRUE
    if(progress) cat("\n")

    err = readLines(errfile, warn=FALSE)
    if(!is.null(status) && status != 0) stop(paste(c("bcftools query failed:", err), collapse="\n"))

    cat(sprintf("Finished: read lines=%d, parsed lines=%d, matched variants=%d/%d\n", nread, nparse, nhit, length(variants)))
    dose
}
