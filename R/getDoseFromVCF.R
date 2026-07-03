#' Read genotype dosages from an indexed VCF/BCF file
#'
#' `getDoseFromVCF()` reads genotype dosages from an indexed VCF or BCF file using
#' `bcftools query`. Either `variants` or `region` must be supplied.
#'
#' If `variants` is supplied, variants are specified as `CHR:POS:REF:ALT`. The
#' function first extracts records overlapping the requested positions and then
#' keeps only records whose chromosome, position, reference allele, and alternate
#' allele exactly match the requested variants. Missing variants are returned as
#' rows filled with `NA`.
#'
#' If `region` is supplied, all variants in the region are read and returned.
#'
#' @param vcf Character scalar. Path to a bgzip-compressed and indexed VCF file
#'   (`.vcf.gz`) or an indexed BCF file (`.bcf`).
#' @param variants Optional character vector. Variants to extract, specified as
#'   `CHR:POS:REF:ALT`, for example `"chr1:12345:A:G"`.
#' @param region Optional character vector. Genomic region(s) to extract, for
#'   example `"chr1:100000-200000"`. Multiple regions are collapsed with commas
#'   and passed to `bcftools query -r`.
#' @param samples Optional character vector. Sample IDs to extract. If `NULL`,
#'   all samples in the VCF/BCF are extracted, unless `sample_file` is supplied.
#' @param sample_file Optional character scalar. Path to a file containing sample
#'   IDs to extract, one sample per line. If supplied, `samples` must be `NULL`.
#' @param BCFTOOLS Character scalar. Path to the `bcftools` executable, or the
#'   command name if `bcftools` is available in the system `PATH`.
#' @param field Character scalar. FORMAT field to read. `"DS"` reads imputed
#'   dosages. `"GT"` converts genotypes to alternate-allele dosages.
#' @param tmpdir Character scalar. Parent directory for temporary files. A unique
#'   subdirectory is created for each function call, so parallel calls do not
#'   share temporary files.
#' @param nline Integer. Number of `bcftools query` output lines read from the
#'   pipe at a time. Smaller values reduce peak memory usage, while larger values
#'   may improve speed.
#' @param progress Logical. If `TRUE`, print progress information while reading
#'   and parsing records.
#'
#' @return A numeric dosage matrix with variants in rows and samples in columns.
#'   In `variants` mode, row names are the requested `CHR:POS:REF:ALT` IDs. In
#'   `region` mode, row names are the observed `CHR:POS:REF:ALT` IDs. Column names
#'   are sample IDs. Variant metadata are stored in the `"variant_info"` attribute.
#'
#' @details
#' This function requires `bcftools` and an index for the input VCF/BCF file.
#' For `.vcf.gz` files, a tabix index such as `.tbi` or `.csi` is expected.
#' For `.bcf` files, a `.csi` index is expected.
#'
#' @examples
#' \dontrun{
#' variants <- c("chr1:12345:A:G", "chr1:67890:C:T")
#'
#' G1 <- getDoseFromVCF(
#'   vcf = "imputed.vcf.gz",
#'   variants = variants,
#'   BCFTOOLS = "bcftools"
#' )
#'
#' G2 <- getDoseFromVCF(
#'   vcf = "imputed.vcf.gz",
#'   region = "chr1:100000-200000",
#'   samples = c("ID001", "ID002"),
#'   BCFTOOLS = "bcftools"
#' )
#'
#' attr(G1, "variant_info")
#' attr(G2, "variant_info")
#' }
#'
#' @export
getDoseFromVCF = function(vcf, variants=NULL, region=NULL, samples=NULL, sample_file=NULL, BCFTOOLS="bcftools", field=c("DS","GT"), tmpdir=tempdir(), nline=20, progress=TRUE){
    field = match.arg(field)
    if(is.null(variants) == is.null(region)) stop("Specify exactly one of 'variants' or 'region'.")
    if(!is.null(samples) && !is.null(sample_file)) stop("Specify only one of 'samples' or 'sample_file'.")
    mode = if(!is.null(variants)) "variants" else "region"

    wdir = tempfile(pattern=paste0("getDoseFromVCF_", Sys.getpid(), "_"), tmpdir=tmpdir)
    dir.create(wdir, recursive=TRUE, showWarnings=FALSE)
    if(!dir.exists(wdir)) stop("Failed to create temporary directory: ", wdir)
    closed = TRUE
    on.exit({
        if(exists("con", inherits=FALSE) && !closed) try(close(con), silent=TRUE)
        unlink(wdir, recursive=TRUE, force=TRUE)
    }, add=TRUE)

    if(!is.null(sample_file)){
        sample_path = sample_file
    }else if(!is.null(samples)){
        sample_path = file.path(wdir, "samples.txt")
        writeLines(as.character(samples), sample_path)
    }else{
        sample_path = NULL
    }

    eh = file.path(wdir, "header.err")
    hargs = if(is.null(sample_path)) c("view","-h",vcf) else c("view","-h","-S",sample_path,vcf)
    header = system2(BCFTOOLS, hargs, stdout=TRUE, stderr=eh)
    if(!is.null(attr(header, "status"))) stop(paste(c("bcftools view -h failed:", readLines(eh, warn=FALSE)), collapse="\n"))

    tag = paste0("##FORMAT=<ID=", field, ",")
    if(!any(grepl(tag, header, fixed=TRUE))) stop("FORMAT/", field, " is not defined in the VCF/BCF header.")

    hline = grep("^#CHROM", header, value=TRUE)
    if(length(hline)==0) stop("Could not find #CHROM header line.")
    hsp = strsplit(hline[length(hline)], "\t", fixed=TRUE)[[1]]
    sample_ids = if(length(hsp) > 9) hsp[-(1:9)] else character(0)
    if(length(sample_ids)==0) stop("No samples were found or selected.")
    ns = length(sample_ids)

    chr_rank = function(x){
        y = sub("^chr", "", x, ignore.case=TRUE)
        z = match(y, c(as.character(1:22), "X", "Y", "XY", "MT", "M"))
        z[is.na(z)] = 1000 + match(y[is.na(z)], unique(y[is.na(z)]))
        z
    }

    if(mode=="variants"){
        variants = as.character(variants)
        sp = strsplit(variants, ":", fixed=TRUE)
        if(any(lengths(sp)!=4)) stop("'variants' must be specified as CHR:POS:REF:ALT.")
        vdat = data.frame(CHR=sapply(sp, `[`, 1), POS=as.integer(sapply(sp, `[`, 2)), REF=sapply(sp, `[`, 3), ALT=sapply(sp, `[`, 4), VID=variants, stringsAsFactors=FALSE)
        if(any(is.na(vdat$POS))) stop("POS in 'variants' must be integer.")
        if(anyDuplicated(variants)) warning("Duplicated variants were supplied; only the first matching row is filled.")

        reg = unique(vdat[,c("CHR","POS")])
        reg = reg[order(chr_rank(reg$CHR), reg$POS),,drop=FALSE]
        bed = file.path(wdir, "target.bed")
        write.table(data.frame(reg$CHR, reg$POS-1L, reg$POS), bed, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

        dose = matrix(NA_real_, nrow=length(variants), ncol=ns, dimnames=list(variants, sample_ids))
        idx = seq_along(variants); names(idx) = variants
        seen = rep(FALSE, length(variants))
        query_arg = paste("-R", shQuote(bed))
        total = length(variants)
    }else{
        region = as.character(region)
        if(length(region)==0 || any(!nzchar(region))) stop("'region' must be a non-empty character vector.")
        dose = NULL
        ids = list()
        mats = list()
        nblk = 0L
        query_arg = paste("-r", shQuote(paste(region, collapse=",")))
        total = NA_integer_
    }

    sarg = if(is.null(sample_path)) "" else paste("-S", shQuote(sample_path))
    fmt = if(field=="DS") "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DS]\\n" else "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    errfile = file.path(wdir, "bcftools.err")
    cmd = paste(shQuote(BCFTOOLS), "query", query_arg, sarg, "-f", shQuote(fmt), shQuote(vcf), "2>", shQuote(errfile))

    bar = function(x, n, width=30){
        p = if(is.na(n) || n <= 0) 0 else min(1, x/n)
        k = floor(width*p)
        paste0("[", paste(rep("=", k), collapse=""), paste(rep(" ", width-k), collapse=""), "]")
    }
    draw = function(nread, nparse, nhit){
        if(!progress) return(invisible())
        if(mode=="variants"){
            cat(sprintf("\rhit %s %d/%d | read lines %d | parsed lines %d", bar(nhit, total), nhit, total, nread, nparse))
        }else{
            cat(sprintf("\rvariants read: %d | parsed lines: %d", nread, nparse))
        }
        flush.console()
    }
    gt2dose = function(x){
        x = gsub("\\|", "/", x)
        y = rep(NA_real_, length(x))
        ok2 = !is.na(x) & grepl("^[0-9]+/[0-9]+$", x)
        if(any(ok2)) y[ok2] = vapply(strsplit(x[ok2], "/", fixed=TRUE), function(g) sum(as.integer(g)==1L), numeric(1))
        ok1 = !is.na(x) & grepl("^[0-9]+$", x)
        if(any(ok1)) y[ok1] = as.numeric(as.integer(x[ok1])==1L)
        y
    }

    con = pipe(cmd, open="r")
    closed = FALSE
    nread = 0L
    nparse = 0L
    nhit = 0L
    draw(0L, 0L, 0L)

    repeat{
        lines = readLines(con, n=nline, warn=FALSE)
        if(length(lines)==0) break
        nread = nread + length(lines)
        z = strsplit(lines, "\t", fixed=TRUE)

        if(mode=="region"){
            id = character(length(z))
            mat = matrix(NA_real_, nrow=length(z), ncol=ns)
        }

        for(j in seq_along(z)){
            a = z[[j]]
            if(length(a)!=(4L+ns)) stop("Unexpected number of columns in bcftools output: got ", length(a), ", expected ", 4L+ns, ".")
            vid = paste(a[1], a[2], a[3], a[4], sep=":")
            x = a[-(1:4)]
            if(field=="DS"){
                x[x=="."] = NA_character_
                val = as.numeric(x)
            }else{
                val = gt2dose(x)
            }

            if(mode=="variants"){
                i = unname(idx[vid])
                if(!is.na(i)){
                    dose[i,] = val
                    if(!seen[i]){
                        seen[i] = TRUE
                        nhit = nhit + 1L
                    }
                }
            }else{
                id[j] = vid
                mat[j,] = val
            }
            nparse = nparse + 1L
        }

        if(mode=="region"){
            nblk = nblk + 1L
            ids[[nblk]] = id
            mats[[nblk]] = mat
        }

        draw(nread, nparse, nhit)
    }

    status = close(con)
    closed = TRUE
    if(progress) cat("\n")

    err = readLines(errfile, warn=FALSE)
    bad = !is.null(status) && length(status)>0 && !is.na(status) && status!=0
    if(bad) stop(paste(c("bcftools query failed:", err), collapse="\n"))

    if(mode=="region"){
        if(length(mats)==0){
            dose = matrix(numeric(0), nrow=0, ncol=ns)
            colnames(dose) = sample_ids
            attr(dose, "variant_info") = data.frame(VID=character(0), CHR=character(0), POS=integer(0), REF=character(0), ALT=character(0), stringsAsFactors=FALSE)
            return(dose)
        }
        dose = do.call(rbind, mats)
        rownames(dose) = unlist(ids, use.names=FALSE)
        colnames(dose) = sample_ids
        sp2 = strsplit(rownames(dose), ":", fixed=TRUE)
        attr(dose, "variant_info") = data.frame(VID=rownames(dose), CHR=sapply(sp2, `[`, 1), POS=as.integer(sapply(sp2, `[`, 2)), REF=sapply(sp2, `[`, 3), ALT=sapply(sp2, `[`, 4), stringsAsFactors=FALSE)
        if(progress) cat(sprintf("Finished: read variants=%d\n", nrow(dose)))
    }else{
        attr(dose, "variant_info") = data.frame(VID=variants, CHR=vdat$CHR, POS=vdat$POS, REF=vdat$REF, ALT=vdat$ALT, found=seen, stringsAsFactors=FALSE)
        attr(dose, "missing_variants") = variants[!seen]
        if(progress) cat(sprintf("Finished: read lines=%d, parsed lines=%d, matched variants=%d/%d\n", nread, nparse, nhit, length(variants)))
    }

    dose
}
