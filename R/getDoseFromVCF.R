#' Read genotype dosages from an indexed VCF file
#'
#' `getDoseFromVCF()` extracts variants from a bgzip-compressed and tabix-indexed
#' VCF file and returns a dosage matrix. If the VCF contains a `DS` field, dosage
#' values are read from `DS`. Otherwise, genotypes in the `GT` field are converted
#' to alternate-allele dosages.
#'
#' @param vcf Path to a bgzip-compressed VCF file indexed by tabix.
#' @param region Character vector of genomic regions passed to tabix, for example
#'   `"chr1:100000-200000"`. At least one region must be supplied.
#' @param tabix Path or command name for the tabix executable.
#' @param chunk_size Integer. Number of VCF lines read per chunk.
#'
#' @return A numeric dosage matrix with variants in rows and samples in columns.
#'   Row names are constructed as `CHROM:POS:REF:ALT`. Variant metadata are stored
#'   in the `"variant_info"` attribute.
#'
#' @examples
#' \dontrun{
#' Gall <- getDoseFromVCF(
#'   vcf = "imputed.vcf.gz",
#'   region = "chr1:100000-200000"
#' )
#' attr(Gall, "variant_info")
#' }
#'
#' @export
getDoseFromVCF <- function(vcf, region, tabix="tabix", chunk_size=1000){
  if(length(region)<1) stop("region must contain at least one region.")
  
  cmd <- paste(c(shQuote(tabix), "-h", shQuote(vcf), shQuote(region)), collapse=" ")
  con <- pipe(cmd, open="r")
  on.exit(close(con))
  
  header <- NULL
  sample_ids <- NULL
  D_list <- list()
  V_list <- list()
  k <- 0L
  
  get_field <- function(x, idx){
    y <- strsplit(x, ":", fixed=TRUE)
    vapply(y, function(z){
      if(length(z) < idx) NA_character_ else z[[idx]]
    }, character(1))
  }
  
  gt_to_dosage <- function(gt){
    vapply(gt, function(g){
      if(is.na(g) || g %in% c(".", "./.", ".|.", "")) return(NA_real_)
      a <- strsplit(g, "[/|]")[[1]]
      if(any(a==".")) return(NA_real_)
      sum(as.integer(a) > 0)
    }, numeric(1))
  }
  
  process_chunk <- function(lines, header, sample_ids){
    f <- strsplit(lines, "\t", fixed=TRUE)
    n <- length(f)
    ns <- length(sample_ids)
    
    chrom <- vapply(f, `[`, character(1), 1)
    pos   <- vapply(f, `[`, character(1), 2)
    id    <- vapply(f, `[`, character(1), 3)
    ref   <- vapply(f, `[`, character(1), 4)
    alt   <- vapply(f, `[`, character(1), 5)
    qual  <- vapply(f, `[`, character(1), 6)
    filt  <- vapply(f, `[`, character(1), 7)
    info  <- vapply(f, `[`, character(1), 8)
    fmt   <- vapply(f, `[`, character(1), 9)
    
    S <- matrix(
      unlist(lapply(f, function(z) z[10:length(z)]), use.names=FALSE),
      nrow=n, ncol=ns, byrow=TRUE
    )
    
    D <- matrix(NA_real_, nrow=n, ncol=ns)
    
    fmt_groups <- split(seq_len(n), fmt)
    
    for(ii in fmt_groups){
      ff <- strsplit(fmt[ii[1]], ":", fixed=TRUE)[[1]]
      ds_idx <- match("DS", ff)
      gt_idx <- match("GT", ff)
      
      x <- as.vector(t(S[ii,,drop=FALSE]))
      
      if(!is.na(ds_idx)){
        val <- get_field(x, ds_idx)
        val[val %in% c(".", "")] <- NA_character_
        D[ii,] <- matrix(as.numeric(val), nrow=length(ii), byrow=TRUE)
      } else if(!is.na(gt_idx)){
        gt <- get_field(x, gt_idx)
        D[ii,] <- matrix(gt_to_dosage(gt), nrow=length(ii), byrow=TRUE)
      }
    }
    
    rownames(D) <- paste(chrom, pos, ref, alt, sep=":")
    colnames(D) <- sample_ids
    
    vinfo <- data.frame(
      CHROM=chrom, POS=pos, ID=id, REF=ref, ALT=alt,
      QUAL=qual, FILTER=filt, INFO=info,
      stringsAsFactors=FALSE,
      check.names=FALSE
    )
    
    list(D=D, vinfo=vinfo)
  }
  
  repeat{
    z <- readLines(con, n=chunk_size, warn=FALSE)
    if(length(z)==0) break
    
    if(is.null(header)){
      h <- grep("^#CHROM", z)
      if(length(h)>0){
        header <- strsplit(z[h[1]], "\t", fixed=TRUE)[[1]]
        sample_ids <- header[10:length(header)]
      }
    }
    
    body <- z[!grepl("^#", z)]
    if(length(body)==0) next
    
    if(is.null(header)) stop("VCF header line '#CHROM' was not found.")
    
    out <- process_chunk(body, header, sample_ids)
    k <- k + 1L
    D_list[[k]] <- out$D
    V_list[[k]] <- out$vinfo
  }
  
  if(k==0L) stop("No variants returned from tabix.")
  
  D <- do.call(rbind, D_list)
  vinfo <- do.call(rbind, V_list)
  
  keep <- !duplicated(rownames(D))
  D <- D[keep,,drop=FALSE]
  vinfo <- vinfo[keep,,drop=FALSE]
  
  attr(D, "variant_info") <- vinfo
  D
}
