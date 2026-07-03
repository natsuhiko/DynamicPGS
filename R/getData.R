#' Print a DynamicPGS object
#'
#' Prints a short summary of a `DynamicPGS` object, including the number of
#' observations, number of individuals, maximum family size, support of the
#' continuous index, and covariate structure.
#'
#' @param x A `DynamicPGS` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.DynamicPGS <- function(x, ...) {
  cat("DynamicPGS object\n")
  cat("Number of samples: ", x$N, "\n", sep="")
  cat("Number of individuals: ", x$Nd, "\n", sep="")
  cat("Max family size: ", x$max_family_size, "\n", sep="")
  cat("Support: ", min(x$x), " to ", max(x$x), "\n", sep="")
  cat("Covariates:\n")
  print(data.frame(degrees_of_freedom=x$nh))
  if(!is.null(x$rho)){cat("Step 1 optimization completed.\n")}
  if(!is.null(x$delta2d)){cat("Step 2 optimization completed.\n")}
  if(!is.null(x$Beta)){cat("Ready for dynamic PGS computation\n")}
  invisible(x)
}


#' Prepare longitudinal phenotype data for DynamicPGS analysis
#'
#' `getData()` reads longitudinal phenotype data, optional covariates, and optional
#' KING relatedness results, and constructs a `DynamicPGS` analysis object.
#'
#' The input phenotype table must contain at least three columns: `IID`, `x`, and
#' `y`, where `IID` is the individual ID, `x` is the continuous index such as age
#' or time, and `y` is the phenotype. Rows with missing `IID`, `x`, `y`, or
#' covariates are removed. If a KING file is provided, closely related individuals
#' are grouped into family blocks and an additive relationship matrix is converted
#' into a block-wise Cholesky representation used in later mixed-model steps.
#'
#' @param Data A data.frame or path to a tab-delimited file containing at least
#'   `IID`, `x`, and `y`.
#' @param Covariates Optional data.frame or path to a covariate file. Numeric
#'   covariates are standardised internally; character or factor covariates are
#'   expanded into dummy variables.
#' @param king A data.frame or an optional path to a KING pairwise relatedness result file.
#'   The file is expected to contain columns `ID1`, `ID2`, and `Kinship`.
#' @param inducing_points Optional numeric vector of inducing points on the `x`
#'   axis. If `NULL`, approximately 11 equally spaced points between the minimum
#'   and maximum observed `x` are used.
#' @param forced Logical. If `TRUE`, suppresses interactive prompts about merging
#'   adjacent numeric covariate groups.
#'
#' @return A list of class `DynamicPGS` containing cleaned phenotype data,
#'   design matrices, inducing points, covariance initial values, individual and
#'   family mappings, and relatedness information.
#'
#' @examples
#' \dontrun{
#' adata <- getData(
#'   Data = "phenotype.tsv",
#'   Covariates = "covariates.tsv",
#'   king = "king.kin0",
#'   inducing_points = seq(0, 60, by = 6)
#' )
#' }
#'
#' @export
getData = function(Data="/path/to/your/data_body.tsv.gz", Covariates=NULL, king=NULL, inducing_points=NULL, forced=T, Verbose=F){
    
    if(is.data.frame(Data)){
        body = Data
    }else{
        body = read.table(Data, header=T, sep="\t")
    }
    iid0=body$IID
    N0=nrow(body)
    
    narm = is.na(body$y)+is.na(body$x)
    X0=NULL
    if(!is.null(Covariates)){
        if(is.data.frame(Covariates)){
            X0 = Covariates
        }else{
            X0 = read.table(Covariates, header=T)
            if(nrow(X0)!=nrow(body)){return("N rows inconsistent...")}
            narm = narm + rowSums(is.na(X0))
        }
    }
    y = body$y[narm==0]
    x = body$x[narm==0]
    iid = body$IID[narm==0]
    if(!is.null(X0)){ X0 = X0[narm==0,] }
    uid = unique(iid)
    if(is.null(inducing_points)){
        ta = seq(min(x),max(x),diff(range(x))/10)
    }else{
        ta = inducing_points
    }
    
    # King
    if(!is.null(king)){
        kingres = make_arm(king,uid)
        lmatres = makeL(kingres$A, kingres$family_df)
        Lmat = cbind(lmatres[[2]][,1:2],lmatres[[1]])
        max_family_size = max(lmatres$max_family_size)
    }else{
        Lmat = data.frame(FID=uid,IID=uid,chol=rep(1.0,length(uid)))
        max_family_size = 1
    }
    
    ord = order(match(iid,Lmat$IID))
    
    y = y[ord]
    x = x[ord]
    X0 = X0[ord,]
    iid = iid[ord] # Sort by order of appearance in Lmat$IID
    mapid = match(iid, unique(iid))

    N = length(y)
    M = length(ta)
    tid = table(mapid)
    Nd = length(tid)
    fid = Lmat$FID[match(iid,Lmat$IID)]
    Nf = length(unique(fid))
    
    if(Verbose){
    cat("Nd: ");cat(Nd);cat("\n")
    cat("Nf: ");cat(Nd);cat("\n")
    cat("N: ");cat(N);cat("\n")
    cat("M: ");cat(M);cat("\n")
    }
    
    # covariates
    nh = 1
    X = rep(1,N)
    isnum = 0
    delta2 = mean(y)^2
    for(i in seq(ncol(X0))){
        if(is.character(X0[,i]) || is.factor(X0[,i])){
            isnum = c(isnum, 0)
            Z1 = array(0,c(N,length(table(X0[,i]))))
            Z1[!is.na(X0[,i]),] = model.matrix(~0+X0[,i])
            delta2=c(delta2,var(coef(lm(y~Z1-1))))
        }else{
            isnum = c(isnum, 1)
            Z1 = matrix(scale(as.numeric(X0[,i])),N)
            Z1[is.na(Z1)]=0
            delta2=c(delta2,coef(lm(y~Z1))[2]^2)
        }
        X  = cbind(X, Z1)
        nh = c(nh, ncol(Z1))
    }
    names(nh)=c("Intercept", names(X0))
    isnum = cumsum(1-isnum)*isnum
    if(sum(isnum)){for(i in seq(max(isnum))){
        if(sum(isnum==i)>1){
            if(!forced){
                mflag = readline(paste("Do you want to merge ", paste(names(nh)[isnum==i],collapse=","), " factors? [N/y]: ",sep=""))
            }else{
                mflag="N"
            }
            if(sum(mflag%in%c("y","Y"))){
                nh    = Collapse(nh, isnum==i)
                isnum = Collapse(isnum, isnum==i)
            }
        }
    }}
    P = length(nh)
    Q = ncol(X)
    if(Verbose)print(nh)
    adata = list(y=y, x=x, ta=ta, X=X, nh=nh, iid=iid, fid=fid, P=P, Q=Q, N=N, M=M, Nf=Nf, Nd=Nd, delta2=delta2, Lmat=Lmat, max_family_size=max_family_size, support_x=range(x))
    class(adata)="DynamicPGS"
    adata
}



makeL <- function(
    A,
    family_df,
    iid_col = "IID",
    fid_col = "FID",
    jitter = 1e-8,
    return_A_sorted = TRUE
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install the Matrix package.")
  }
  
  if (is.null(rownames(A)) || is.null(colnames(A))) {
    stop("A must have rownames and colnames.")
  }
  
  if (!identical(rownames(A), colnames(A))) {
    stop("rownames(A) and colnames(A) must be identical.")
  }
  
  if (!all(c(iid_col, fid_col) %in% colnames(family_df))) {
    stop("family_df must contain IID and FID columns.")
  }
  
  family_df[[iid_col]] <- as.character(family_df[[iid_col]])
  family_df[[fid_col]] <- as.character(family_df[[fid_col]])
  
  ids_in_A <- rownames(A)
  
  if (!all(ids_in_A %in% family_df[[iid_col]])) {
    stop("Some IDs in A are missing from family_df.")
  }
  
  if (!all(family_df[[iid_col]] %in% ids_in_A)) {
    stop("Some IDs in family_df are missing from A.")
  }
  
  #----------------------------
  # 1. Sort by FID
  #----------------------------
  family_df <- family_df[match(ids_in_A, family_df[[iid_col]]), ]
  
  ord <- order(family_df[[fid_col]], family_df[[iid_col]])
  
  family_df_sorted <- family_df[ord, ]
  sorted_ids <- family_df_sorted[[iid_col]]
  
  A_sorted <- A[sorted_ids, sorted_ids]
  
  N <- length(sorted_ids)
  
  family_sizes <- table(family_df_sorted[[fid_col]])
  max_family_size <- max(as.integer(family_sizes))
  
  #----------------------------
  # 2. Output matrix
  #----------------------------
  L_left <- matrix(
    0,
    nrow = N,
    ncol = max_family_size
  )
  
  rownames(L_left) <- sorted_ids
  colnames(L_left) <- paste0("chol", seq_len(max_family_size))
  
  #----------------------------
  # 3. Cholesky by family block
  #----------------------------
  row_start <- 1
  
  block_info <- data.frame(
    FID = names(family_sizes),
    row_start = NA_integer_,
    row_end = NA_integer_,
    family_size = as.integer(family_sizes),
    stringsAsFactors = FALSE
  )
  
  for (b in seq_along(family_sizes)) {
    fid <- names(family_sizes)[b]
    m <- as.integer(family_sizes[b])
    
    rows <- row_start:(row_start + m - 1)
    
    A_block <- as.matrix(A_sorted[rows, rows])
    
    # Ensure exact symmetry
    A_block <- (A_block + t(A_block)) / 2
    
    # Cholesky in R returns upper triangular R:
    # A = t(R) %*% R
    # Therefore L = t(R) gives A = L %*% t(L)
    R <- tryCatch(
      chol(A_block),
      error = function(e) {
        chol(A_block + diag(jitter, m))
      }
    )
    
    L <- t(R)
    
    # left-align the family block
    L_left[rows, seq_len(m)] <- L
    
    block_info$row_start[b] <- min(rows)
    block_info$row_end[b] <- max(rows)
    
    row_start <- row_start + m
  }
  
  out <- list(
    L_left = L_left,
    family_df_sorted = family_df_sorted,
    block_info = block_info,
    max_family_size = max_family_size
  )
  
  if (return_A_sorted) {
    out$A_sorted <- A_sorted
  }
  
  return(out)
}




make_arm <- function(
    king,
    given_ids,
    id1_col = "ID1",
    id2_col = "ID2",
    kinship_col = "Kinship",
    kin_cutoff = 0.0442,
    fid_prefix = "FAM",
    unrelated_prefix = "UNR"
) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install the Matrix package.")
  }
  
  #----------------------------
  # 1. Read KING result
  #----------------------------
  if(!is.data.frame(king)){
      king <- read.table(
        king,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
  }
  
  
  required_cols <- c(id1_col, id2_col, kinship_col)
  if (!all(required_cols %in% colnames(king))) {
    stop("KING file does not contain required columns.")
  }
  
  ID1 <- as.character(king[[id1_col]])
  ID2 <- as.character(king[[id2_col]])
  Kinship <- as.numeric(king[[kinship_col]])
  
  ok <- !is.na(ID1) & !is.na(ID2) & !is.na(Kinship)
  ID1 <- ID1[ok]
  ID2 <- ID2[ok]
  Kinship <- Kinship[ok]
  
  #----------------------------
  # 2. Prepare given IDs
  #----------------------------
  given_ids <- unique(as.character(given_ids))
  given_ids <- sort(given_ids)
  
  n <- length(given_ids)
  
  id_index <- seq_len(n)
  names(id_index) <- given_ids
  
  #----------------------------
  # 3. Keep KING pairs within given IDs
  #----------------------------
  keep <- (
    Kinship > kin_cutoff &
      ID1 %in% given_ids &
      ID2 %in% given_ids
  )
  
  ID1_keep <- ID1[keep]
  ID2_keep <- ID2[keep]
  Kinship_keep <- Kinship[keep]
  
  #----------------------------
  # 4. Make sparse additive relationship matrix
  #
  # KING Kinship = phi
  # Additive relationship A_ij = 2 * phi_ij
  # Diagonal A_ii = 1
  #----------------------------
  i_diag <- seq_len(n)
  j_diag <- seq_len(n)
  x_diag <- rep(1, n)
  
  if (length(Kinship_keep) > 0) {
    i_off <- id_index[ID1_keep]
    j_off <- id_index[ID2_keep]
    x_off <- 2 * Kinship_keep
    
    i <- c(i_diag, i_off, j_off)
    j <- c(j_diag, j_off, i_off)
    x <- c(x_diag, x_off, x_off)
  } else {
    i <- i_diag
    j <- j_diag
    x <- x_diag
  }
  
  A <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(n, n),
    dimnames = list(given_ids, given_ids)
  )
  
  #----------------------------
  # 5. Assign FID using union-find
  #----------------------------
  parent <- seq_len(n)
  
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  
  union_sets <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (ra != rb) {
      parent[rb] <<- ra
    }
  }
  
  if (length(Kinship_keep) > 0) {
    for (k in seq_along(Kinship_keep)) {
      i <- id_index[ID1_keep[k]]
      j <- id_index[ID2_keep[k]]
      union_sets(i, j)
    }
  }
  
  roots <- integer(n)
  for (i in seq_len(n)) {
    roots[i] <- find_root(i)
  }
  
  component <- match(roots, unique(roots))
  family_size <- tabulate(component)[component]
  
  #----------------------------
  # 6. Create readable FID
  #----------------------------
  FID <- character(n)
  
  fam_components <- sort(unique(component[family_size > 1]))
  fam_number <- seq_along(fam_components)
  names(fam_number) <- fam_components
  
  unrelated_counter <- 1
  
  for (i in seq_len(n)) {
    if (family_size[i] > 1) {
      FID[i] <- paste0(
        fid_prefix,
        sprintf("%07d", fam_number[as.character(component[i])])
      )
    } else {
      FID[i] <- paste0(
        unrelated_prefix,
        sprintf("%07d", unrelated_counter)
      )
      unrelated_counter <- unrelated_counter + 1
    }
  }
  
  family_df <- data.frame(
    FID = FID,
    IID = given_ids,
    family_component = component,
    family_size = family_size,
    stringsAsFactors = FALSE
  )
  
  #----------------------------
  # 7. Return
  #----------------------------
  return(list(
    A = A,
    family_df = family_df,
    king_used = data.frame(
      ID1 = ID1_keep,
      ID2 = ID2_keep,
      Kinship = Kinship_keep,
      AdditiveRelationship = 2 * Kinship_keep,
      stringsAsFactors = FALSE
    )
  ))
}
