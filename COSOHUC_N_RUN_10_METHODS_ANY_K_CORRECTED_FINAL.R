# ================================================================
# G-RUN K-MEANS CLUSTERING COMPARISON
# ================================================================
#
# METHODS
#
# 1. COSOHUC
# 2. Farthest First
# 3. Canopy
# 4. K-Means++
# 5. Hartigan-Wong
# 6. MacQueen
# 7. Lloyd
# 8. Forgy
# 9. ECKM (Empty Circles based K-means)
# 10. Mini-Batch K-Means++
#
# ================================================================
#
# IMPORTANT
# ---------------------------------------------------------------
# COSOHUC:
#   - Uses the original deterministic initialization.
#   - Uses the original kmeans() call.
#   - Is calculated ONLY ONCE.
#   - The same result is copied to every run.
#
# STOCHASTIC METHODS:
#   - Use their original random-center implementation.
#   - Use a different fixed seed for every run.
#
# ECKM:
#   - Faithful R translation of the supplied published/official procedure.
#   - Original paper implementation is Python; R/deldir may differ slightly.
#   - Deterministic for a fixed dataset and k.
#   - Calculated ONLY ONCE and copied to every run.
#
# MINI-BATCH K-MEANS++:
#   - Implements the uploaded paper-based Mini-Batch K-Means++ initialization.
#   - Stochastic: first center is randomly selected.
#   - A different fixed seed is used for every run.
#   - The selected centers are passed to standard Lloyd K-means.
#
# ================================================================


# ================================================================
# 1. REQUIRED PACKAGES
# ================================================================

required_packages <- c(
  "BiocGenerics",
  "deldir",
  "clusterSim",
  "fossil",
  "factoextra",
  "NbClust",
  "clusterCrit",
  "cluster"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
}

library(BiocGenerics)
library(clusterSim)
library(fossil)
library(factoextra)
library(NbClust)
library(clusterCrit)
library(cluster)


# ================================================================
# 2. RAND INDEX
# ================================================================

calculate_rand_index <- function(
    true_labels,
    predicted_labels
) {
  
  if (
    length(true_labels) !=
    length(predicted_labels)
  ) {
    
    stop(
      "True labels and predicted labels must have the same length."
    )
    
  }
  
  true_labels <- as.factor(
    true_labels
  )
  
  predicted_labels <- as.factor(
    predicted_labels
  )
  
  n <- length(
    true_labels
  )
  
  if (
    n < 2
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  total_pairs <- choose(
    n,
    2
  )
  
  true_same <- outer(
    true_labels,
    true_labels,
    FUN = "=="
  )
  
  predicted_same <- outer(
    predicted_labels,
    predicted_labels,
    FUN = "=="
  )
  
  agreement <- (
    true_same ==
      predicted_same
  )
  
  agreement_pairs <- agreement[
    upper.tri(
      agreement
    )
  ]
  
  RI <- sum(
    agreement_pairs
  ) / total_pairs
  
  return(
    as.numeric(RI)
  )
}


# ================================================================
# 3. BROWSE DATASET
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "SELECT DATASET FILE\n"
)

cat(
  "============================================================\n"
)

file_path <- file.choose()

if (
  is.null(file_path) ||
  length(file_path) == 0
) {
  
  stop(
    "No dataset file selected."
  )
  
}

file_path <- normalizePath(
  file_path,
  winslash = "/",
  mustWork = TRUE
)

cat(
  "\nSelected dataset:\n",
  file_path,
  "\n"
)


# ================================================================
# 4. READ DATASET
# ================================================================

p <- read.csv(
  file = file_path,
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ================================================================
# 5. DATASET CHECK
# ================================================================

if (
  nrow(p) < 2
) {
  
  stop(
    "Dataset must contain at least 2 observations."
  )
  
}

if (
  ncol(p) < 1
) {
  
  stop(
    "Dataset contains no columns."
  )
  
}

cat(
  "\n============================================================\n"
)

cat(
  "DATASET INFORMATION\n"
)

cat(
  "============================================================\n"
)

cat(
  "Rows    :",
  nrow(p),
  "\n"
)

cat(
  "Columns :",
  ncol(p),
  "\n"
)

cat(
  "\nDataset columns:\n"
)

for (
  i in seq_len(ncol(p))
) {
  
  cat(
    i,
    ":",
    names(p)[i],
    "\n"
  )
  
}


# ================================================================
# 6. RAND INDEX Y/N
# ================================================================

repeat {
  
  ri_answer <- toupper(
    trimws(
      readline(
        prompt =
          "\nDo you want to calculate Rand Index? (Y/N): "
      )
    )
  )
  
  if (
    ri_answer %in%
    c("Y", "N")
  ) {
    
    break
    
  }
  
  cat(
    "\nINVALID INPUT. Please enter Y or N.\n"
  )
  
}

calculate_RI <- (
  ri_answer == "Y"
)


# ================================================================
# 7. GROUND-TRUTH VARIABLES
# ================================================================

true_labels <- NULL

truth_column <- NULL

truth_name <- "Not used"

number_of_true_classes <- NA_integer_


# ================================================================
# 8. GROUND-TRUTH COLUMN
# ================================================================

if (
  calculate_RI
) {
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "GROUND-TRUTH COLUMN SELECTION\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  for (
    i in seq_len(ncol(p))
  ) {
    
    cat(
      i,
      ":",
      names(p)[i],
      "\n"
    )
    
  }
  
  repeat {
    
    truth_input <- trimws(
      readline(
        prompt =
          "\nEnter the ground-truth column number: "
      )
    )
    
    if (
      grepl(
        "^[0-9]+$",
        truth_input
      )
    ) {
      
      truth_column <- as.integer(
        truth_input
      )
      
    } else {
      
      truth_column <- NA_integer_
      
    }
    
    if (
      !is.na(truth_column) &&
      truth_column >= 1 &&
      truth_column <= ncol(p)
    ) {
      
      break
      
    }
    
    cat(
      "\nINVALID GROUND-TRUTH COLUMN.\n"
    )
    
  }
  
  true_labels <- p[
    ,
    truth_column
  ]
  
  truth_name <- names(p)[
    truth_column
  ]
  
  number_of_true_classes <- length(
    unique(
      true_labels
    )
  )
  
  cat(
    "\nSelected ground-truth column:",
    truth_name,
    "\n"
  )
  
  cat(
    "Number of true classes =",
    number_of_true_classes,
    "\n"
  )
  
  # --------------------------------------------------------------
  # Remove ground-truth column from clustering data
  # --------------------------------------------------------------
  
  d <- p[
    ,
    -truth_column,
    drop = FALSE
  ]
  
} else {
  
  d <- p
  
}


# ================================================================
# 9. CHECK NUMERIC CLUSTERING COLUMNS
# ================================================================

non_numeric <- !sapply(
  d,
  is.numeric
)

if (
  any(non_numeric)
) {
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "NON-NUMERIC COLUMNS FOUND\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  print(
    names(d)[non_numeric]
  )
  
  stop(
    paste0(
      "\nK-means requires numeric clustering attributes.\n",
      "Non-numeric clustering columns were found."
    )
  )
  
}


if (
  ncol(d) < 1
) {
  
  stop(
    "No numeric clustering attributes are available."
  )
  
}


# ================================================================
# 10. NUMERIC MATRIX
# ================================================================

X <- as.matrix(
  d
)

storage.mode(X) <- "numeric"


if (
  anyNA(X)
) {
  
  stop(
    "Dataset contains missing values."
  )
  
}


if (
  any(!is.finite(X))
) {
  
  stop(
    "Dataset contains infinite values."
  )
  
}


# ================================================================
# 11. USER ENTERS k
# ================================================================

repeat {
  
  k_input <- trimws(
    readline(
      prompt =
        paste0(
          "\nEnter number of clusters (k) [2-",
          nrow(X) - 1,
          "]: "
        )
    )
  )
  
  if (
    !grepl(
      "^[0-9]+$",
      k_input
    )
  ) {
    
    cat(
      "\nINVALID k.\n"
    )
    
    next
    
  }
  
  k <- suppressWarnings(
    as.integer(
      k_input
    )
  )
  
  if (
    !is.na(k) &&
    k >= 2 &&
    k < nrow(X)
  ) {
    
    break
    
  }
  
  cat(
    "\nINVALID k.\n"
  )
  
}


# ================================================================
# 12. USER ENTERS MINI-BATCH SIZE (b)
# ================================================================
#
# The Mini-Batch K-Means++ paper defines n <= B*b.
# B is calculated as ceiling(n/b).
#
# ================================================================

repeat {
  b_input <- trimws(
    readline(
      prompt =
        paste0(
          "\nEnter Mini-Batch K-Means++ batch size (b) [1-",
          nrow(X),
          "]: "
        )
    )
  )

  if (!grepl("^[0-9]+$", b_input)) {
    cat("\nINVALID batch size. Enter an integer.\n")
    next
  }

  b <- suppressWarnings(as.integer(b_input))

  if (!is.na(b) && b >= 1 && b <= nrow(X)) {
    break
  }

  cat("\nINVALID batch size.\n")
}

B <- ceiling(nrow(X) / b)


# ================================================================
# 13. USER ENTERS N
# ================================================================

repeat {
  
  N_input <- trimws(
    readline(
      prompt =
        "\nEnter number of runs (N): "
    )
  )
  
  if (
    !grepl(
      "^[0-9]+$",
      N_input
    )
  ) {
    
    cat(
      "\nINVALID N. Enter a positive integer.\n"
    )
    
    next
    
  }
  
  N <- suppressWarnings(
    as.integer(
      N_input
    )
  )
  
  if (
    !is.na(N) &&
    N >= 1
  ) {
    
    break
    
  }
  
  cat(
    "\nINVALID N. N must be >= 1.\n"
  )
  
}


# ================================================================
# 14. TRUE CLASS CHECK
# ================================================================

if (
  calculate_RI
) {
  
  if (number_of_true_classes != k) {

    cat(
      paste0(
        "NOTE: true classes = ",
        number_of_true_classes,
        ", requested k = ",
        k,
        ". Rand Index will still be calculated.
"
      )
    )
  }
  
}


# ================================================================
# 15. DISPLAY EXPERIMENT INFORMATION
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT INFORMATION\n"
)

cat(
  "============================================================\n"
)

cat(
  "Dataset file          :",
  basename(file_path),
  "\n"
)

cat(
  "Observations          :",
  nrow(X),
  "\n"
)

cat(
  "Clustering attributes :",
  ncol(X),
  "\n"
)

cat(
  "Clusters (k)          :",
  k,
  "\n"
)

cat(
  "Number of runs (N)    :",
  N,
  "\n"
)

cat(
  "Mini-batch size (b)   :",
  b,
  "\n"
)

cat(
  "Number of batches (B) :",
  B,
  "\n"
)

cat(
  "Rand Index            :",
  ifelse(
    calculate_RI,
    "YES",
    "NO"
  ),
  "\n"
)

if (
  calculate_RI
) {
  
  cat(
    "Ground-truth column   :",
    truth_name,
    "\n"
  )
  
}


# ================================================================
# 15. COSOHUC INITIALIZATION
# ================================================================
#
# THIS IS THE DETERMINISTIC COSOHUC CODE.
#
# IMPORTANT:
# - No set.seed()
# - No sample()
# - No random center selection
# - No run-dependent operation
#
# ================================================================

COSOHUC_initial_centers <- function(
    X,
    k
) {
  
  gf <- X
  
  cc <- ncol(
    gf
  )
  
  rr <- nrow(
    gf
  )
  
  thre <- rr / k
  
  md <- matrix(
    data = NA_real_,
    ncol = cc,
    nrow = k
  )
  
  colnames(md) <-
    colnames(gf)
  
  
  for (
    j in 1:cc
  ) {
    
    ee <- gf[
      ,
      j
    ]
    
    strth <- 1
    
    endth <- thre
    
    
    for (
      ce in 1:k
    ) {
      
      # ----------------------------------------------------------
      # EXACT ORDER OPERATION
      # ----------------------------------------------------------
      
      ord <- order(
        gf[
          ,
          j
        ],
        decreasing = FALSE
      )[1:rr]
      
      
      # ----------------------------------------------------------
      # EXACT INDEX CALCULATION
      # ----------------------------------------------------------
      
      start_index <- ceiling(
        strth
      )
      
      
      if (
        ce < k
      ) {
        
        end_index <- floor(
          endth
        )
        
      } else {
        
        end_index <- rr
        
      }
      
      
      # ----------------------------------------------------------
      # EXACT VALUE SELECTION
      # ----------------------------------------------------------
      
      selected_values <- ee[
        ord[
          start_index:end_index
        ]
      ]
      
      
      # ----------------------------------------------------------
      # EXACT INITIAL CENTER CALCULATION
      # ----------------------------------------------------------
      
      md[
        ce,
        j
      ] <- mean(
        selected_values,
        na.rm = TRUE
      )
      
      
      # ----------------------------------------------------------
      # UPDATE RANGE
      # ----------------------------------------------------------
      
      strth <- endth
      
      endth <-
        endth + thre
      
    }
    
  }
  
  return(
    md
  )
}


# ================================================================
# 16. CALCULATE COSOHUC CENTERS ONLY ONCE
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "CALCULATING COSOHUC INITIAL CENTERS\n"
)

cat(
  "============================================================\n"
)


COSOHUC_initial <- COSOHUC_initial_centers(
  X,
  k
)


# ================================================================
# 17. COSOHUC K-MEANS
# ================================================================
#
# VERY IMPORTANT:
#
# This is the SAME kmeans() call used in your verified COSOHUC
# code.
#
# There is NO:
#
#       algorithm = "Lloyd"
#
# here.
#
# COSOHUC uses the default R kmeans algorithm.
#
# ================================================================

COSOHUC_result <- kmeans(
  X,
  centers = COSOHUC_initial,
  iter.max = 100,
  nstart = 1
)


COSOHUC_result$initial.centers <-
  COSOHUC_initial


# ================================================================
# 18. DUNN INDEX
# ================================================================

calculate_dunn <- function(
    X,
    cluster
) {
  
  D <- as.matrix(
    dist(X)
  )
  
  clusters <- sort(
    unique(cluster)
  )
  
  if (
    length(clusters) < 2
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  maximum_intra <- 0
  
  for (
    cl in clusters
  ) {
    
    members <- which(
      cluster == cl
    )
    
    if (
      length(members) > 1
    ) {
      
      intra <- D[
        members,
        members,
        drop = FALSE
      ]
      
      maximum_intra <-
        max(
          maximum_intra,
          max(intra)
        )
      
    }
    
  }
  
  minimum_inter <- Inf
  
  if (
    length(clusters) > 1
  ) {
    
    for (
      a in 1:(length(clusters) - 1)
    ) {
      
      for (
        b in (a + 1):length(clusters)
      ) {
        
        members_a <- which(
          cluster == clusters[a]
        )
        
        members_b <- which(
          cluster == clusters[b]
        )
        
        inter <- D[
          members_a,
          members_b,
          drop = FALSE
        ]
        
        minimum_inter <-
          min(
            minimum_inter,
            min(inter)
          )
        
      }
      
    }
    
  }
  
  if (
    maximum_intra <= 0 ||
    !is.finite(minimum_inter)
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  return(
    minimum_inter /
      maximum_intra
  )
}


# ================================================================
# 19. SILHOUETTE
# ================================================================

calculate_silhouette <- function(
    X,
    cluster
) {
  
  tryCatch(
    
    {
      
      sil <- cluster::silhouette(
        as.integer(cluster),
        dist(X)
      )
      
      mean(
        sil[
          ,
          "sil_width"
        ],
        na.rm = TRUE
      )
      
    },
    
    error = function(e) {
      
      NA_real_
      
    }
    
  )
}


# ================================================================
# 20. DAVIES-BOULDIN
# ================================================================

calculate_davies_bouldin <- function(
    X,
    cluster
) {
  
  clusters <- sort(
    unique(cluster)
  )
  
  k_current <- length(
    clusters
  )
  
  if (
    k_current < 2
  ) {
    
    return(
      NA_real_
    )
    
  }
  
  centroids <- matrix(
    NA_real_,
    nrow = k_current,
    ncol = ncol(X)
  )
  
  scatter <- numeric(
    k_current
  )
  
  for (
    i in seq_len(k_current)
  ) {
    
    members <- which(
      cluster == clusters[i]
    )
    
    cluster_data <- X[
      members,
      ,
      drop = FALSE
    ]
    
    centroid <- colMeans(
      cluster_data
    )
    
    centroids[
      i,
    ] <- centroid
    
    differences <- sweep(
      cluster_data,
      2,
      centroid,
      "-"
    )
    
    scatter[i] <-
      mean(
        sqrt(
          rowSums(
            differences^2
          )
        )
      )
    
  }
  
  centroid_distances <- as.matrix(
    dist(
      centroids
    )
  )
  
  db_values <- numeric(
    k_current
  )
  
  for (
    i in seq_len(k_current)
  ) {
    
    ratios <- numeric(
      0
    )
    
    for (
      j in seq_len(k_current)
    ) {
      
      if (
        i != j &&
        centroid_distances[
          i,
          j
        ] > 0
      ) {
        
        ratios <- c(
          ratios,
          (
            scatter[i] +
              scatter[j]
          ) /
            centroid_distances[
              i,
              j
            ]
        )
        
      }
      
    }
    
    if (
      length(ratios) == 0
    ) {
      
      return(
        NA_real_
      )
      
    }
    
    db_values[i] <-
      max(
        ratios
      )
    
  }
  
  return(
    mean(
      db_values
    )
  )
}


# ================================================================
# 21. CALCULATE ALL INDICES
# ================================================================

calculate_all_indices <- function(
    X,
    result
) {
  
  dunn_value <-
    calculate_dunn(
      X,
      result$cluster
    )
  
  silhouette_value <-
    calculate_silhouette(
      X,
      result$cluster
    )
  
  dbi_value <-
    calculate_davies_bouldin(
      X,
      result$cluster
    )
  
  bss_percentage <-
    
    if (!is.null(result$totss) && length(result$totss) == 1 && is.finite(result$totss) && result$totss > 0) {
      
      100 *
        result$betweenss /
        result$totss
      
    } else {
      
      NA_real_
      
    }
  
  return(
    list(
      
      Iteration =
        as.numeric(
          result$iter
        ),
      
      BSS =
        as.numeric(
          result$betweenss
        ),
      
      TWSS =
        as.numeric(
          result$tot.withinss
        ),
      
      TotalSS =
        as.numeric(
          result$totss
        ),
      
      BSS_Percentage =
        as.numeric(
          bss_percentage
        ),
      
      Dunn =
        as.numeric(
          dunn_value
        ),
      
      DBI =
        as.numeric(
          dbi_value
        ),
      
      Silhouette =
        as.numeric(
          silhouette_value
        )
      
    )
  )
}


# ================================================================
# 22. DISPLAY COSOHUC RESULT
# ================================================================

COSOHUC_indices <-
  calculate_all_indices(
    X,
    COSOHUC_result
  )


cat(
  "\n============================================================\n"
)

cat(
  "COSOHUC RESULT\n"
)

cat(
  "============================================================\n"
)


cat(
  "Iteration      =",
  COSOHUC_indices$Iteration,
  "\n"
)

cat(
  "BSS            =",
  COSOHUC_indices$BSS,
  "\n"
)

cat(
  "TWSS           =",
  COSOHUC_indices$TWSS,
  "\n"
)

cat(
  "TotalSS        =",
  COSOHUC_indices$TotalSS,
  "\n"
)

cat(
  "BSS Percentage =",
  COSOHUC_indices$BSS_Percentage,
  "\n"
)

cat(
  "Dunn           =",
  COSOHUC_indices$Dunn,
  "\n"
)

cat(
  "DBI            =",
  COSOHUC_indices$DBI,
  "\n"
)

cat(
  "Silhouette     =",
  COSOHUC_indices$Silhouette,
  "\n"
)


# ================================================================
# 23. ORIGINAL FARTHest FIRST
# ================================================================
#
# Original behavior:
#   - seed is passed to function
#   - set.seed(seed) is inside function
#   - first center is random
#
# ================================================================

farthest_first_centers <- function(
    X,
    k,
    seed
) {

  set.seed(seed)

  unique_X <- unique(X)

  if (nrow(unique_X) < k) {
    stop(
      paste0(
        "Farthest First requires ",
        k,
        " distinct initial centers, but only ",
        nrow(unique_X),
        " unique vectors are available."
      )
    )
  }

  n <- nrow(unique_X)

  selected <- integer(0)

  first <- sample(
    seq_len(n),
    1
  )

  selected <- c(
    selected,
    first
  )

  while (length(selected) < k) {

    min_distances <- rep(
      Inf,
      n
    )

    for (i in seq_len(n)) {

      if (i %in% selected) {

        min_distances[i] <- -Inf

      } else {

        distances <- sapply(
          selected,
          function(s) {
            sqrt(
              sum(
                (
                  unique_X[i, ] -
                    unique_X[s, ]
                )^2
              )
            )
          }
        )

        min_distances[i] <- min(
          distances
        )
      }
    }

    next_point <- which.max(
      min_distances
    )

    selected <- c(
      selected,
      next_point
    )
  }

  centers <- unique_X[
    selected,
    ,
    drop = FALSE
  ]

  colnames(centers) <- colnames(X)

  centers
}


# ================================================================
# 24. ORIGINAL CANOPY
# ================================================================
#
# Original behavior:
#   - seed is passed to function
#   - set.seed(seed) is inside
#   - original Canopy thresholds
#   - original Farthest First fallback with same seed
#
# ================================================================

canopy_initial <- function(
    X,
    k,
    seed
) {
  
  set.seed(
    seed
  )
  
  n <- nrow(X)
  
  D <- as.matrix(
    dist(X)
  )
  
  positive_distances <- D[
    D > 0
  ]
  
  positive_distances <-
    positive_distances[
      is.finite(
        positive_distances
      )
    ]
  
  if (
    length(positive_distances) == 0
  ) {
    
    stop(
      "Cannot calculate Canopy distances."
    )
    
  }
  
  T2 <- as.numeric(
    quantile(
      positive_distances,
      probs = 0.25,
      na.rm = TRUE
    )
  )
  
  T1 <- as.numeric(
    quantile(
      positive_distances,
      probs = 0.50,
      na.rm = TRUE
    )
  )
  
  if (
    T1 <= T2
  ) {
    
    T1 <-
      T2 * 1.5
    
  }
  
  cat(
    "\nCanopy T1 =",
    T1,
    "\n"
  )
  
  cat(
    "Canopy T2 =",
    T2,
    "\n"
  )
  
  remaining <-
    seq_len(n)
  
  selected <-
    integer(
      0
    )
  
  while (
    length(remaining) > 0 &&
    length(selected) < k
  ) {
    
    center_id <-
      remaining[1]
    
    selected <- c(
      selected,
      center_id
    )
    
    distances <-
      D[
        center_id,
        remaining
      ]
    
    remove_ids <-
      remaining[
        distances <= T2
      ]
    
    remaining <-
      setdiff(
        remaining,
        remove_ids
      )
    
  }
  
  
  # --------------------------------------------------------------
  # ORIGINAL FALLBACK
  # --------------------------------------------------------------
  
  if (
    length(selected) < k
  ) {
    
    ff <-
      farthest_first_centers(
        X,
        k,
        seed
      )
    
    for (
      i in 1:nrow(ff)
    ) {
      
      if (
        length(selected) >= k
      ) {
        
        break
        
      }
      
      distances <- apply(
        X,
        1,
        function(row) {
          
          sum(
            (
              row -
                ff[i, ]
            )^2
          )
          
        }
      )
      
      candidate <-
        which.min(
          distances
        )
      
      if (
        !(candidate %in% selected)
      ) {
        
        selected <- c(
          selected,
          candidate
        )
        
      }
      
    }
    
  }
  
  
  if (
    length(selected) < k
  ) {
    
    stop(
      "Canopy could not generate k initial centers."
    )
    
  }
  
  selected <-
    selected[
      1:k
    ]
  
  centers <- X[
    selected,
    ,
    drop = FALSE
  ]
  
  colnames(centers) <-
    colnames(X)
  
  return(
    centers
  )
}


# ================================================================
# 25. ORIGINAL K-MEANS++
# ================================================================
#
# Original behavior:
#   - seed is passed to function
#   - set.seed(seed) is inside
#   - random first center
#   - probability-based subsequent centers
#
# ================================================================

kmeans_plus_plus_centers <- function(
    X,
    k,
    seed
) {

  set.seed(seed)

  unique_X <- unique(X)

  if (nrow(unique_X) < k) {
    stop(
      paste0(
        "K-Means++ requires ",
        k,
        " distinct initial centers, but only ",
        nrow(unique_X),
        " unique vectors are available."
      )
    )
  }

  n <- nrow(unique_X)

  selected <- integer(0)

  first <- sample(
    seq_len(n),
    1
  )

  selected <- c(
    selected,
    first
  )

  while (length(selected) < k) {

    candidates <- setdiff(
      seq_len(n),
      selected
    )

    distances_squared <- numeric(
      length(candidates)
    )

    for (a in seq_along(candidates)) {

      i <- candidates[a]

      distances <- sapply(
        selected,
        function(s) {
          sum(
            (
              unique_X[i, ] -
                unique_X[s, ]
            )^2
          )
        }
      )

      distances_squared[a] <- min(
        distances
      )
    }

    total <- sum(
      distances_squared
    )

    if (
      total <= 0 ||
      !is.finite(total)
    ) {

      next_point <- sample(
        candidates,
        1
      )

    } else {

      probabilities <-
        distances_squared / total

      next_point <- sample(
        candidates,
        1,
        prob = probabilities
      )
    }

    selected <- c(
      selected,
      next_point
    )
  }

  centers <- unique_X[
    selected,
    ,
    drop = FALSE
  ]

  colnames(centers) <- colnames(X)

  centers
}


# ================================================================
# 26. OTHER METHODS: K-MEANS FROM GIVEN CENTERS
# ================================================================
#
# This helper is ONLY for:
#   Farthest First
#   Canopy
#   K-Means++
#
# COSOHUC DOES NOT USE THIS FUNCTION.
#
# ================================================================

run_from_centers <- function(
    X,
    centers
) {
  
  result <- kmeans(
    X,
    centers = centers,
    iter.max = 100,
    nstart = 1,
    algorithm = "Lloyd"
  )
  
  result$initial.centers <-
    centers
  
  return(
    result
  )
}




# ================================================================
# 27. ECKM - FAITHFUL R TRANSLATION OF THE PUBLISHED PROCEDURE
# ================================================================
#
# Reference:
# T. K. Biswas, K. Giri, S. Roy,
# "ECKM: An improved K-means clustering based on computational
# geometry", Expert Systems with Applications, 212, 118862, 2023.
# DOI: 10.1016/j.eswa.2022.118862
#
# The authors' implementation is Python. This is a faithful R
# translation of the supplied published/official procedure.
#
# Sequence preserved:
#   PCA (when p > 2)
#   -> unique geometry points
#   -> Delaunay/Voronoi vertices
#   -> interior vertices
#   -> empty-circle radii
#   -> decreasing-radius LEC order
#   -> ceiling(k/3) circles initially
#   -> circumference points
#   -> round(abs(distance-radius), 3) < epsilon
#   -> unique() preserving first occurrence order
#   -> increase LEC count until > k points
#   -> take first k points
#   -> map to distinct original observations
#   -> standard K-means, nstart=1
#
# Note: R/deldir can differ slightly from the authors' Python geometry
# library because of floating-point and tie-breaking behavior.
# ================================================================

point_inside_convex_hull <- function(point, hull_points) {

  hull_index <- chull(
    hull_points[, 1],
    hull_points[, 2]
  )

  hull <- hull_points[
    c(hull_index, hull_index[1]),
    ,
    drop = FALSE
  ]

  x <- point[1]
  y <- point[2]
  inside <- FALSE

  n <- nrow(hull) - 1
  j <- n

  if (n < 3) {
    return(FALSE)
  }

  for (i in seq_len(n)) {

    xi <- hull[i, 1]
    yi <- hull[i, 2]

    xj <- hull[j, 1]
    yj <- hull[j, 2]

    intersects <-
      ((yi > y) != (yj > y)) &&
      (
        x <
          (
            (xj - xi) *
              (y - yi) /
              (yj - yi + .Machine$double.eps)
          ) +
          xi
      )

    if (intersects) {
      inside <- !inside
    }

    j <- i
  }

  inside
}


calculate_voronoi_vertices_eckm <- function(points) {

  if (nrow(points) < 3) {
    stop(
      "ECKM requires at least 3 unique 2-D points."
    )
  }

  delaunay <- deldir::deldir(
    points[, 1],
    points[, 2],
    suppressMsge = TRUE
  )

  tiles <- deldir::tile.list(
    delaunay
  )

  voronoi_vertices <- list()
  counter <- 0L

  for (tile in tiles) {

    tile_x <- tile$x
    tile_y <- tile$y

    if (length(tile_x) < 3) {
      next
    }

    ax <- tile_x[1]
    ay <- tile_y[1]

    bx <- tile_x[2]
    by <- tile_y[2]

    cx <- tile_x[3]
    cy <- tile_y[3]

    denominator <-
      2 *
      (
        ax * (by - cy) +
        bx * (cy - ay) +
        cx * (ay - by)
      )

    if (
      abs(denominator) <
      .Machine$double.eps
    ) {
      next
    }

    ux <-
      (
        (ax^2 + ay^2) * (by - cy) +
        (bx^2 + by^2) * (cy - ay) +
        (cx^2 + cy^2) * (ay - by)
      ) /
      denominator

    uy <-
      (
        (ax^2 + ay^2) * (cx - bx) +
        (bx^2 + by^2) * (ax - cx) +
        (cx^2 + cy^2) * (bx - ax)
      ) /
      denominator

    counter <- counter + 1L

    voronoi_vertices[[counter]] <-
      c(
        ux,
        uy
      )
  }

  if (
    length(voronoi_vertices) == 0
  ) {
    stop(
      "No Voronoi vertices could be calculated."
    )
  }

  vertices <- do.call(
    rbind,
    voronoi_vertices
  )

  vertices <- unique(
    round(
      vertices,
      12
    )
  )

  vertices
}


prepare_eckm_geometry <- function(
    X,
    use_pca = TRUE
) {

  if (ncol(X) < 2) {
    stop(
      "ECKM requires at least two attributes."
    )
  }

  if (
    ncol(X) > 2 &&
    use_pca
  ) {

    pca_result <- prcomp(
      X,
      center = TRUE,
      scale. = FALSE
    )

    geometry_data <- pca_result$x[
      ,
      1:2,
      drop = FALSE
    ]

  } else {

    geometry_data <- X[
      ,
      1:2,
      drop = FALSE
    ]
  }

  colnames(geometry_data) <- c(
    "PC1",
    "PC2"
  )

  duplicate_geometry <- duplicated(
    geometry_data
  )

  geometry_data_unique <-
    geometry_data[
      !duplicate_geometry,
      ,
      drop = FALSE
    ]

  if (
    nrow(geometry_data_unique) < 3
  ) {
    stop(
      "Insufficient unique points for Voronoi construction."
    )
  }

  list(
    geometry_data = geometry_data,
    geometry_data_unique =
      geometry_data_unique
  )
}


get_eckm_interior_vertices <- function(
    geometry_data_unique
) {

  voronoi_vertices <-
    calculate_voronoi_vertices_eckm(
      geometry_data_unique
    )

  hull_index <- chull(
    geometry_data_unique[, 1],
    geometry_data_unique[, 2]
  )

  hull_points <- geometry_data_unique[
    hull_index,
    ,
    drop = FALSE
  ]

  interior_vertices <-
    matrix(
      numeric(0),
      ncol = 2
    )

  for (
    i in seq_len(
      nrow(voronoi_vertices)
    )
  ) {

    point <- voronoi_vertices[
      i,
      ]

    if (
      point_inside_convex_hull(
        point,
        hull_points
      )
    ) {

      interior_vertices <-
        rbind(
          interior_vertices,
          point
        )
    }
  }

  if (
    nrow(interior_vertices) == 0
  ) {
    stop(
      "No interior Voronoi vertices were found."
    )
  }

  unique(
    round(
      interior_vertices,
      12
    )
  )
}


calculate_radii_eckm <- function(
    centers,
    points
) {

  radius_values <- numeric(
    nrow(centers)
  )

  for (
    i in seq_len(
      nrow(centers)
    )
  ) {

    distances <- sqrt(

      (
        points[, 1] -
          centers[i, 1]
      )^2 +

      (
        points[, 2] -
          centers[i, 2]
      )^2
    )

    radius_values[i] <-
      min(distances)
  }

  radius_values
}


find_circumference_points_eckm <- function(
    circle_centers,
    circle_radii,
    data_points,
    k,
    epsilon = 0.001
) {

  number_requested <-
    ceiling(
      k / 3
    )

  repeat {

    if (
      number_requested >
      length(circle_radii)
    ) {

      stop(
        paste0(
          "Not enough Voronoi circles available to construct ",
          k,
          " initial centers."
        )
      )
    }

    current_centers <-
      circle_centers[
        seq_len(number_requested),
        ,
        drop = FALSE
      ]

    current_radii <-
      circle_radii[
        seq_len(number_requested)
      ]

    all_indices <- integer(
      0
    )

    for (
      i in seq_len(
        nrow(current_centers)
      )
    ) {

      distances <- sqrt(

        (
          data_points[, 1] -
            current_centers[i, 1]
        )^2 +

        (
          data_points[, 2] -
            current_centers[i, 2]
        )^2
      )

      difference <- abs(
        distances -
          current_radii[i]
      )

      difference_rounded <- round(
        difference,
        3
      )

      point_indices <- which(
        difference_rounded <
          epsilon
      )

      if (
        length(point_indices) > 0
      ) {

        all_indices <- c(
          all_indices,
          point_indices
        )
      }
    }

    # unique() preserves first occurrence order in R.
    all_indices <- unique(
      all_indices
    )

    # Preserve the supplied implementation's stopping rule.
    if (
      length(all_indices) > k
    ) {
      break
    }

    number_requested <-
      number_requested + 1L
  }

  list(
    selected_count =
      number_requested,

    point_indices =
      all_indices
  )
}


map_eckm_geometry_to_original <- function(
    geometry_data,
    eckm_geometry_centers,
    k
) {

  selected_original_indices <-
    integer(0)

  for (
    i in seq_len(
      nrow(eckm_geometry_centers)
    )
  ) {

    distances <- sqrt(

      (
        geometry_data[, 1] -
          eckm_geometry_centers[i, 1]
      )^2 +

      (
        geometry_data[, 2] -
          eckm_geometry_centers[i, 2]
      )^2
    )

    candidate_order <- order(
      distances
    )

    candidate <-
      candidate_order[1]

    if (
      candidate %in%
      selected_original_indices
    ) {

      found_new <- FALSE

      for (
        cand in candidate_order
      ) {

        if (
          !(cand %in%
            selected_original_indices)
        ) {

          candidate <- cand
          found_new <- TRUE
          break
        }
      }

      if (!found_new) {

        stop(
          "ECKM could not map selected geometry centers to distinct original observations."
        )
      }
    }

    selected_original_indices <-
      c(
        selected_original_indices,
        candidate
      )
  }

  if (
    length(
      unique(
        selected_original_indices
      )
    ) < k
  ) {

    stop(
      "ECKM could not map k distinct geometry centers to data observations."
    )
  }

  selected_original_indices[
    seq_len(k)
  ]
}


ECKM_initial_centers <- function(
    X,
    k,
    epsilon = 0.001,
    use_pca = TRUE
) {

  if (
    nrow(X) < k
  ) {

    stop(
      "Number of observations must be >= k."
    )
  }

  geometry <-
    prepare_eckm_geometry(
      X,
      use_pca =
        use_pca
    )

  geometry_data <-
    geometry$geometry_data

  geometry_data_unique <-
    geometry$geometry_data_unique

  interior_vertices <-
    get_eckm_interior_vertices(
      geometry_data_unique
    )

  radius_values <-
    calculate_radii_eckm(
      interior_vertices,
      geometry_data_unique
    )

  # Largest Empty Circles: descending radius.
  lec_order <-
    order(
      radius_values,
      decreasing = TRUE
    )

  lec_centers <-
    interior_vertices[
      lec_order,
      ,
      drop = FALSE
    ]

  lec_radii <-
    radius_values[
      lec_order
    ]

  lec_result <-
    find_circumference_points_eckm(

      circle_centers =
        lec_centers,

      circle_radii =
        lec_radii,

      data_points =
        geometry_data_unique,

      k =
        k,

      epsilon =
        epsilon
    )

  lec_indices <-
    lec_result$point_indices

  if (
    length(lec_indices) < k
  ) {

    stop(
      "ECKM could not obtain enough circumference points."
    )
  }

  # Corresponds to circumpts_all[:noc].
  lec_indices <-
    lec_indices[
      seq_len(k)
    ]

  eckm_geometry_centers <-
    geometry_data_unique[
      lec_indices,
      ,
      drop = FALSE
    ]

  selected_original_indices <-
    map_eckm_geometry_to_original(

      geometry_data =
        geometry_data,

      eckm_geometry_centers =
        eckm_geometry_centers,

      k =
        k
    )

  initial_centers <-
    X[
      selected_original_indices,
      ,
      drop = FALSE
    ]

  colnames(initial_centers) <-
    colnames(X)

  list(

    centers =
      initial_centers,

    selected_rows =
      selected_original_indices,

    geometry =
      geometry_data,

    geometry_unique =
      geometry_data_unique,

    interior_vertices =
      interior_vertices,

    radius_values =
      radius_values,

    lec_order =
      lec_order,

    lec_centers =
      lec_centers,

    lec_radii =
      lec_radii,

    selected_lec_count =
      lec_result$selected_count,

    circumference_indices =
      lec_indices
  )
}


run_ECKM <- function(
    X,
    k,
    epsilon = 0.001,
    use_pca = TRUE
) {

  initialization <-
    ECKM_initial_centers(

      X =
        X,

      k =
        k,

      epsilon =
        epsilon,

      use_pca =
        use_pca
    )

  # Authors' Python code uses n_init=1.
  # R equivalent is nstart=1.
  result <-
    kmeans(

      X,

      centers =
        initialization$centers,

      nstart =
        1,

      iter.max =
        300,

      algorithm =
        "Lloyd"
    )

  result$initial.centers <-
    initialization$centers

  result$ECKM_selected_rows <-
    initialization$selected_rows

  result$ECKM_selected_lec_count <-
    initialization$selected_lec_count

  result$ECKM_circumference_indices <-
    initialization$circumference_indices

  result$ECKM_radius_values <-
    initialization$radius_values

  result$ECKM_geometry <-
    initialization$geometry

  result
}


# ================================================================
# 28. CALCULATE ECKM ONCE
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "CALCULATING DETERMINISTIC ECKM\n"
)

cat(
  "============================================================\n"
)

ECKM_result <-
  run_ECKM(
    X,
    k
  )

ECKM_indices <-
  calculate_all_indices(
    X,
    ECKM_result
  )

cat(
  "\nECKM selected original rows:\n"
)

print(
  ECKM_result$ECKM_selected_rows
)

cat(
  "\nECKM LEC count used:\n"
)

print(
  ECKM_result$ECKM_selected_lec_count
)

cat(
  "\nECKM result:\n"
)

cat(
  "Iteration =",
  ECKM_indices$Iteration,
  "\n"
)

cat(
  "BSS       =",
  ECKM_indices$BSS,
  "\n"
)

cat(
  "TWSS      =",
  ECKM_indices$TWSS,
  "\n"
)

cat(
  "TotalSS   =",
  ECKM_indices$TotalSS,
  "\n"
)

cat(
  "Dunn      =",
  ECKM_indices$Dunn,
  "\n"
)

cat(
  "DBI       =",
  ECKM_indices$DBI,
  "\n"
)

cat(
  "Silhouette =",
  ECKM_indices$Silhouette,
  "\n"
)

# ================================================================
# 53. MINI-BATCH K-MEANS++ IMPLEMENTATION
# ================================================================
#
# Integrated from the standalone Mini-Batch K-Means++ implementation.
# The existing COSOHUC, ECKM, and other initialization methods are not
# modified by this addition.
#
# ================================================================

precompute_batch_terms <- function(
    X,
    b
) {

  n <- nrow(X)
  d <- ncol(X)

  B <- ceiling(
    n / b
  )

  cumulative_coordinate_sums <- matrix(
    0,
    nrow = B,
    ncol = d
  )

  cumulative_squared_sums <- numeric(
    B
  )

  batch_start <- integer(B)
  batch_end <- integer(B)

  running_coordinate_sum <- rep(
    0,
    d
  )

  running_squared_sum <- 0

  for (
    batch_id in seq_len(B)
  ) {

    start_id <- (
      batch_id - 1
    ) * b + 1

    end_id <- min(
      batch_id * b,
      n
    )

    batch_data <- X[
      start_id:end_id,
      ,
      drop = FALSE
    ]

    current_coordinate_sum <- colSums(
      batch_data
    )

    current_squared_sum <- sum(
      batch_data^2
    )

    running_coordinate_sum <-
      running_coordinate_sum +
      current_coordinate_sum

    running_squared_sum <-
      running_squared_sum +
      current_squared_sum

    cumulative_coordinate_sums[
      batch_id,
      ] <-
      running_coordinate_sum

    cumulative_squared_sums[
      batch_id
    ] <-
      running_squared_sum

    batch_start[
      batch_id
    ] <-
      start_id

    batch_end[
      batch_id
    ] <-
      end_id
  }

  list(
    B = B,
    b = b,
    cumulative_coordinate_sums =
      cumulative_coordinate_sums,
    cumulative_squared_sums =
      cumulative_squared_sums,
    batch_start =
      batch_start,
    batch_end =
      batch_end
  )
}


cumulative_distance_sum_to_batch <- function(
    selected_centers,
    cumulative_coordinate_sum,
    cumulative_squared_sum,
    cumulative_n
) {

  if (
    is.null(
      dim(selected_centers)
    )
  ) {

    selected_centers <- matrix(
      selected_centers,
      nrow = 1
    )
  }

  number_centers <- nrow(
    selected_centers
  )

  center_squared_sum <- sum(
    selected_centers^2
  )

  centers_coordinate_sum <- colSums(
    selected_centers
  )

  value <-
    cumulative_n *
    center_squared_sum -
    2 *
    sum(
      cumulative_coordinate_sum *
      centers_coordinate_sum
    ) +
    number_centers *
    cumulative_squared_sum

  if (
    value < 0 &&
    value > -1e-10
  ) {
    value <- 0
  }

  max(
    0,
    value
  )
}


find_target_batch <- function(
    R,
    selected_centers,
    batch_terms
) {

  B <- batch_terms$B

  if (R <= 0) {
    return(1L)
  }

  left <- 1L
  right <- B

  while (left < right) {

    middle <- floor(
      (left + right) / 2
    )

    cumulative_n <- batch_terms$batch_end[
      middle
    ]

    cumulative_coordinate_sum <-
      batch_terms$cumulative_coordinate_sums[
        middle,
        ,
        drop = TRUE
      ]

    cumulative_squared_sum <-
      batch_terms$cumulative_squared_sums[
        middle
      ]

    cumulative_distance <-
      cumulative_distance_sum_to_batch(

        selected_centers =
          selected_centers,

        cumulative_coordinate_sum =
          cumulative_coordinate_sum,

        cumulative_squared_sum =
          cumulative_squared_sum,

        cumulative_n =
          cumulative_n
      )

    if (
      cumulative_distance >= R
    ) {

      right <- middle

    } else {

      left <- middle + 1L
    }
  }

  left
}


find_point_inside_batch <- function(
    X,
    selected_centers,
    batch_start,
    batch_end,
    R,
    distance_before_batch
) {

  batch_ids <- batch_start:batch_end

  target_inside_batch <-
    R -
    distance_before_batch

  cumulative_inside_batch <- 0

  for (
    local_id in seq_along(batch_ids)
  ) {

    global_id <- batch_ids[
      local_id
    ]

    point <- X[
      global_id,
      ,
      drop = TRUE
    ]

    point_distance <- 0

    for (
      center_id in seq_len(
        nrow(selected_centers)
      )
    ) {

      difference <-
        point -
        selected_centers[
          center_id,
          ,
          drop = TRUE
        ]

      point_distance <-
        point_distance +
        sum(
          difference^2
        )
    }

    cumulative_inside_batch <-
      cumulative_inside_batch +
      point_distance

    if (
      cumulative_inside_batch >=
      target_inside_batch
    ) {

      return(
        global_id
      )
    }
  }

  # Floating-point boundary fallback.
  batch_ids[
    length(batch_ids)
  ]
}


mini_batch_kmeans_plus_plus_centers <- function(
    X,
    k,
    b,
    seed
) {

  set.seed(seed)

  n <- nrow(X)

  if (k > n) {
    stop(
      "k cannot exceed number of observations."
    )
  }

  batch_terms <- precompute_batch_terms(
    X,
    b
  )

  # --------------------------------------------------------------
  # First center
  # --------------------------------------------------------------

  first_index <- sample(
    seq_len(n),
    size = 1
  )

  selected_indices <- first_index

  selected_centers <- X[
    first_index,
    ,
    drop = FALSE
  ]

  # --------------------------------------------------------------
  # Select remaining centers
  # --------------------------------------------------------------

  while (
    length(selected_indices) < k
  ) {

    total_distance <-
      cumulative_distance_sum_to_batch(

        selected_centers =
          selected_centers,

        cumulative_coordinate_sum =
          batch_terms$
          cumulative_coordinate_sums[
            batch_terms$B,
            ,
            drop = TRUE
          ],

        cumulative_squared_sum =
          batch_terms$
          cumulative_squared_sums[
            batch_terms$B
          ],

        cumulative_n =
          n
      )

    # ------------------------------------------------------------
    # Degenerate case:
    # all distances are zero.
    # ------------------------------------------------------------

    if (
      total_distance <= 0 ||
      !is.finite(total_distance)
    ) {

      remaining <- setdiff(
        seq_len(n),
        selected_indices
      )

      next_index <- sample(
        remaining,
        size = 1
      )

      selected_indices <- c(
        selected_indices,
        next_index
      )

      selected_centers <- rbind(
        selected_centers,
        X[
          next_index,
          ,
          drop = FALSE
        ]
      )

      next
    }

    # ------------------------------------------------------------
    # Distance-proportional random number R.
    # ------------------------------------------------------------

    R <- runif(
      1,
      min = 0,
      max = total_distance
    )

    # ------------------------------------------------------------
    # Locate target batch.
    # ------------------------------------------------------------

    target_batch <- find_target_batch(
      R =
        R,
      selected_centers =
        selected_centers,
      batch_terms =
        batch_terms
    )

    # ------------------------------------------------------------
    # Distance accumulated before target batch.
    # ------------------------------------------------------------

    if (
      target_batch == 1
    ) {

      distance_before_batch <- 0

    } else {

      cumulative_n_before <-
        batch_terms$batch_end[
          target_batch - 1
        ]

      cumulative_coordinate_sum_before <-
        batch_terms$cumulative_coordinate_sums[
          target_batch - 1,
          ,
          drop = TRUE
        ]

      cumulative_squared_sum_before <-
        batch_terms$cumulative_squared_sums[
          target_batch - 1
        ]

      distance_before_batch <-
        cumulative_distance_sum_to_batch(

          selected_centers =
            selected_centers,

          cumulative_coordinate_sum =
            cumulative_coordinate_sum_before,

          cumulative_squared_sum =
            cumulative_squared_sum_before,

          cumulative_n =
            cumulative_n_before
        )
    }

    # ------------------------------------------------------------
    # Search inside selected batch.
    # ------------------------------------------------------------

    batch_start <- batch_terms$batch_start[
      target_batch
    ]

    batch_end <- batch_terms$batch_end[
      target_batch
    ]

    next_index <- find_point_inside_batch(

      X =
        X,

      selected_centers =
        selected_centers,

      batch_start =
        batch_start,

      batch_end =
        batch_end,

      R =
        R,

      distance_before_batch =
        distance_before_batch
    )

    # ------------------------------------------------------------
    # Never duplicate a selected center.
    # This should occur only because of floating-point boundaries.
    # ------------------------------------------------------------

    if (
      next_index %in%
      selected_indices
    ) {

      remaining_in_batch <- setdiff(
        batch_start:batch_end,
        selected_indices
      )

      if (
        length(remaining_in_batch) > 0
      ) {

        scores <- numeric(
          length(remaining_in_batch)
        )

        for (
          z in seq_along(
            remaining_in_batch
          )
        ) {

          point <- X[
            remaining_in_batch[z],
            ,
            drop = TRUE
          ]

          point_score <- 0

          for (
            center_id in seq_len(
              nrow(selected_centers)
            )
          ) {

            difference <-
              point -
              selected_centers[
                center_id,
                ,
                drop = TRUE
              ]

            point_score <-
              point_score +
              sum(
                difference^2
              )
          }

          scores[z] <- point_score
        }

        next_index <-
          remaining_in_batch[
            which.max(scores)
          ]

      } else {

        remaining <- setdiff(
          seq_len(n),
          selected_indices
        )

        next_index <- sample(
          remaining,
          size = 1
        )
      }
    }

    selected_indices <- c(
      selected_indices,
      next_index
    )

    selected_centers <- rbind(
      selected_centers,
      X[
        next_index,
        ,
        drop = FALSE
      ]
    )
  }

  colnames(
    selected_centers
  ) <- colnames(X)

  list(
    indices =
      selected_indices,

    centers =
      selected_centers,

    batch_terms =
      batch_terms
  )
}


run_mini_batch_kmeans_plus_plus <- function(
    X,
    k,
    b,
    seed
) {

  initialization_start <-
    proc.time()[["elapsed"]]

  initialization <-
    mini_batch_kmeans_plus_plus_centers(

      X =
        X,

      k =
        k,

      b =
        b,

      seed =
        seed
    )

  initialization_time <-
    proc.time()[["elapsed"]] -
    initialization_start

  # --------------------------------------------------------------
  # Standard K-means after MBKM++ initialization.
  # --------------------------------------------------------------

  clustering_start <-
    proc.time()[["elapsed"]]

  result <- kmeans(

    X,

    centers =
      initialization$centers,

    nstart =
      1,

    iter.max =
      100,

    algorithm =
      "Lloyd"
  )

  clustering_time <-
    proc.time()[["elapsed"]] -
    clustering_start

  result$initial.centers <-
    initialization$centers

  result$initial.indices <-
    initialization$indices

  result$initialization.time <-
    initialization_time

  result$clustering.time <-
    clustering_time

  result$total.time <-
    initialization_time +
    clustering_time

  result$batch.size <-
    b

  result$number.batches <-
    initialization$batch_terms$B

  result
}


# ================================================================
# METHOD KEYS - DEFINE BEFORE RUN LOOP
# ================================================================

method_names <- c(
  "Farthest_First",
  "Canopy",
  "K_Means_Plus_Plus",
  "Hartigan_Wong",
  "MacQueen",
  "Lloyd",
  "Forgy",
  "ECKM",
  "Mini_Batch_K_Means_Plus_Plus",
  "COSOHUC"
)

method_display_names <- c(
  "Farthest First",
  "Canopy",
  "K-Means++",
  "Hartigan-Wong",
  "MacQueen",
  "Lloyd",
  "Forgy",
  "ECKM",
  "Mini-Batch K-Means++",
  "COSOHUC"
)


# ================================================================
# 54. RUN N EXPERIMENTS
# ================================================================
# ================================================================
#
# COSOHUC and ECKM are deterministic and are NOT recalculated.
# The stochastic methods, including Mini-Batch K-Means++, use a different fixed seed for every run.
#
# ================================================================
#
# COSOHUC is NOT recalculated.
#
# ================================================================

all_runs <-
  vector(
    "list",
    N
  )


for (
  run_number in seq_len(N)
) {
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "RUN ",
    run_number,
    " OF ",
    N,
    "\n",
    sep = ""
  )
  
  cat(
    "============================================================\n"
  )
  
  
  # --------------------------------------------------------------
  # Seed for stochastic methods
  # --------------------------------------------------------------
  
  seed <-
    1000 +
    run_number
  
  
  # ==============================================================
  # MINI-BATCH K-MEANS++
  # ==============================================================
  
  MBKMPP_result <-
    run_mini_batch_kmeans_plus_plus(
      X,
      k,
      b,
      seed
    )
  
  
  # ==============================================================
  # FARTHEST FIRST
  # ==============================================================
  
  farthest_initial <-
    farthest_first_centers(
      X,
      k,
      seed
    )
  
  result_farthest <-
    run_from_centers(
      X,
      farthest_initial
    )
  
  
  # ==============================================================
  # CANOPY
  # ==============================================================
  
  canopy_initial_centers <-
    canopy_initial(
      X,
      k,
      seed
    )
  
  result_canopy <-
    run_from_centers(
      X,
      canopy_initial_centers
    )
  
  
  # ==============================================================
  # K-MEANS++
  # ==============================================================
  
  kmeanspp_initial <-
    kmeans_plus_plus_centers(
      X,
      k,
      seed
    )
  
  result_kmeanspp <-
    run_from_centers(
      X,
      kmeanspp_initial
    )
  
  
  # ==============================================================
  # HARTIGAN-WONG
  # ==============================================================
  
  set.seed(
    seed
  )
  
  result_HW <-
    kmeans(
      X,
      centers = k,
      iter.max = 100,
      nstart = 1,
      algorithm = "Hartigan-Wong"
    )
  
  
  # ==============================================================
  # MACQUEEN
  # ==============================================================
  
  set.seed(
    seed
  )
  
  result_MacQueen <-
    kmeans(
      X,
      centers = k,
      iter.max = 100,
      nstart = 1,
      algorithm = "MacQueen"
    )
  
  
  # ==============================================================
  # LLOYD
  # ==============================================================
  
  set.seed(
    seed
  )
  
  result_Lloyd <-
    kmeans(
      X,
      centers = k,
      iter.max = 100,
      nstart = 1,
      algorithm = "Lloyd"
    )
  
  
  # ==============================================================
  # FORGY
  # ==============================================================

  set.seed(
    seed
  )

  unique_X <- unique(X)

  if (nrow(unique_X) < k) {
    stop(
      paste0(
        "RUN ",
        run_number,
        ": Forgy cannot select ",
        k,
        " distinct initial centers because only ",
        nrow(unique_X),
        " unique observation vectors are available."
      )
    )
  }

  forgy_indices <- sample(
    seq_len(nrow(unique_X)),
    size = k,
    replace = FALSE
  )

  forgy_initial <- unique_X[
    forgy_indices,
    ,
    drop = FALSE
  ]

  result_Forgy <- kmeans(
    X,
    centers = forgy_initial,
    iter.max = 100,
    nstart = 1,
    algorithm = "Lloyd"
  )

  result_Forgy$initial.centers <- forgy_initial


# ==============================================================
  # STORE ALL METHODS
  # ==============================================================
  
  all_runs[[run_number]] <-
    list(
      
      # ------------------------------------------------------------
      # VERY IMPORTANT:
      # This is the already-calculated deterministic COSOHUC result.
      # It is NOT recalculated.
      # ------------------------------------------------------------
      
      COSOHUC =
        COSOHUC_result,
      
      Farthest_First =
        result_farthest,
      
      Canopy =
        result_canopy,
      
      K_Means_Plus_Plus =
        result_kmeanspp,
      
      Hartigan_Wong =
        result_HW,
      
      MacQueen =
        result_MacQueen,
      
      Lloyd =
        result_Lloyd,
      
      Forgy =
        result_Forgy,
      
      ECKM =
        ECKM_result,
      
      Mini_Batch_K_Means_Plus_Plus =
        MBKMPP_result
      
    )
  
  
  # ==============================================================
  # DISPLAY ITERATIONS
  # ==============================================================
  
  cat(
    "\nIterations:\n"
  )
  
  cat(
    "COSOHUC        =",
    COSOHUC_result$iter,
    "\n"
  )
  
  cat(
    "Farthest First =",
    result_farthest$iter,
    "\n"
  )
  
  cat(
    "Canopy         =",
    result_canopy$iter,
    "\n"
  )
  
  cat(
    "K-Means++      =",
    result_kmeanspp$iter,
    "\n"
  )
  
  cat(
    "Hartigan-Wong  =",
    result_HW$iter,
    "\n"
  )
  
  cat(
    "MacQueen       =",
    result_MacQueen$iter,
    "\n"
  )
  
  cat(
    "Lloyd          =",
    result_Lloyd$iter,
    "\n"
  )
  
  cat(
    "Forgy          =",
    result_Forgy$iter,
    "\n"
  )
  
  cat(
    "ECKM           =",
    ECKM_result$iter,
    "\n"
  )

  cat(
    "Mini-Batch K-Means++ =",
    MBKMPP_result$iter,
    "\n"
  )
  
}


# ================================================================
# 29. RAW RESULTS TABLE
# ================================================================

raw_results <- data.frame(
  
  Run = integer(),
  
  Method = character(),
  
  Iteration = numeric(),
  
  BSS = numeric(),
  
  TWSS = numeric(),
  
  TotalSS = numeric(),
  
  BSS_Percentage = numeric(),
  
  Dunn = numeric(),
  
  DBI = numeric(),
  
  Silhouette = numeric(),
  
  Rand_Index = numeric(),
  
  stringsAsFactors = FALSE
  
)


# ================================================================
# 30. EXTRACT RESULTS
# ================================================================

for (
  run_number in seq_len(N)
) {
  
  current_results <-
    all_runs[[run_number]]
  
  
  for (
    method_name in method_names
  ) {
    
    # ------------------------------------------------------------
    # IMPORTANT:
    # Correct R list extraction syntax.
    # ------------------------------------------------------------
    
    result <-
      current_results[[method_name]]
    
    
    indices <-
      calculate_all_indices(
        X,
        result
      )
    
    
    # ------------------------------------------------------------
    # RAND INDEX
    # ------------------------------------------------------------
    
    if (
      calculate_RI
    ) {
      
      rand_value <-
        calculate_rand_index(
          true_labels,
          result$cluster
        )
      
    } else {
      
      rand_value <-
        NA_real_
      
    }
    
    
    # ------------------------------------------------------------
    # CREATE RESULT ROW
    # ------------------------------------------------------------
    
    new_row <- data.frame(
      
      Run =
        run_number,
      
      Method =
        method_name,
      
      Iteration =
        indices$Iteration,
      
      BSS =
        indices$BSS,
      
      TWSS =
        indices$TWSS,
      
      TotalSS =
        indices$TotalSS,
      
      BSS_Percentage =
        indices$BSS_Percentage,
      
      Dunn =
        indices$Dunn,
      
      DBI =
        indices$DBI,
      
      Silhouette =
        indices$Silhouette,
      
      Rand_Index =
        rand_value,
      
      stringsAsFactors = FALSE
      
    )
    
    
    raw_results <-
      rbind(
        raw_results,
        new_row
      )
    
  }
  
}


# ================================================================
# 31. METHOD DISPLAY NAMES
# ================================================================

raw_results$Method <- factor(
  
  raw_results$Method,
  
  levels = method_names,
  
  labels = method_display_names
  
)


raw_results <-
  raw_results[
    order(
      raw_results$Run,
      raw_results$Method
    ),
    ,
    drop = FALSE
  ]


raw_results$Method <-
  as.character(
    raw_results$Method
  )


# ================================================================
# 32. COSOHUC DETERMINISM CHECK
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "COSOHUC DETERMINISM CHECK\n"
)

cat(
  "============================================================\n"
)


cosohuc_rows <-
  raw_results[
    raw_results$Method ==
      "COSOHUC",
    ,
    drop = FALSE
  ]


if (
  nrow(cosohuc_rows) > 1
) {
  
  check_columns <- c(
    
    "Iteration",
    
    "BSS",
    
    "TWSS",
    
    "TotalSS",
    
    "BSS_Percentage",
    
    "Dunn",
    
    "DBI",
    
    "Silhouette"
    
  )
  
  
  first_values <-
    as.numeric(
      cosohuc_rows[
        1,
        check_columns,
        drop = FALSE
      ]
    )
  
  
  deterministic <- TRUE
  
  
  for (
    i in 2:nrow(cosohuc_rows)
  ) {
    
    current_values <-
      as.numeric(
        cosohuc_rows[
          i,
          check_columns,
          drop = FALSE
        ]
      )
    
    
    different_NA <-
      is.na(first_values) !=
      is.na(current_values)
    
    
    different_values <-
      !is.na(first_values) &
      !is.na(current_values) &
      (
        abs(
          first_values -
            current_values
        ) > 1e-12
      )
    
    
    if (
      any(different_NA) ||
      any(different_values)
    ) {
      
      deterministic <- FALSE
      
      break
      
    }
    
  }
  
  
  if (
    !deterministic
  ) {
    
    stop(
      paste0(
        "ERROR: COSOHUC changed between runs.\n",
        "The COSOHUC result must be identical for every run."
      )
    )
    
  }
  
}


cat(
  "PASS: COSOHUC is identical in all ",
  N,
  " runs.\n",
  sep = ""
)


# ================================================================
# ECKM DETERMINISM CHECK
# ================================================================

eckm_rows <- raw_results[
  raw_results$Method == "ECKM",
  ,
  drop = FALSE
]

if (nrow(eckm_rows) > 1) {
  eckm_check_columns <- c(
    "Iteration",
    "BSS",
    "TWSS",
    "TotalSS",
    "BSS_Percentage",
    "Dunn",
    "DBI",
    "Silhouette"
  )
  
  first_values <- as.numeric(
    eckm_rows[1, eckm_check_columns, drop = FALSE]
  )
  
  for (i in 2:nrow(eckm_rows)) {
    current_values <- as.numeric(
      eckm_rows[i, eckm_check_columns, drop = FALSE]
    )
    
    different_na <-
      is.na(first_values) != is.na(current_values)
    
    different_values <-
      !is.na(first_values) &
      !is.na(current_values) &
      abs(first_values - current_values) > 1e-12
    
    if (any(different_na) || any(different_values)) {
      stop(
        "ERROR: ECKM changed between runs. The ECKM result must be identical for a fixed dataset and k."
      )
    }
  }
}

cat(
  "PASS: ECKM is identical in all ",
  N,
  " runs.\n",
  sep = ""
)


# ================================================================
# 35. DISPLAY RAW RESULTS
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "RAW RESULTS - ALL ",
  N,
  " RUNS\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)


print(
  raw_results,
  row.names = FALSE
)


# ================================================================
# 34. NUMERIC VARIABLES
# ================================================================

numeric_variables <- c(
  
  "Iteration",
  
  "BSS",
  
  "TWSS",
  
  "TotalSS",
  
  "BSS_Percentage",
  
  "Dunn",
  
  "DBI",
  
  "Silhouette",
  
  "Rand_Index"
  
)


# ================================================================
# 35. AVERAGE RESULTS
# ================================================================

average_results <-
  aggregate(
    
    raw_results[
      ,
      numeric_variables,
      drop = FALSE
    ],
    
    by = list(
      Method =
        raw_results$Method
    ),
    
    FUN = function(x) {
      
      if (
        all(is.na(x))
      ) {
        
        return(
          NA_real_
        )
        
      }
      
      mean(
        x,
        na.rm = TRUE
      )
      
    }
    
  )


colnames(
  average_results
)[1] <-
  "Method"


# ================================================================
# 36. STANDARD DEVIATION
# ================================================================

sd_results <-
  aggregate(
    
    raw_results[
      ,
      numeric_variables,
      drop = FALSE
    ],
    
    by = list(
      Method =
        raw_results$Method
    ),
    
    FUN = function(x) {
      
      valid_values <-
        x[
          !is.na(x)
        ]
      
      
      if (
        length(valid_values) < 2
      ) {
        
        return(
          0
        )
        
      }
      
      
      sd(
        valid_values
      )
      
    }
    
  )


colnames(
  sd_results
)[1] <-
  "Method"


# ================================================================
# 37. FORCE METHOD ORDER
# ================================================================

average_results$Method <-
  factor(
    average_results$Method,
    levels =
      method_display_names
  )


average_results <-
  average_results[
    order(
      average_results$Method
    ),
    ,
    drop = FALSE
  ]


average_results$Method <-
  as.character(
    average_results$Method
  )


sd_results$Method <-
  factor(
    sd_results$Method,
    levels =
      method_display_names
  )


sd_results <-
  sd_results[
    order(
      sd_results$Method
    ),
    ,
    drop = FALSE
  ]


sd_results$Method <-
  as.character(
    sd_results$Method
  )


# ================================================================
# 38. ROUND RESULTS
# ================================================================

average_results[
  ,
  -1
] <-
  round(
    average_results[
      ,
      -1
    ],
    6
  )


sd_results[
  ,
  -1
] <-
  round(
    sd_results[
      ,
      -1
    ],
    6
  )


# ================================================================
# 39. MEAN +/- SD
# ================================================================

mean_sd_results <-
  data.frame(
    
    Method =
      average_results$Method,
    
    stringsAsFactors = FALSE
    
  )


for (
  variable in numeric_variables
) {
  
  mean_values <-
    average_results[
      ,
      variable
    ]
  
  
  sd_values <-
    sd_results[
      ,
      variable
    ]
  
  
  mean_sd_results[
    ,
    variable
  ] <-
    paste0(
      mean_values,
      " +/- ",
      sd_values
    )
  
}


# ================================================================
# 40. DISPLAY AVERAGE
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "AVERAGE RESULTS OVER ",
  N,
  " RUNS\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

print(
  average_results,
  row.names = FALSE
)


# ================================================================
# 41. DISPLAY SD
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "STANDARD DEVIATION OVER ",
  N,
  " RUNS\n",
  sep = ""
)

cat(
  "============================================================\n\n"
)

print(
  sd_results,
  row.names = FALSE
)


# ================================================================
# 42. DISPLAY MEAN +/- SD
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "MEAN +/- SD RESULTS\n"
)

cat(
  "============================================================\n\n"
)

print(
  mean_sd_results,
  row.names = FALSE
)


# ================================================================
# 43. CREATE RESULTS FOLDER
# ================================================================

dataset_directory <-
  dirname(
    file_path
  )


results_folder <-
  file.path(
    dataset_directory,
    "Results"
  )


if (
  !dir.exists(
    results_folder
  )
) {
  
  dir.create(
    results_folder,
    recursive = TRUE
  )
  
  
  cat(
    "\nResults folder created:\n",
    results_folder,
    "\n"
  )
  
} else {
  
  cat(
    "\nResults folder already exists:\n",
    results_folder,
    "\n"
  )
  
}


# ================================================================
# 44. DATASET NAME
# ================================================================

dataset_name <-
  tools::file_path_sans_ext(
    basename(
      file_path
    )
  )


# ================================================================
# 45. OUTPUT FILE NAMES
# ================================================================

all_results_file <-
  file.path(
    
    results_folder,
    
    paste0(
      dataset_name,
      "_k",
      k,
      "_N",
      N,
      "_all_runs.csv"
    )
    
  )


average_file <-
  file.path(
    
    results_folder,
    
    paste0(
      dataset_name,
      "_k",
      k,
      "_N",
      N,
      "_average.csv"
    )
    
  )


sd_file <-
  file.path(
    
    results_folder,
    
    paste0(
      dataset_name,
      "_k",
      k,
      "_N",
      N,
      "_SD.csv"
    )
    
  )


mean_sd_file <-
  file.path(
    
    results_folder,
    
    paste0(
      dataset_name,
      "_k",
      k,
      "_N",
      N,
      "_Mean_SD.csv"
    )
    
  )


# ================================================================
# 46. SAVE ALL RUN RESULTS
# ================================================================

write.csv(
  raw_results,
  all_results_file,
  row.names = FALSE
)


# ================================================================
# 47. SAVE AVERAGE
# ================================================================

write.csv(
  average_results,
  average_file,
  row.names = FALSE
)


# ================================================================
# 48. SAVE SD
# ================================================================

write.csv(
  sd_results,
  sd_file,
  row.names = FALSE
)


# ================================================================
# 49. SAVE MEAN +/- SD
# ================================================================

write.csv(
  mean_sd_results,
  mean_sd_file,
  row.names = FALSE
)


# ================================================================
# 50. FILE VERIFICATION
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "RESULTS SAVED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nAll runs file :",
  file.exists(
    all_results_file
  ),
  "\n"
)

cat(
  "Average file  :",
  file.exists(
    average_file
  ),
  "\n"
)

cat(
  "SD file       :",
  file.exists(
    sd_file
  ),
  "\n"
)

cat(
  "Mean-SD file  :",
  file.exists(
    mean_sd_file
  ),
  "\n"
)


# ================================================================
# 51. FINAL SUMMARY
# ================================================================

cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT COMPLETED SUCCESSFULLY\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nDataset        :",
  basename(file_path),
  "\n"
)


cat(
  "Clusters (k)   :",
  k,
  "\n"
)


cat(
  "Number of runs :",
  N,
  "\n"
)


cat(
  "Rand Index     :",
  ifelse(
    calculate_RI,
    "Calculated",
    "Not calculated"
  ),
  "\n"
)


if (
  calculate_RI
) {
  
  cat(
    "Ground truth   :",
    truth_name,
    "\n"
  )
  
}


cat(
  "\nMETHOD ORDER:\n"
)


for (
  i in seq_along(
    method_display_names
  )
) {
  
  cat(
    i,
    ".",
    method_display_names[i],
    "\n"
  )
  
}


cat(
  "\nCOSOHUC:\n"
)

cat(
  "Deterministic = YES\n"
)

cat(
  "Calculated once = YES\n"
)

cat(
  "Same result copied to all runs = YES\n"
)

cat(
  "COSOHUC uses its own original kmeans() call = YES\n"
)

cat(
  "ECKM uses published empty-circle LEC initialization logic = YES\n"
)

cat(
  "ECKM is deterministic and calculated once = YES\n"
)

cat(
  "Mini-Batch K-Means++ uses paper-based batch initialization = YES\n"
)

cat(
  "Mini-Batch K-Means++ is stochastic and re-run with fixed seeds = YES\n"
)


cat(
  "\nResults saved in:\n",
  results_folder,
  "\n"
)


cat(
  "\n============================================================\n"
)
