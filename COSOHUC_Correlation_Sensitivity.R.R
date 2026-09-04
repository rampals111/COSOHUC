# ========================================================================
# COSOHUC + 10-METHOD FEATURE-CORRELATION SENSITIVITY EXPERIMENT
# ========================================================================
#
# PURPOSE
# ------------------------------------------------------------------------
# Evaluate the same 10 clustering methods used in the main G-RUN experiment
# over controlled feature-correlation levels.
#
# METHODS
# ------------------------------------------------------------------------
# 1. Farthest First
# 2. Canopy
# 3. K-Means++
# 4. Hartigan-Wong
# 5. MacQueen
# 6. Lloyd
# 7. Forgy
# 8. ECKM
# 9. Mini-Batch K-Means++
# 10. COSOHUC
#
#
# CANOPY REPEATED-RUN IMPLEMENTATION
# ------------------------------------------------------------------------
# For each fixed correlation dataset, the pairwise distance matrix and
# Canopy thresholds T1/T2 are computed once because the dataset is unchanged
# across repeated runs. The Canopy initialization/fallback and final
# clustering are still evaluated within every run using that precomputed
# information.
#
# EXPERIMENT DESIGN
# ------------------------------------------------------------------------
# For EACH requested correlation level:
#
#   - ONE fixed synthetic dataset is generated.
#   - The SAME dataset is used for all N runs.
#   - Stochastic methods use a different fixed seed per run.
#   - COSOHUC and ECKM are computed ONCE for that fixed dataset and copied
#     to all runs.
#
# USER CAN SELECT:
#   - observations
#   - number of features
#   - k (any valid integer)
#   - number of runs
#   - mini-batch size
#   - correlation levels
#
# IMPORTANT:
# The requested rho controls the covariance/correlation structure of the
# synthetic noise. The program ALSO reports the realized mean within-cluster
# feature correlation so the experiment can be verified numerically.
#
# ========================================================================


# ========================================================================
# 1. REQUIRED PACKAGE
# ========================================================================

if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}

if (!requireNamespace("deldir", quietly = TRUE)) {
  install.packages("deldir")
}

library(cluster)
library(deldir)



# ========================================================================
# 2. USER ENTERS EXPERIMENT SETTINGS
# ========================================================================

repeat {

  input_value <-
    trimws(
      readline(
        prompt =
          "\nEnter number of observations/rows: "
      )
    )

  N_OBSERVATIONS <-
    suppressWarnings(
      as.integer(input_value)
    )

  if (
    !is.na(N_OBSERVATIONS) &&
    N_OBSERVATIONS >= 6
  ) {
    break
  }

  cat(
    "\nINVALID VALUE. Enter an integer >= 6.\n"
  )
}


repeat {

  input_value <-
    trimws(
      readline(
        prompt =
          "\nEnter number of features/columns: "
      )
    )

  N_FEATURES <-
    suppressWarnings(
      as.integer(input_value)
    )

  if (
    !is.na(N_FEATURES) &&
    N_FEATURES >= 2
  ) {
    break
  }

  cat(
    "\nINVALID VALUE. Enter an integer >= 2.\n"
  )
}


repeat {

  input_value <-
    trimws(
      readline(
        prompt =
          paste0(
            "\nEnter number of clusters (k) [2-",
            N_OBSERVATIONS - 1,
            "]: "
          )
      )
    )

  K <-
    suppressWarnings(
      as.integer(input_value)
    )

  if (
    !is.na(K) &&
    K >= 2 &&
    K < N_OBSERVATIONS
  ) {

    # ECKM and standard k-means need k distinct centers. A dataset with
    # fewer than k unique rows cannot support k distinct initial centers.
    # Synthetic continuous data will normally have N_OBSERVATIONS unique rows.
    break
  }

  cat(
    "\nINVALID k.\n"
  )
}


repeat {

  input_value <-
    trimws(
      readline(
        prompt =
          "\nEnter Mini-Batch size (b): "
      )
    )

  b <-
    suppressWarnings(
      as.integer(input_value)
    )

  if (
    !is.na(b) &&
    b >= 1 &&
    b <= N_OBSERVATIONS
  ) {
    break
  }

  cat(
    "\nINVALID batch size.\n"
  )
}


repeat {

  input_value <-
    trimws(
      readline(
        prompt =
          "\nEnter number of repeated runs (N): "
      )
    )

  N <-
    suppressWarnings(
      as.integer(input_value)
    )

  if (
    !is.na(N) &&
    N >= 1
  ) {
    break
  }

  cat(
    "\nINVALID N. Enter an integer >= 1.\n"
  )
}


cat(
  "\n============================================================\n"
)

cat(
  "CORRELATION LEVEL SELECTION\n"
)

cat(
  "============================================================\n"
)

cat(
  "\n1 = Default levels: 0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95\n"
)

cat(
  "2 = Enter custom levels\n"
)


repeat {

  choice <-
    trimws(
      readline(
        prompt =
          "\nEnter 1 or 2: "
      )
    )

  if (
    choice %in% c("1", "2")
  ) {
    break
  }

  cat(
    "\nINVALID choice.\n"
  )
}


if (
  choice == "1"
) {

  CORRELATION_LEVELS <- c(
    0.00,
    0.20,
    0.40,
    0.60,
    0.80,
    0.90,
    0.95
  )

} else {

  repeat {

    correlation_input <-
      trimws(
        readline(
          prompt =
            "\nEnter correlation levels separated by commas: "
        )
      )

    values <-
      suppressWarnings(
        as.numeric(
          trimws(
            unlist(
              strsplit(
                correlation_input,
                ","
              )
            )
          )
        )
      )

    if (
      length(values) == 0 ||
      anyNA(values) ||
      any(values < 0 | values >= 1)
    ) {

      cat(
        "\nINVALID levels. Each value must be >= 0 and < 1.\n"
      )

      next
    }

    CORRELATION_LEVELS <-
      sort(
        unique(
          values
        )
      )

    break
  }
}


B <-
  ceiling(
    N_OBSERVATIONS / b
  )


cat(
  "\n============================================================\n"
)

cat(
  "FINAL SETTINGS\n"
)

cat(
  "============================================================\n"
)

cat(
  "Observations :",
  N_OBSERVATIONS,
  "\n"
)

cat(
  "Features     :",
  N_FEATURES,
  "\n"
)

cat(
  "k            :",
  K,
  "\n"
)

cat(
  "Runs         :",
  N,
  "\n"
)

cat(
  "Batch size   :",
  b,
  "\n"
)

cat(
  "Batches      :",
  B,
  "\n"
)

cat(
  "Correlation  :",
  paste(
    CORRELATION_LEVELS,
    collapse = ", "
  ),
  "\n"
)


# ========================================================================
# 3. CORRELATION DATA GENERATION
# ========================================================================

create_cluster_means <- function(
    k,
    p,
    separation = 6
) {

  if (
    k < 2 ||
    p < 2
  ) {

    stop(
      "k must be >= 2 and p must be >= 2."
    )
  }

  # --------------------------------------------------------------
  # All features carry cluster information.
  #
  # The cluster means are evenly spaced from negative to positive
  # values. This gives a controlled, interpretable separation while
  # allowing arbitrary k.
  # --------------------------------------------------------------

  cluster_axis <-
    seq(
      -separation,
      separation,
      length.out = k
    )

  means <-
    matrix(
      0,
      nrow = k,
      ncol = p
    )

  for (
    cl in seq_len(k)
  ) {

    means[
      cl,
      1:p
    ] <-
      cluster_axis[
        cl
      ]
  }

  colnames(means) <-
    paste0(
      "Feature_",
      seq_len(p)
    )

  means
}


create_correlation_matrix <- function(
    p,
    rho
) {

  if (
    rho < 0 ||
    rho >= 1
  ) {

    stop(
      "rho must satisfy 0 <= rho < 1."
    )
  }

  Sigma <-
    matrix(
      rho,
      nrow =
        p,
      ncol =
        p
    )

  diag(
    Sigma
  ) <-
    1

  Sigma
}


generate_synthetic_dataset <- function(
    n,
    p,
    k,
    rho,
    seed,
    separation = 6
) {

  set.seed(
    seed
  )

  if (
    n < k
  ) {

    stop(
      "Number of observations must be >= k."
    )
  }

  if (
    p < 2
  ) {

    stop(
      "Number of features must be >= 2."
    )
  }

  # --------------------------------------------------------------
  # Balanced cluster sizes.
  # --------------------------------------------------------------

  cluster_sizes <-
    rep(
      floor(
        n / k
      ),
      k
    )

  remainder <-
    n %% k

  if (
    remainder > 0
  ) {

    cluster_sizes[
      seq_len(
        remainder
      )
    ] <-
      cluster_sizes[
        seq_len(
          remainder
        )
      ] +
      1
  }

  # --------------------------------------------------------------
  # Cluster mean structure.
  # --------------------------------------------------------------

  means <-
    create_cluster_means(
      k =
        k,
      p =
        p,
      separation =
        separation
    )

  # --------------------------------------------------------------
  # Feature-correlation structure.
  # --------------------------------------------------------------

  Sigma <-
    create_correlation_matrix(
      p =
        p,
      rho =
        rho
    )

  chol_Sigma <-
    chol(
      Sigma
    )

  X_list <-
    vector(
      "list",
      k
    )

  true_labels <-
    integer(
      0
    )

  for (
    cl in seq_len(
      k
    )
  ) {

    n_cl <-
      cluster_sizes[
        cl
      ]

    Z <-
      matrix(
        rnorm(
          n_cl *
          p
        ),
        nrow =
          n_cl,
        ncol =
          p
      )

    noise <-
      Z %*%
      chol_Sigma

    cluster_data <-
      sweep(
        noise,
        2,
        means[
          cl,
          ],
        "+"
      )

    X_list[[cl]] <-
      cluster_data

    true_labels <-
      c(
        true_labels,
        rep(
          cl,
          n_cl
        )
      )
  }

  X <-
    do.call(
      rbind,
      X_list
    )

  colnames(X) <-
    paste0(
      "Feature_",
      seq_len(p)
    )

  list(
    X =
      X,

    true_labels =
      true_labels
  )
}


calculate_within_cluster_correlation <- function(
    X,
    true_labels
) {

  values <-
    numeric(
      0
    )

  for (
    cl in sort(
      unique(
        true_labels
      )
    )
  ) {

    cluster_data <-
      X[
        true_labels == cl,
        ,
        drop = FALSE
      ]

    if (
      nrow(cluster_data) < 2
    ) {
      next
    }

    C <-
      cor(
        cluster_data
      )

    upper <-
      upper.tri(
        C
      )

    values <-
      c(
        values,
        C[
          upper
        ]
      )
  }

  if (
    length(values) == 0
  ) {
    return(NA_real_)
  }

  mean(
    values,
    na.rm =
      TRUE
  )
}


# ========================================================================
# 4. METHOD ORDER
# ========================================================================

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


run_COSOHUC <- function(
    X,
    k
) {

  initial_centers <-
    COSOHUC_initial_centers(
      X,
      k
    )

  if (
    nrow(
      unique(
        initial_centers
      )
    ) != k
  ) {

    stop(
      "COSOHUC produced non-distinct initial center vectors."
    )
  }

  result <-
    tryCatch(

      kmeans(
        X,
        centers =
          initial_centers,
        iter.max =
          100,
        nstart =
          1,
        algorithm =
          "Lloyd"
      ),

      warning = function(w) {

        # Convert empty-cluster warning into an explicit diagnostic error.
        message_text <-
          conditionMessage(
            w
          )

        if (
          grepl(
            "empty cluster",
            message_text,
            ignore.case =
              TRUE
          )
        ) {

          invokeRestart(
            "muffleWarning"
          )

          stop(
            paste0(
              "COSOHUC final K-means produced an empty cluster: ",
              message_text
            )
          )
        }

        invokeRestart(
          "muffleWarning"
        )

        warning(
          message_text
        )
      },

      error = function(e) {
        stop(
          conditionMessage(e)
        )
      }
    )

  cluster_sizes <-
    table(
      result$cluster
    )

  if (
    length(
      cluster_sizes
    ) != k ||
    any(
      cluster_sizes <= 0
    )
  ) {

    stop(
      paste0(
        "COSOHUC final clustering contains fewer than ",
        k,
        " non-empty clusters."
      )
    )
  }

  result$initial.centers <-
    initial_centers

  result
}


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

# ========================================================================
# ROBUST RANDOM-INITIALIZATION K-MEANS
# ========================================================================
#
# Some synthetic datasets/initial seeds can make base R's random
# initialization produce:
#
#   "empty cluster: try a better set of initial centers"
#
# This wrapper retries the SAME requested k-means algorithm with a
# deterministic sequence of seeds until a valid non-empty-cluster result
# is obtained. The clustering data are never modified and no jitter is
# added.
#
# The first attempt uses the requested seed. Later attempts use seed+1,
# seed+2, ... only when the previous attempt generated an empty-cluster
# warning or an error.
#
# ========================================================================

safe_random_kmeans <- function(
    X,
    k,
    algorithm,
    seed,
    iter.max = 100,
    max_attempts = 50
) {

  last_message <- NULL

  for (
    attempt in seq_len(max_attempts)
  ) {

    current_seed <-
      seed +
      attempt -
      1L

    warning_message <- NULL

    result <- tryCatch(

      withCallingHandlers(

        {

          set.seed(
            current_seed
          )

          kmeans(
            X,
            centers =
              k,
            iter.max =
              iter.max,
            nstart =
              1,
            algorithm =
              algorithm
          )

        },

        warning = function(w) {

          warning_message <<-
            conditionMessage(
              w
            )

          invokeRestart(
            "muffleWarning"
          )
        }

      ),

      error = function(e) {

        last_message <<-
          conditionMessage(
            e
          )

        NULL

      }

    )

    if (
      !is.null(result)
    ) {

      cluster_sizes <-
        table(
          result$cluster
        )

      # A valid clustering must contain all k clusters.
      if (
        length(cluster_sizes) == k &&
        all(cluster_sizes > 0)
      ) {

        result$initialization_seed <-
          current_seed

        result$attempt <-
          attempt

        result$empty_cluster_warning <-
          if (
            !is.null(warning_message)
          ) {
            warning_message
          } else {
            NA_character_
          }

        return(
          result
        )
      }

      last_message <-
        "The attempted clustering contained an empty cluster."
    }
  }

  stop(
    paste0(
      "k-means failed to obtain ",
      k,
      " non-empty clusters after ",
      max_attempts,
      " deterministic attempts. Last message: ",
      last_message
    )
  )
}


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

canopy_precompute <- function(
    X,
    D = NULL
) {

  n <- nrow(X)

  if (
    is.null(D)
  ) {

    D <-
      as.matrix(
        dist(
          X
        )
      )
  }

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
      "Cannot calculate Canopy distances because all observations are identical."
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

  list(
    D = D,
    T1 = T1,
    T2 = T2,
    n = n
  )
}


canopy_initial_precomputed <- function(
    X,
    k,
    seed,
    canopy_terms
) {

  set.seed(
    seed
  )

  remaining <-
    seq_len(
      canopy_terms$n
    )

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

    selected <-
      c(
        selected,
        center_id
      )

    distances <-
      canopy_terms$D[
        center_id,
        remaining
      ]

    remove_ids <-
      remaining[
        distances <=
          canopy_terms$T2
      ]

    remaining <-
      setdiff(
        remaining,
        remove_ids
      )
  }

  # Same fallback principle as the existing Canopy implementation.
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
      i in seq_len(
        nrow(ff)
      )
    ) {

      if (
        length(selected) >= k
      ) {
        break
      }

      distances <-
        apply(
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

        selected <-
          c(
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
      paste0(
        "Canopy could not generate ",
        k,
        " initial centers."
      )
    )
  }

  selected <-
    selected[
      seq_len(k)
    ]

  centers <-
    X[
      selected,
      ,
      drop = FALSE
    ]

  if (
    nrow(
      unique(
        centers
      )
    ) < k
  ) {

    centers <-
      farthest_first_centers(
        X,
        k,
        seed
      )
  }

  if (
    nrow(
      unique(
        centers
      )
    ) < k
  ) {

    stop(
      paste0(
        "Canopy could not generate ",
        k,
        " distinct initial center vectors."
      )
    )
  }

  colnames(centers) <-
    colnames(X)

  centers
}


# Backward-compatible wrapper. This retains the same public function name.
canopy_initial <- function(
    X,
    k,
    seed
) {

  canopy_terms <-
    canopy_precompute(
      X
    )

  canopy_initial_precomputed(
    X =
      X,
    k =
      k,
    seed =
      seed,
    canopy_terms =
      canopy_terms
  )
}


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

calculate_dunn <- function(
    X = NULL,
    cluster,
    D = NULL
) {

  if (
    is.null(D)
  ) {

    if (
      is.null(X)
    ) {

      stop(
        "Either X or precomputed distance matrix D must be supplied."
      )
    }

    D <-
      as.matrix(
        dist(
          X
        )
      )
  }

  clusters <-
    sort(
      unique(
        cluster
      )
    )

  if (
    length(clusters) < 2
  ) {

    return(
      NA_real_
    )
  }

  maximum_intra <-
    0

  for (
    cl in clusters
  ) {

    members <-
      which(
        cluster == cl
      )

    if (
      length(members) > 1
    ) {

      intra <-
        D[
          members,
          members,
          drop = FALSE
        ]

      maximum_intra <-
        max(
          maximum_intra,
          max(
            intra
          )
        )
    }
  }

  minimum_inter <-
    Inf

  for (
    a in seq_len(
      length(clusters) - 1
    )
  ) {

    for (
      b in (a + 1):
        length(clusters)
    ) {

      members_a <-
        which(
          cluster ==
            clusters[a]
        )

      members_b <-
        which(
          cluster ==
            clusters[b]
        )

      inter <-
        D[
          members_a,
          members_b,
          drop = FALSE
        ]

      minimum_inter <-
        min(
          minimum_inter,
          min(
            inter
          )
        )
    }
  }

  if (
    maximum_intra <= 0 ||
    !is.finite(
      minimum_inter
    )
  ) {

    return(
      NA_real_
    )
  }

  minimum_inter /
    maximum_intra
}


calculate_silhouette <- function(
    X = NULL,
    cluster,
    D = NULL
) {

  tryCatch(

    {

      if (
        is.null(D)
      ) {

        if (
          is.null(X)
        ) {

          stop(
            "Either X or precomputed distance matrix D must be supplied."
          )
        }

        D <-
          as.matrix(
            dist(
              X
            )
          )
      }

      sil <-
        cluster::silhouette(
          as.integer(
            cluster
          ),
          D
        )

      mean(
        sil[
          ,
          "sil_width"
        ],
        na.rm =
          TRUE
      )

    },

    error = function(e) {

      NA_real_

    }
  )
}


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


calculate_all_indices <- function(
    X,
    result,
    D = NULL
) {

  if (
    is.null(result) ||
    is.null(result$cluster)
  ) {

    stop(
      "Clustering result does not contain cluster assignments."
    )
  }

  cluster <-
    as.integer(
      result$cluster
    )

  if (
    length(cluster) !=
      nrow(X) ||
    anyNA(cluster)
  ) {

    stop(
      "Invalid cluster assignments returned by the clustering method."
    )
  }

  if (
    is.null(D)
  ) {

    D <-
      as.matrix(
        dist(
          X
        )
      )
  }

  TWSS <-
    if (
      !is.null(
        result$tot.withinss
      ) &&
      length(
        result$tot.withinss
      ) == 1 &&
      is.finite(
        result$tot.withinss
      )
    ) {

      as.numeric(
        result$tot.withinss
      )

    } else {

      0
    }

  if (
    TWSS == 0
  ) {

    for (
      cl in sort(
        unique(
          cluster
        )
      )
    ) {

      members <-
        which(
          cluster == cl
        )

      cluster_data <-
        X[
          members,
          ,
          drop = FALSE
        ]

      center <-
        colMeans(
          cluster_data
        )

      TWSS <-
        TWSS +
        sum(
          rowSums(
            sweep(
              cluster_data,
              2,
              center,
              "-"
            )^2
          )
        )
    }
  }

  grand_mean <-
    colMeans(
      X
    )

  TotalSS <-
    sum(
      rowSums(
        sweep(
          X,
          2,
          grand_mean,
          "-"
        )^2
      )
    )

  BSS <-
    TotalSS -
    TWSS

  if (
    abs(BSS) <
      1e-10
  ) {

    BSS <- 0
  }

  BSS_Percentage <-
    if (
      is.finite(
        TotalSS
      ) &&
      TotalSS > 0
    ) {

      value <-
        100 *
        BSS /
        TotalSS

      if (
        abs(value) <
          1e-10
      ) {

        value <-
          0
      }

      value

    } else {

      NA_real_
    }

  list(

    Iteration =
      if (
        !is.null(
          result$iter
        ) &&
        length(
          result$iter
        ) == 1 &&
        is.finite(
          result$iter
        )
      ) {

        as.numeric(
          result$iter
        )

      } else {

        NA_real_
      },

    BSS =
      as.numeric(
        BSS
      ),

    TWSS =
      as.numeric(
        TWSS
      ),

    TotalSS =
      as.numeric(
        TotalSS
      ),

    BSS_Percentage =
      as.numeric(
        BSS_Percentage
      ),

    Dunn =
      as.numeric(
        calculate_dunn(
          cluster =
            cluster,
          D =
            D
        )
      ),

    DBI =
      as.numeric(
        calculate_davies_bouldin(
          X,
          cluster
        )
      ),

    Silhouette =
      as.numeric(
        calculate_silhouette(
          cluster =
            cluster,
          D =
            D
        )
      )
  )
}





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


# ========================================================================
# 5. MAIN CORRELATION-SENSITIVITY EXPERIMENT
# ========================================================================

raw_results <-
  data.frame(

    Correlation_Level =
      numeric(),

    Actual_Within_Cluster_Correlation =
      numeric(),

    Run =
      integer(),

    Method =
      character(),

    Iteration =
      numeric(),

    BSS =
      numeric(),

    TWSS =
      numeric(),

    TotalSS =
      numeric(),

    BSS_Percentage =
      numeric(),

    Dunn =
      numeric(),

    DBI =
      numeric(),

    Silhouette =
      numeric(),

    Rand_Index =
      numeric(),

    stringsAsFactors =
      FALSE
  )


add_result_row <- function(
    raw_results,
    rho,
    actual_rho,
    run_number,
    method_display,
    result,
    X,
    true_labels,
    D = NULL
) {

  metrics <-
    calculate_all_indices(
      X,
      result,
      D =
        D
    )

  RI <-
    calculate_rand_index(
      true_labels,
      result$cluster
    )

  rbind(
    raw_results,
    data.frame(

      Correlation_Level =
        rho,

      Actual_Within_Cluster_Correlation =
        actual_rho,

      Run =
        run_number,

      Method =
        method_display,

      Iteration =
        metrics$Iteration,

      BSS =
        metrics$BSS,

      TWSS =
        metrics$TWSS,

      TotalSS =
        metrics$TotalSS,

      BSS_Percentage =
        metrics$BSS_Percentage,

      Dunn =
        metrics$Dunn,

      DBI =
        metrics$DBI,

      Silhouette =
        metrics$Silhouette,

      Rand_Index =
        RI,

      stringsAsFactors =
        FALSE
    )
  )
}


for (
  rho_index in seq_along(
    CORRELATION_LEVELS
  )
) {

  rho <-
    CORRELATION_LEVELS[
      rho_index
    ]

  fixed_data <-
    generate_synthetic_dataset(

      n =
        N_OBSERVATIONS,

      p =
        N_FEATURES,

      k =
        K,

      rho =
        rho,

      seed =
        50000L +
        rho_index
    )

  X_fixed <-
    fixed_data$X

  true_labels_fixed <-
    fixed_data$true_labels

  actual_rho <-
    calculate_within_cluster_correlation(
      X_fixed,
      true_labels_fixed
    )

  # ==============================================================
  # FIXED DISTANCE MATRIX FOR THIS CORRELATION LEVEL
  # ==============================================================
  #
  # X_fixed is unchanged throughout all N repeated runs for this
  # correlation level. Therefore pairwise distances are calculated
  # exactly ONCE and reused by Canopy, Dunn and Silhouette.
  #
  # This is a performance optimization only; it does not change
  # any clustering algorithm or metric definition.
  #
  # ==============================================================

  D_fixed <-
    as.matrix(
      dist(
        X_fixed
      )
    )

  # Verify that all requested clusters are represented in the generated data.
  generated_class_count <-
    length(
      unique(
        true_labels_fixed
      )
    )

  if (
    generated_class_count != K
  ) {

    stop(
      paste0(
        "Generated dataset contains ",
        generated_class_count,
        " true classes, but k = ",
        K,
        "."
      )
    )
  }

  cat(
    "\n============================================================\n"
  )

  cat(
    "CORRELATION LEVEL: ",
    sprintf(
      "%.2f",
      rho
    ),
    "\n",
    sep = ""
  )

  cat(
    "REALIZED WITHIN-CLUSTER CORRELATION: ",
    sprintf(
      "%.6f",
      actual_rho
    ),
    "\n",
    sep = ""
  )

  cat(
    "============================================================\n"
  )


  # --------------------------------------------------------------
  # Canopy dataset-dependent quantities: compute ONCE per correlation
  # level because X_fixed does not change across the 40 runs.
  # --------------------------------------------------------------

  canopy_terms <-
    canopy_precompute(
      X =
        X_fixed,
      D =
        D_fixed
    )

  # --------------------------------------------------------------
  # Deterministic COSOHUC
  # --------------------------------------------------------------

  COSOHUC_result <-
    run_COSOHUC(
      X_fixed,
      K
    )

  COSOHUC_metrics <-
    calculate_all_indices(
      X_fixed,
      COSOHUC_result,
      D =
        D_fixed
    )

  COSOHUC_RI <-
    calculate_rand_index(
      true_labels_fixed,
      COSOHUC_result$cluster
    )

  # --------------------------------------------------------------
  # Deterministic ECKM
  # --------------------------------------------------------------

  ECKM_result <-
    run_ECKM(
      X_fixed,
      K
    )

  ECKM_metrics <-
    calculate_all_indices(
      X_fixed,
      ECKM_result,
      D =
        D_fixed
    )

  ECKM_RI <-
    calculate_rand_index(
      true_labels_fixed,
      ECKM_result$cluster
    )


  for (
    run_number in seq_len(N)
  ) {

    seed <-
      100000L +
      rho_index * 10000L +
      run_number

    cat(
      "Run ",
      run_number,
      " / ",
      N,
      "\n",
      sep = ""
    )


    # ============================================================
    # FARTHest FIRST
    # ============================================================

    farthest_initial <-
      farthest_first_centers(
        X_fixed,
        K,
        seed
      )

    result_farthest <-
      run_from_centers(
        X_fixed,
        farthest_initial
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "Farthest First",
        result_farthest,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # CANOPY
    # ============================================================

    canopy_initial_centers <-
      canopy_initial_precomputed(
        X =
          X_fixed,
        k =
          K,
        seed =
          seed,
        canopy_terms =
          canopy_terms
      )

    result_canopy <-
      run_from_centers(
        X_fixed,
        canopy_initial_centers
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "Canopy",
        result_canopy,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # K-MEANS++
    # ============================================================

    kmeanspp_initial <-
      kmeans_plus_plus_centers(
        X_fixed,
        K,
        seed
      )

    result_kmeanspp <-
      run_from_centers(
        X_fixed,
        kmeanspp_initial
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "K-Means++",
        result_kmeanspp,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # HARTIGAN-WONG
    # ============================================================

    result_HW <-
      safe_random_kmeans(
        X =
          X_fixed,
        k =
          K,
        algorithm =
          "Hartigan-Wong",
        seed =
          seed,
        iter.max =
          100
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "Hartigan-Wong",
        result_HW,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # MACQUEEN
    # ============================================================

    result_MacQueen <-
      safe_random_kmeans(
        X =
          X_fixed,
        k =
          K,
        algorithm =
          "MacQueen",
        seed =
          seed,
        iter.max =
          100
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "MacQueen",
        result_MacQueen,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # LLOYD
    # ============================================================

    result_Lloyd <-
      safe_random_kmeans(
        X =
          X_fixed,
        k =
          K,
        algorithm =
          "Lloyd",
        seed =
          seed,
        iter.max =
          100
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "Lloyd",
        result_Lloyd,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # FORGY
    # ============================================================

    set.seed(
      seed
    )

    unique_X <-
      unique(
        X_fixed
      )

    if (
      nrow(unique_X) < K
    ) {

      stop(
        paste0(
          "Forgy cannot select ",
          K,
          " distinct centers: only ",
          nrow(unique_X),
          " unique observations exist."
        )
      )
    }

    forgy_indices <-
      sample(
        seq_len(
          nrow(unique_X)
        ),
        size =
          K,
        replace =
          FALSE
      )

    forgy_initial <-
      unique_X[
        forgy_indices,
        ,
        drop =
          FALSE
      ]

    # --------------------------------------------------------------
    # Run Lloyd from the Forgy-selected centers.
    # Retry Forgy itself if the selected centers yield an empty cluster.
    # --------------------------------------------------------------

    forgy_result <- NULL

    for (
      forgy_attempt in seq_len(50)
    ) {

      if (
        forgy_attempt > 1
      ) {

        set.seed(
          seed +
            forgy_attempt -
            1L
        )

        forgy_indices <-
          sample(
            seq_len(
              nrow(unique_X)
            ),
            size =
              K,
            replace =
              FALSE
          )

        forgy_initial <-
          unique_X[
            forgy_indices,
            ,
            drop =
              FALSE
          ]
      }

      warning_message <- NULL

      candidate <- tryCatch(

        withCallingHandlers(

          kmeans(
            X_fixed,
            centers =
              forgy_initial,
            iter.max =
              100,
            nstart =
              1,
            algorithm =
              "Lloyd"
          ),

          warning = function(w) {

            warning_message <<-
              conditionMessage(w)

            invokeRestart(
              "muffleWarning"
            )
          }

        ),

        error = function(e) {
          NULL
        }
      )

      if (
        !is.null(candidate) &&
        length(
          unique(
            candidate$cluster
          )
        ) == K &&
        all(
          table(candidate$cluster) > 0
        )
      ) {

        forgy_result <-
          candidate

        break
      }
    }

    if (
      is.null(forgy_result)
    ) {

      stop(
        paste0(
          "Forgy failed to obtain ",
          K,
          " non-empty clusters after 50 deterministic attempts."
        )
      )
    }

    result_Forgy <-
      forgy_result

    result_Forgy$initial.centers <-
      forgy_initial

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "Forgy",
        result_Forgy,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # ECKM - REUSE DETERMINISTIC RESULT
    # ============================================================

    raw_results <-
      rbind(
        raw_results,
        data.frame(

          Correlation_Level =
            rho,

          Actual_Within_Cluster_Correlation =
            actual_rho,

          Run =
            run_number,

          Method =
            "ECKM",

          Iteration =
            ECKM_metrics$Iteration,

          BSS =
            ECKM_metrics$BSS,

          TWSS =
            ECKM_metrics$TWSS,

          TotalSS =
            ECKM_metrics$TotalSS,

          BSS_Percentage =
            ECKM_metrics$BSS_Percentage,

          Dunn =
            ECKM_metrics$Dunn,

          DBI =
            ECKM_metrics$DBI,

          Silhouette =
            ECKM_metrics$Silhouette,

          Rand_Index =
            ECKM_RI,

          stringsAsFactors =
            FALSE
        )
      )


    # ============================================================
    # MINI-BATCH K-MEANS++
    # ============================================================

    MBKMPP_result <-
      run_mini_batch_kmeans_plus_plus(
        X_fixed,
        K,
        b,
        seed
      )

    raw_results <-
      add_result_row(
        raw_results,
        rho,
        actual_rho,
        run_number,
        "Mini-Batch K-Means++",
        MBKMPP_result,
        X_fixed,
        true_labels_fixed,
        D = D_fixed
      )


    # ============================================================
    # COSOHUC - REUSE DETERMINISTIC RESULT
    # ============================================================

    raw_results <-
      rbind(
        raw_results,
        data.frame(

          Correlation_Level =
            rho,

          Actual_Within_Cluster_Correlation =
            actual_rho,

          Run =
            run_number,

          Method =
            "COSOHUC",

          Iteration =
            COSOHUC_metrics$Iteration,

          BSS =
            COSOHUC_metrics$BSS,

          TWSS =
            COSOHUC_metrics$TWSS,

          TotalSS =
            COSOHUC_metrics$TotalSS,

          BSS_Percentage =
            COSOHUC_metrics$BSS_Percentage,

          Dunn =
            COSOHUC_metrics$Dunn,

          DBI =
            COSOHUC_metrics$DBI,

          Silhouette =
            COSOHUC_metrics$Silhouette,

          Rand_Index =
            COSOHUC_RI,

          stringsAsFactors =
            FALSE
        )
      )
  }
}


# ========================================================================
# 6. FINAL ORDERING
# ========================================================================

raw_results$Method <-
  factor(
    raw_results$Method,
    levels =
      method_display_names
  )

raw_results <-
  raw_results[
    order(
      raw_results$Correlation_Level,
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


# ========================================================================
# 7. AVERAGE
# ========================================================================

numeric_variables <- c(

  "Actual_Within_Cluster_Correlation",

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


average_results <-
  aggregate(

    raw_results[
      ,
      numeric_variables,
      drop = FALSE
    ],

    by =
      list(

        Correlation_Level =
          raw_results$Correlation_Level,

        Method =
          raw_results$Method

      ),

    FUN =
      function(x) {

        if (
          all(
            is.na(x)
          )
        ) {

          return(
            NA_real_
          )
        }

        mean(
          x,
          na.rm =
            TRUE
        )
      }
  )

colnames(
  average_results
)[1:2] <-
  c(
    "Correlation_Level",
    "Method"
  )


# ========================================================================
# 8. SD
# ========================================================================

sd_results <-
  aggregate(

    raw_results[
      ,
      numeric_variables,
      drop = FALSE
    ],

    by =
      list(

        Correlation_Level =
          raw_results$Correlation_Level,

        Method =
          raw_results$Method

      ),

    FUN =
      function(x) {

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
)[1:2] <-
  c(
    "Correlation_Level",
    "Method"
  )


# ========================================================================
# 9. METHOD ORDER FOR SUMMARY TABLES
# ========================================================================

average_results$Method <-
  factor(
    average_results$Method,
    levels =
      method_display_names
  )

average_results <-
  average_results[
    order(
      average_results$Correlation_Level,
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
      sd_results$Correlation_Level,
      sd_results$Method
    ),
    ,
    drop = FALSE
  ]

sd_results$Method <-
  as.character(
    sd_results$Method
  )


# ========================================================================
# 10. MEAN +/- SD
# ========================================================================

mean_sd_results <-
  data.frame(

    Correlation_Level =
      average_results$Correlation_Level,

    Method =
      average_results$Method,

    stringsAsFactors =
      FALSE
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

      sprintf(
        "%.6f",
        mean_values
      ),

      " +/- ",

      sprintf(
        "%.6f",
        sd_values
      )
    )
}


# ========================================================================
# 11. COSOHUC SUMMARY
# ========================================================================

COSOHUC_summary <-
  average_results[
    average_results$Method ==
      "COSOHUC",
    ,
    drop = FALSE
  ]


# ========================================================================
# 12. CORRELATION VERIFICATION
# ========================================================================

correlation_verification <-
  unique(

    raw_results[
      ,
      c(
        "Correlation_Level",
        "Actual_Within_Cluster_Correlation"
      ),
      drop = FALSE
    ]
  )


correlation_verification <-
  correlation_verification[
    order(
      correlation_verification$Correlation_Level
    ),
    ,
    drop = FALSE
  ]


# ========================================================================
# 13. DISPLAY
# ========================================================================

cat(
  "\n============================================================\n"
)

cat(
  "AVERAGE RESULTS\n"
)

cat(
  "============================================================\n\n"
)

print(
  average_results,
  row.names =
    FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "STANDARD DEVIATION RESULTS\n"
)

cat(
  "============================================================\n\n"
)

print(
  sd_results,
  row.names =
    FALSE
)


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
  row.names =
    FALSE
)


cat(
  "\n============================================================\n"
)

cat(
  "CORRELATION VERIFICATION\n"
)

cat(
  "============================================================\n\n"
)

print(
  correlation_verification,
  row.names =
    FALSE
)


# ========================================================================
# 14. OUTPUT FOLDER
# ========================================================================
#
# Requested root folder:
#
# D:\kmean\DATA SET\data\usefull dataset\dataset my paper\correlation_sensitivity
#
# IMPORTANT:
# Windows commonly gives "Permission denied" when a CSV with the same
# filename is already OPEN in Excel. To avoid that problem, every run is
# saved in a NEW timestamped subfolder.
#
# Example:
# correlation_sensitivity/
#   Run_20260901_183000/
#       Correlation_Sensitivity_Raw_Results.csv
#       Correlation_Sensitivity_Average.csv
#       ...
#
# ========================================================================

output_root <-
  "D:/kmean/DATA SET/data/usefull dataset/dataset my paper/correlation_sensitivity"


if (
  !dir.exists(
    output_root
  )
) {

  dir.create(
    output_root,
    recursive =
      TRUE,
    showWarnings =
      FALSE
  )
}


if (
  !dir.exists(
    output_root
  )
) {

  stop(
    paste0(
      "Cannot create/access output folder:\n",
      output_root,
      "\n\nCheck that the folder exists and you have write permission."
    )
  )
}


# --------------------------------------------------------------
# Check whether the root directory is writable.
# --------------------------------------------------------------

test_file <-
  file.path(
    output_root,
    ".write_test.tmp"
  )

write_test_ok <-
  tryCatch(

    {

      con <-
        file(
          test_file,
          open =
            "w"
        )

      writeLines(
        "write test",
        con
      )

      close(
        con
      )

      unlink(
        test_file
      )

      TRUE
    },

    error = function(e) {

      FALSE
    }
  )


if (
  !write_test_ok
) {

  stop(
    paste0(
      "Permission denied for the Results folder:\n",
      output_root,
      "\n\n",
      "Close RStudio/Excel files using this folder and verify that the ",
      "folder is writable."
    )
  )
}


# --------------------------------------------------------------
# Create a NEW run-specific folder.
# This prevents conflicts with old CSV files that may still be open.
# --------------------------------------------------------------

run_stamp <-
  format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )


output_directory <-
  file.path(
    output_root,
    paste0(
      "Run_",
      run_stamp,
      "_Rows",
      N_OBSERVATIONS,
      "_Features",
      N_FEATURES,
      "_K",
      K,
      "_Runs",
      N
    )
  )


dir.create(
  output_directory,
  recursive =
    TRUE,
  showWarnings =
    FALSE
)


if (
  !dir.exists(
    output_directory
  )
) {

  stop(
    paste0(
      "Could not create run-specific output folder:\n",
      output_directory
    )
  )
}


cat(
  "\nResults will be saved to:\n",
  output_directory,
  "\n"
)


# ========================================================================
# 15. OUTPUT FILES
# ========================================================================

raw_file <-
  file.path(
    output_directory,
    "Correlation_Sensitivity_Raw_Results.csv"
  )

average_file <-
  file.path(
    output_directory,
    "Correlation_Sensitivity_Average.csv"
  )

sd_file <-
  file.path(
    output_directory,
    "Correlation_Sensitivity_SD.csv"
  )

mean_sd_file <-
  file.path(
    output_directory,
    "Correlation_Sensitivity_Mean_SD.csv"
  )

verification_file <-
  file.path(
    output_directory,
    "Correlation_Verification.csv"
  )

cosohuc_summary_file <-
  file.path(
    output_directory,
    "COSOHUC_Correlation_Summary.csv"
  )


# ========================================================================
# 16. SAVE
# ========================================================================
#
# All files are written into a NEW run-specific folder, so an existing
# Excel file from an older run cannot block the current write.
#
# ========================================================================

safe_write_csv <- function(
    data,
    filename
) {

  result <-
    tryCatch(

      {

        write.csv(
          data,
          filename,
          row.names =
            FALSE
        )

        TRUE
      },

      error = function(e) {

        cat(
          "\nSAVE ERROR:\n",
          conditionMessage(e),
          "\n",
          sep = ""
        )

        FALSE
      }
    )

  if (
    !result
  ) {

    stop(
      paste0(
        "\nCould not save file:\n",
        filename,
        "\n\n",
        "Close any Excel/R application using this file and run again."
      )
    )
  }

  invisible(
    TRUE
  )
}


safe_write_csv(
  raw_results,
  raw_file
)

safe_write_csv(
  average_results,
  average_file
)

safe_write_csv(
  sd_results,
  sd_file
)

safe_write_csv(
  mean_sd_results,
  mean_sd_file
)

safe_write_csv(
  correlation_verification,
  verification_file
)

safe_write_csv(
  COSOHUC_summary,
  cosohuc_summary_file
)


cat(
  "\nAll result files saved successfully.\n"
)

# ========================================================================
# 17. FINAL CHECKS
# ========================================================================

expected_rows <-
  length(
    CORRELATION_LEVELS
  ) *
  N *
  length(
    method_names
  )

if (
  nrow(raw_results) !=
  expected_rows
) {

  stop(
    paste0(
      "ERROR: expected ",
      expected_rows,
      " raw result rows but obtained ",
      nrow(raw_results),
      "."
    )
  )
}


# Check method completeness.
for (
  method in method_display_names
) {

  if (
    !any(
      raw_results$Method ==
        method
    )
  ) {

    stop(
      paste0(
        "ERROR: missing method in output: ",
        method
      )
    )
  }
}


# Check COSOHUC determinism at every rho.
for (
  rho in CORRELATION_LEVELS
) {

  rows <-
    raw_results[
      raw_results$Correlation_Level ==
        rho &
        raw_results$Method ==
        "COSOHUC",
      ,
      drop = FALSE
    ]

  for (
    variable in c(
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
  ) {

    values <-
      rows[
        ,
        variable
      ]

    if (
      all(
        is.na(values)
      )
    ) {

      next
    }

    if (
      any(
        abs(
          values -
            values[1]
        ) >
        1e-12
      )
    ) {

      stop(
        paste0(
          "ERROR: COSOHUC is not deterministic for rho = ",
          rho,
          " in metric ",
          variable
        )
      )
    }
  }
}


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
  "Expected raw rows = ",
  expected_rows,
  "\n",
  sep = ""
)

cat(
  "Actual raw rows   = ",
  nrow(raw_results),
  "\n",
  sep = ""
)

cat(
  "\nCOSOHUC determinism check: PASSED\n"
)

cat(
  "Method completeness check: PASSED\n"
)

cat(
  "\nOutput directory:\n",
  output_directory,
  "\n"
)

cat(
  "\nFiles created:\n"
)

cat(
  "1. ",
  basename(raw_file),
  "\n",
  sep = ""
)

cat(
  "2. ",
  basename(average_file),
  "\n",
  sep = ""
)

cat(
  "3. ",
  basename(sd_file),
  "\n",
  sep = ""
)

cat(
  "4. ",
  basename(mean_sd_file),
  "\n",
  sep = ""
)

cat(
  "5. ",
  basename(verification_file),
  "\n",
  sep = ""
)

cat(
  "6. ",
  basename(cosohuc_summary_file),
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n"
)
