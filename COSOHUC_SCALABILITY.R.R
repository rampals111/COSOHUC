# ======================================================================
# COSOHUC SCALABILITY EXPERIMENT
# FINAL CLEAN VERSION - ALL 7 INITIALIZATION METHODS
# ======================================================================
#
# Initialization methods:
#   1. Farthest First
#   2. Canopy
#   3. K-Means++
#   4. Forgy
#   5. ECKM
#   6. Mini-Batch K-Means++
#   7. COSOHUC
#
# IMPORTANT:
#   - Same dataset is used for all methods within each condition/run.
#   - Same k is used for all methods.
#   - SAME Lloyd refinement is used after every initialization.
#   - Hartigan-Wong, MacQueen and standalone Lloyd are not treated as
#     initialization methods here.
#   - COSOHUC remains deterministic.
#   - ECKM remains deterministic for fixed X and k.
#   - Farthest First, K-Means++, Forgy and Mini-Batch K-Means++ use
#     run-specific fixed seeds.
#   - Dunn/DBI/Silhouette/RI are NOT calculated during profiling.
#   - No parallel implementation is used or claimed.
#
# NOTE ON CANOPY:
#   Your uploaded Canopy implementation constructs a full n x n
#   distance matrix. Therefore the common all-7 experiment uses
#   moderate n values. An optional extended experiment reports larger
#   n for the other six methods without silently changing Canopy.
# ======================================================================


# ======================================================================
# 1. PACKAGE
# ======================================================================

if (!requireNamespace("deldir", quietly = TRUE)) {
  install.packages("deldir", dependencies = TRUE)
}

if (!requireNamespace("peakRAM", quietly = TRUE)) {
  install.packages("peakRAM", dependencies = TRUE)
}


# ======================================================================
# 2. SETTINGS
# ======================================================================

OUTPUT_DIR <- "D:/kmean/DATA SET/data/usefull dataset/dataset my paper/scalability"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

N_RUNS <- 10

TIMING_REPEATS <- 3

K <- 3

MAX_ITER <- 100

BASE_SEED <- 20260902

# Common experiment: all 7 methods
OBSERVATION_SIZES_COMMON <- c(
  1000,
  2000,
  3000,
  5000
)

OBS_P <- 10

# Feature scaling: all 7 methods
FEATURE_N <- 5000

FEATURE_SIZES <- c(
  10,
  25,
  50,
  100
)

# Extended n: Canopy excluded because of exact n x n distance matrix
RUN_EXTENDED <- TRUE

OBSERVATION_SIZES_EXTENDED <- c(
  10000,
  20000
)

# Mini-Batch K-Means++ batch size = fraction of n
MINIBATCH_FRACTION <- 0.05

METHODS <- c(
  "Farthest First",
  "Canopy",
  "K-Means++",
  "Forgy",
  "ECKM",
  "Mini-Batch K-Means++",
  "COSOHUC"
)


# ======================================================================
# 3. WRITE-PERMISSION CHECK
# ======================================================================

test_file <- file.path(
  OUTPUT_DIR,
  "write_test.tmp"
)

write_ok <- tryCatch(
  {
    writeLines("test", test_file)
    TRUE
  },
  error = function(e) {
    FALSE
  }
)

if (file.exists(test_file)) {
  unlink(test_file)
}

if (!write_ok) {
  stop(
    paste0(
      "Cannot write to output directory:\n",
      OUTPUT_DIR,
      "\nPlease check folder permissions."
    )
  )
}


# ======================================================================
# 4. SYNTHETIC DATA GENERATOR
# ======================================================================
#
# Same X is supplied to every method for a given condition/run.
# All numerical features contribute to cluster separation.
# ======================================================================

generate_scalability_data <- function(
    n,
    p,
    k,
    seed
) {

  set.seed(seed)

  cluster_sizes <- rep(
    floor(n / k),
    k
  )

  remainder <- n - sum(cluster_sizes)

  if (remainder > 0) {
    cluster_sizes[seq_len(remainder)] <-
      cluster_sizes[seq_len(remainder)] + 1
  }

  center_spacing <- 6

  X_list <- vector(
    "list",
    k
  )

  for (cluster_id in seq_len(k)) {

    current_n <- cluster_sizes[cluster_id]

    X_cluster <- matrix(
      rnorm(
        current_n * p,
        mean = 0,
        sd = 1
      ),
      nrow = current_n,
      ncol = p
    )

    cluster_mean <-
      (cluster_id - 1) * center_spacing

    X_cluster <- sweep(
      X_cluster,
      2,
      cluster_mean,
      "+"
    )

    X_list[cluster_id] <- list(X_cluster)
  }

  X <- do.call(
    rbind,
    X_list
  )

  shuffle_index <- sample.int(
    nrow(X)
  )

  X <- X[
    shuffle_index,
    ,
    drop = FALSE
  ]

  colnames(X) <- paste0(
    "Feature_",
    seq_len(p)
  )

  X
}


validate_data <- function(
    X,
    k
) {

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("X must be numeric.")
  }

  if (any(!is.finite(X))) {
    stop("X contains NA, NaN or Inf.")
  }

  if (nrow(X) < k) {
    stop("nrow(X) must be >= k.")
  }

  if (ncol(X) < 1) {
    stop("X must contain at least one feature.")
  }

  X
}


# ======================================================================
# 5. COSOHUC INITIALIZATION
# ======================================================================
#
# Feature-wise:
#   sort -> partition into k equal-sized ordered groups -> mean
#
# The corresponding partition means across features form the k
# synthetic multivariate centers.
#
# ======================================================================

COSOHUC_initial_centers <- function(
    X,
    k
) {

  X <- validate_data(
    X,
    k
  )

  n <- nrow(X)
  p <- ncol(X)

  centers <- matrix(
    NA_real_,
    nrow = k,
    ncol = p
  )

  boundaries <- floor(
    seq(
      0,
      n,
      length.out = k + 1
    )
  )

  boundaries[1] <- 0
  boundaries[length(boundaries)] <- n

  for (feature_id in seq_len(p)) {

    sorted_values <- sort(
      X[, feature_id]
    )

    for (cluster_id in seq_len(k)) {

      start_index <-
        boundaries[cluster_id] + 1

      end_index <-
        boundaries[cluster_id + 1]

      centers[
        cluster_id,
        feature_id
      ] <- mean(
        sorted_values[
          start_index:end_index
        ]
      )
    }
  }

  colnames(centers) <- colnames(X)

  centers
}


# ======================================================================
# 6. FARTHEST FIRST
# ======================================================================

farthest_first_initial_centers <- function(
    X,
    k,
    seed
) {

  X <- validate_data(
    X,
    k
  )

  set.seed(seed)

  unique_X <- unique(X)

  if (nrow(unique_X) < k) {
    stop(
      paste0(
        "Farthest First requires ",
        k,
        " distinct observations."
      )
    )
  }

  n <- nrow(unique_X)
  p <- ncol(unique_X)

  selected <- integer(0)

  first_index <- sample.int(
    n,
    1
  )

  selected <- c(
    selected,
    first_index
  )

  min_distances <- rep(
    Inf,
    n
  )

  first_center <- unique_X[
    first_index,
    ,
    drop = TRUE
  ]

  for (i in seq_len(n)) {
    difference <-
      unique_X[
        i,
        ,
        drop = TRUE
      ] -
      first_center

    min_distances[i] <-
      sum(difference^2)
  }

  while (
    length(selected) < k
  ) {

    min_distances[selected] <- -Inf

    next_index <- which.max(
      min_distances
    )

    selected <- c(
      selected,
      next_index
    )

    new_center <- unique_X[
      next_index,
      ,
      drop = TRUE
    ]

    for (i in seq_len(n)) {

      if (!(i %in% selected)) {

        difference <-
          unique_X[
            i,
            ,
            drop = TRUE
          ] -
          new_center

        new_distance <-
          sum(
            difference^2
          )

        if (
          new_distance <
          min_distances[i]
        ) {
          min_distances[i] <-
            new_distance
        }
      }
    }
  }

  centers <- unique_X[
    selected,
    ,
    drop = FALSE
  ]

  colnames(centers) <- colnames(X)

  centers
}


# ======================================================================
# 7. K-MEANS++
# ======================================================================

kmeans_plus_plus_initial_centers <- function(
    X,
    k,
    seed
) {

  X <- validate_data(
    X,
    k
  )

  set.seed(seed)

  unique_X <- unique(X)

  if (nrow(unique_X) < k) {
    stop(
      paste0(
        "K-Means++ requires ",
        k,
        " distinct observations."
      )
    )
  }

  n <- nrow(unique_X)

  selected <- integer(0)

  first_index <- sample.int(
    n,
    1
  )

  selected <- c(
    selected,
    first_index
  )

  while (
    length(selected) < k
  ) {

    candidates <- setdiff(
      seq_len(n),
      selected
    )

    distances_squared <- numeric(
      length(candidates)
    )

    for (a in seq_along(candidates)) {

      i <- candidates[a]

      point <- unique_X[
        i,
        ,
        drop = TRUE
      ]

      distances_to_selected <- vapply(
        selected,
        function(s) {

          difference <-
            point -
            unique_X[
              s,
              ,
              drop = TRUE
            ]

          sum(
            difference^2
          )
        },
        numeric(1)
      )

      distances_squared[a] <-
        min(
          distances_to_selected
        )
    }

    total <- sum(
      distances_squared
    )

    if (
      !is.finite(total) ||
      total <= 0
    ) {

      next_point <- sample(
        candidates,
        1
      )

    } else {

      probabilities <-
        distances_squared /
        total

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


# ======================================================================
# 8. FORGY
# ======================================================================

forgy_initial_centers <- function(
    X,
    k,
    seed
) {

  X <- validate_data(
    X,
    k
  )

  set.seed(seed)

  unique_X <- unique(X)

  if (nrow(unique_X) < k) {
    stop(
      paste0(
        "Forgy requires ",
        k,
        " distinct observations."
      )
    )
  }

  selected <- sample(
    seq_len(
      nrow(unique_X)
    ),
    size = k,
    replace = FALSE
  )

  centers <- unique_X[
    selected,
    ,
    drop = FALSE
  ]

  colnames(centers) <- colnames(X)

  centers
}


# ======================================================================
# 9. CANOPY
# ======================================================================
#
# This retains the exact computational characteristic of your uploaded
# Canopy logic: full pairwise distance matrix.
#
# T2 = 25th percentile of positive distances
# T1 = 50th percentile of positive distances
#
# The selection itself uses T2, followed by Farthest First fallback.
#
# ======================================================================

canopy_initial <- function(
    X,
    k,
    seed
) {

  X <- validate_data(
    X,
    k
  )

  set.seed(seed)

  n <- nrow(X)

  D <- as.matrix(
    dist(X)
  )

  positive_distances <- D[
    D > 0
  ]

  positive_distances <- positive_distances[
    is.finite(
      positive_distances
    )
  ]

  if (
    length(positive_distances) == 0
  ) {
    stop(
      "Canopy cannot be constructed because all distances are zero."
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

  if (T1 <= T2) {
    T1 <- T2 * 1.5
  }

  remaining <- seq_len(n)
  selected <- integer(0)

  while (
    length(remaining) > 0 &&
    length(selected) < k
  ) {

    center_id <- remaining[1]

    selected <- c(
      selected,
      center_id
    )

    distances <- D[
      center_id,
      remaining
    ]

    remove_ids <- remaining[
      distances <= T2
    ]

    remaining <- setdiff(
      remaining,
      remove_ids
    )
  }

  # Farthest First fallback, matching the uploaded structure.
  if (
    length(selected) < k
  ) {

    ff_centers <- farthest_first_initial_centers(
      X,
      k,
      seed
    )

    for (
      i in seq_len(
        nrow(ff_centers)
      )
    ) {

      if (
        length(selected) >= k
      ) {
        break
      }

      point <- ff_centers[
        i,
        ,
        drop = TRUE
      ]

      distances <- rowSums(
        (
          X -
          matrix(
            point,
            nrow = n,
            ncol = ncol(X),
            byrow = TRUE
          )
        )^2
      )

      candidate <- which.min(
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
      "Canopy could not produce k centers."
    )
  }

  selected <- selected[
    seq_len(k)
  ]

  centers <- X[
    selected,
    ,
    drop = FALSE
  ]

  if (
    nrow(
      unique(centers)
    ) < k
  ) {

    centers <-
      farthest_first_initial_centers(
        X,
        k,
        seed
      )
  }

  colnames(centers) <- colnames(X)

  centers
}


# ======================================================================
# 10. ECKM
# ======================================================================
#
# Structure used in your uploaded implementation:
#
#   PCA for p > 2
#      -> 2-D geometry
#      -> unique points
#      -> Voronoi vertices
#      -> interior vertices
#      -> empty-circle radii
#      -> descending radius order
#      -> circumference points
#      -> map to distinct original observations
#
# ======================================================================

point_inside_convex_hull <- function(
    point,
    hull_points
) {

  hull_index <- chull(
    hull_points[, 1],
    hull_points[, 2]
  )

  hull <- hull_points[
    c(
      hull_index,
      hull_index[1]
    ),
    ,
    drop = FALSE
  ]

  n <- nrow(hull) - 1

  if (n < 3) {
    return(FALSE)
  }

  x <- point[1]
  y <- point[2]

  inside <- FALSE

  j <- n

  for (
    i in seq_len(n)
  ) {

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


calculate_voronoi_vertices_eckm <- function(
    points
) {

  if (nrow(points) < 3) {
    stop(
      "ECKM requires at least 3 unique points."
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

  vertices_list <- vector(
    "list",
    length(tiles)
  )

  counter <- 0L

  for (
    tile in tiles
  ) {

    tile_x <- tile$x
    tile_y <- tile$y

    if (
      length(tile_x) < 3
    ) {
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

    vertices_list[counter] <- list(c(ux, uy))
  }

  if (counter == 0) {
    stop(
      "No Voronoi vertices were obtained."
    )
  }

  vertices <- do.call(
    rbind,
    vertices_list[seq_len(counter)]
  )

  unique(
    round(
      vertices,
      12
    )
  )
}


prepare_eckm_geometry <- function(
    X,
    use_pca = TRUE
) {

  if (ncol(X) < 2) {
    stop(
      "ECKM requires at least two features."
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

  geometry_data_unique <- geometry_data[
    !duplicated(
      geometry_data
    ),
    ,
    drop = FALSE
  ]

  if (
    nrow(geometry_data_unique) < 3
  ) {
    stop(
      "Insufficient unique geometry points."
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

  interior_vertices <- matrix(
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
      ,
      drop = TRUE
    ]

    if (
      point_inside_convex_hull(
        point,
        hull_points
      )
    ) {

      interior_vertices <- rbind(
        interior_vertices,
        point
      )
    }
  }

  if (
    nrow(interior_vertices) == 0
  ) {
    stop(
      "No interior Voronoi vertices were obtained."
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

    differences <- sweep(
      points,
      2,
      centers[
        i,
        ],
      "-"
    )

    distances <- sqrt(
      rowSums(
        differences^2
      )
    )

    radius_values[i] <- min(
      distances
    )
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

  number_requested <- ceiling(
    k / 3
  )

  repeat {

    if (
      number_requested >
      length(circle_radii)
    ) {
      stop(
        "Not enough empty circles for ECKM."
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

    all_indices <- integer(0)

    for (
      i in seq_len(
        nrow(current_centers)
      )
    ) {

      difference <- sweep(
        data_points,
        2,
        current_centers[
          i,
          ],
        "-"
      )

      distances <- sqrt(
        rowSums(
          difference^2
        )
      )

      difference_from_circle <- abs(
        distances -
        current_radii[i]
      )

      point_indices <- which(
        round(
          difference_from_circle,
          3
        ) < epsilon
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

    all_indices <- unique(
      all_indices
    )

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

  selected_original_indices <- integer(0)

  for (
    i in seq_len(
      nrow(eckm_geometry_centers)
    )
  ) {

    differences <- sweep(
      geometry_data,
      2,
      eckm_geometry_centers[
        i,
        ],
      "-"
    )

    distances <- sqrt(
      rowSums(
        differences^2
      )
    )

    candidate_order <- order(
      distances
    )

    candidate <- candidate_order[1]

    if (
      candidate %in%
      selected_original_indices
    ) {

      new_candidate <- candidate_order[
        !(candidate_order %in%
            selected_original_indices)
      ]

      if (
        length(new_candidate) == 0
      ) {
        stop(
          "ECKM could not map k distinct observations."
        )
      }

      candidate <- new_candidate[1]
    }

    selected_original_indices <- c(
      selected_original_indices,
      candidate
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

  if (nrow(X) < k) {
    stop(
      "ECKM requires n >= k."
    )
  }

  geometry <- prepare_eckm_geometry(
    X,
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

  lec_order <- order(
    radius_values,
    decreasing = TRUE
  )

  lec_centers <- interior_vertices[
    lec_order,
    ,
    drop = FALSE
  ]

  lec_radii <- radius_values[
    lec_order
  ]

  circumference <- find_circumference_points_eckm(
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

  if (
    length(
      circumference$point_indices
    ) < k
  ) {
    stop(
      "ECKM could not generate k circumference points."
    )
  }

  selected_geometry_indices <-
    circumference$point_indices[
      seq_len(k)
    ]

  eckm_geometry_centers <-
    geometry_data_unique[
      selected_geometry_indices,
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

  initial_centers <- X[
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
    selected_lec_count =
      circumference$selected_count
  )
}


# ======================================================================
# 11. MINI-BATCH K-MEANS++
# ======================================================================

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

    start_id <-
      (batch_id - 1) * b + 1

    end_id <-
      min(
        batch_id * b,
        n
      )

    batch_data <- X[
      start_id:end_id,
      ,
      drop = FALSE
    ]

    running_coordinate_sum <-
      running_coordinate_sum +
      colSums(batch_data)

    running_squared_sum <-
      running_squared_sum +
      sum(batch_data^2)

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
    ] <- start_id

    batch_end[
      batch_id
    ] <- end_id
  }

  list(
    B =
      B,
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

  center_squared_sum <-
    sum(
      selected_centers^2
    )

  centers_coordinate_sum <-
    colSums(
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

  while (
    left < right
  ) {

    middle <- floor(
      (
        left +
        right
      ) / 2
    )

    cumulative_n <- batch_terms$batch_end[
      middle
    ]

    cumulative_coordinate_sum <-
      batch_terms$
      cumulative_coordinate_sums[
        middle,
        ,
        drop = TRUE
      ]

    cumulative_squared_sum <-
      batch_terms$
      cumulative_squared_sums[
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
    global_id in batch_ids
  ) {

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

  batch_ids[
    length(batch_ids)
  ]
}


mini_batch_kmeans_plus_plus_initial_centers <- function(
    X,
    k,
    b,
    seed
) {

  X <- validate_data(
    X,
    k
  )

  set.seed(seed)

  n <- nrow(X)

  if (
    k > n
  ) {
    stop(
      "k cannot exceed n."
    )
  }

  batch_terms <- precompute_batch_terms(
    X,
    b
  )

  first_index <- sample(
    seq_len(n),
    size = 1
  )

  selected_indices <-
    first_index

  selected_centers <-
    X[
      first_index,
      ,
      drop = FALSE
    ]

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

    } else {

      R <- runif(
        1,
        min = 0,
        max = total_distance
      )

      target_batch <-
        find_target_batch(
          R =
            R,
          selected_centers =
            selected_centers,
          batch_terms =
            batch_terms
        )

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
          batch_terms$
          cumulative_coordinate_sums[
            target_batch - 1,
            ,
            drop = TRUE
          ]

        cumulative_squared_sum_before <-
          batch_terms$
          cumulative_squared_sums[
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

      batch_start <-
        batch_terms$batch_start[
          target_batch
        ]

      batch_end <-
        batch_terms$batch_end[
          target_batch
        ]

      next_index <-
        find_point_inside_batch(

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
    }

    if (
      next_index %in%
      selected_indices
    ) {

      remaining <- setdiff(
        seq_len(n),
        selected_indices
      )

      if (
        length(remaining) == 0
      ) {
        stop(
          "No remaining observations for Mini-Batch K-Means++."
        )
      }

      next_index <- sample(
        remaining,
        size = 1
      )
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

  colnames(selected_centers) <-
    colnames(X)

  selected_centers
}


# ======================================================================
# 12. DISPATCHER
# ======================================================================

get_initial_centers <- function(
    X,
    method,
    k,
    seed,
    minibatch_size
) {

  if (
    method ==
    "COSOHUC"
  ) {
    return(
      COSOHUC_initial_centers(
        X,
        k
      )
    )
  }

  if (
    method ==
    "Farthest First"
  ) {
    return(
      farthest_first_initial_centers(
        X,
        k,
        seed
      )
    )
  }

  if (
    method ==
    "Canopy"
  ) {
    return(
      canopy_initial(
        X,
        k,
        seed
      )
    )
  }

  if (
    method ==
    "K-Means++"
  ) {
    return(
      kmeans_plus_plus_initial_centers(
        X,
        k,
        seed
      )
    )
  }

  if (
    method ==
    "Forgy"
  ) {
    return(
      forgy_initial_centers(
        X,
        k,
        seed
      )
    )
  }

  if (
    method ==
    "ECKM"
  ) {

    eckm <- ECKM_initial_centers(
      X,
      k
    )

    return(
      eckm$centers
    )
  }

  if (
    method ==
    "Mini-Batch K-Means++"
  ) {

    return(
      mini_batch_kmeans_plus_plus_initial_centers(

        X =
          X,

        k =
          k,

        b =
          minibatch_size,

        seed =
          seed

      )
    )
  }

  stop(
    paste(
      "Unknown method:",
      method
    )
  )
}


# ======================================================================
# 13. COMMON LLOYD REFINEMENT
# ======================================================================
#
# Every initialization method uses this exact refinement.
# ======================================================================

run_common_lloyd <- function(
    X,
    centers,
    max_iter = 100
) {

  X <- validate_data(
    X,
    nrow(centers)
  )

  if (
    any(
      !is.finite(
        centers
      )
    )
  ) {
    stop(
      "Initial centers contain non-finite values."
    )
  }

  if (
    nrow(
      unique(
        centers
      )
    ) <
    nrow(centers)
  ) {
    stop(
      "Initial centers are duplicated."
    )
  }

  stats::kmeans(
    x =
      X,
    centers =
      centers,
    nstart =
      1,
    iter.max =
      max_iter,
    algorithm =
      "Lloyd"
  )
}


# ======================================================================
# 14. SAFE INITIALIZATION
# ======================================================================

safe_get_initial_centers <- function(
    X,
    method,
    k,
    seed,
    minibatch_size
) {

  tryCatch(

    {

      get_initial_centers(

        X =
          X,

        method =
          method,

        k =
          k,

        seed =
          seed,

        minibatch_size =
          minibatch_size

      )

    },

    error = function(e) {

      stop(
        paste0(
          "Initialization failed: ",
          conditionMessage(e)
        )
      )
    }
  )
}


# ======================================================================
# 15. TIME INITIALIZATION
# ======================================================================
#
# Timing repeats are only for measurement.
# COSOHUC remains deterministic.
# ======================================================================

measure_initialization <- function(
    X,
    method,
    k,
    seed,
    minibatch_size,
    repeats
) {

  elapsed <- numeric(
    repeats
  )

  for (
    r in seq_len(repeats)
  ) {

    current_seed <-
      seed + r

    timing <-
      system.time(

        {

          centers <-
            safe_get_initial_centers(

              X =
                X,

              method =
                method,

              k =
                k,

              seed =
                current_seed,

              minibatch_size =
                minibatch_size

            )

          invisible(
            centers
          )
        }

      )

    elapsed[r] <-
      unname(
        timing["elapsed"]
      ) *
      1000
  }

  finite_elapsed <- elapsed[
    is.finite(elapsed)
  ]

  if (
    length(finite_elapsed) == 0
  ) {
    stop(
      "No valid initialization timing was obtained."
    )
  }

  list(

    mean_ms =
      mean(
        finite_elapsed
      ),

    median_ms =
      median(
        finite_elapsed
      ),

    sd_ms =
      if (
        length(finite_elapsed) > 1
      ) {
        sd(finite_elapsed)
      } else {
        NA_real_
      }

  )
}


# ======================================================================
# 16. TIME TOTAL PIPELINE
# ======================================================================

measure_total_pipeline <- function(
    X,
    method,
    k,
    seed,
    minibatch_size,
    max_iter,
    repeats
) {

  elapsed <- numeric(
    repeats
  )

  for (
    r in seq_len(repeats)
  ) {

    current_seed <-
      seed +
      10000 +
      r

    timing <-
      system.time(

        {

          centers <-
            safe_get_initial_centers(

              X =
                X,

              method =
                method,

              k =
                k,

              seed =
                current_seed,

              minibatch_size =
                minibatch_size

            )

          result <-
            run_common_lloyd(

              X =
                X,

              centers =
                centers,

              max_iter =
                max_iter

            )

          invisible(
            result
          )
        }

      )

    elapsed[r] <-
      unname(
        timing["elapsed"]
      ) *
      1000
  }

  finite_elapsed <- elapsed[
    is.finite(elapsed)
  ]

  if (
    length(finite_elapsed) == 0
  ) {
    stop(
      "No valid total timing was obtained."
    )
  }

  list(

    mean_ms =
      mean(
        finite_elapsed
      ),

    median_ms =
      median(
        finite_elapsed
      ),

    sd_ms =
      if (
        length(finite_elapsed) > 1
      ) {
        sd(finite_elapsed)
      } else {
        NA_real_
      }

  )
}


# ======================================================================
# 17. PEAK MEMORY
# ======================================================================

measure_peak_memory <- function(
    X,
    method,
    k,
    seed,
    minibatch_size,
    max_iter
) {

  holder <-
    new.env(
      parent = emptyenv()
    )

  memory_result <-
    tryCatch(

      peakRAM::peakRAM(

        {

          centers <-
            safe_get_initial_centers(

              X =
                X,

              method =
                method,

              k =
                k,

              seed =
                seed,

              minibatch_size =
                minibatch_size

            )

          holder$result <-
            run_common_lloyd(

              X =
                X,

              centers =
                centers,

              max_iter =
                max_iter

            )

        }

      ),

      error = function(e) {

        NULL
      }
    )

  if (
    is.null(
      memory_result
    )
  ) {
    return(
      NA_real_
    )
  }

  if (
    "Peak_RAM_Used_MiB" %in%
    names(memory_result)
  ) {

    return(
      as.numeric(
        memory_result[
          "Peak_RAM_Used_MiB"
        ][1]
      )
    )
  }

  NA_real_
}


# ======================================================================
# 18. ONE METHOD/CONDITION
# ======================================================================

benchmark_method <- function(
    X,
    method,
    k,
    seed,
    minibatch_size,
    max_iter,
    timing_repeats
) {

  init_time <-
    measure_initialization(

      X =
        X,

      method =
        method,

      k =
        k,

      seed =
        seed,

      minibatch_size =
        minibatch_size,

      repeats =
        timing_repeats

    )


  total_time <-
    measure_total_pipeline(

      X =
        X,

      method =
        method,

      k =
        k,

      seed =
        seed,

      minibatch_size =
        minibatch_size,

      max_iter =
        max_iter,

      repeats =
        timing_repeats

    )


  peak_memory <-
    measure_peak_memory(

      X =
        X,

      method =
        method,

      k =
        k,

      seed =
        seed,

      minibatch_size =
        minibatch_size,

      max_iter =
        max_iter

    )


  # ------------------------------------------------------------
  # One actual final run
  # ------------------------------------------------------------

  final_centers <-
    safe_get_initial_centers(

      X =
        X,

      method =
        method,

      k =
        k,

      seed =
        seed,

      minibatch_size =
        minibatch_size

    )


  final_result <-
    run_common_lloyd(

      X =
        X,

      centers =
        final_centers,

      max_iter =
        max_iter

    )


  list(

    Status =
      "OK",

    Initialization_Mean_ms =
      init_time$mean_ms,

    Initialization_Median_ms =
      init_time$median_ms,

    Initialization_SD_ms =
      init_time$sd_ms,

    Total_Mean_ms =
      total_time$mean_ms,

    Total_Median_ms =
      total_time$median_ms,

    Total_SD_ms =
      total_time$sd_ms,

    Peak_Memory_MiB =
      peak_memory,

    Iterations =
      as.numeric(
        final_result$iter
      ),

    Error =
      ""

  )
}


# ======================================================================
# 19. SAFE WRAPPER
# ======================================================================

benchmark_method_safe <- function(
    X,
    method,
    k,
    seed,
    minibatch_size,
    max_iter,
    timing_repeats
) {

  tryCatch(

    {

      benchmark_method(

        X =
          X,

        method =
          method,

        k =
          k,

        seed =
          seed,

        minibatch_size =
          minibatch_size,

        max_iter =
          max_iter,

        timing_repeats =
          timing_repeats

      )

    },

    error = function(e) {

      list(

        Status =
          "FAILED",

        Initialization_Mean_ms =
          NA_real_,

        Initialization_Median_ms =
          NA_real_,

        Initialization_SD_ms =
          NA_real_,

        Total_Mean_ms =
          NA_real_,

        Total_Median_ms =
          NA_real_,

        Total_SD_ms =
          NA_real_,

        Peak_Memory_MiB =
          NA_real_,

        Iterations =
          NA_real_,

        Error =
          conditionMessage(e)

      )
    }
  )
}


# ======================================================================
# 20. RUN ONE DATASET CONDITION
# ======================================================================

run_condition <- function(
    X,
    experiment_name,
    dataset_label,
    methods_to_run
) {

  results <- data.frame()

  n_current <-
    nrow(X)

  p_current <-
    ncol(X)

  minibatch_size <-
    max(
      10,
      min(
        n_current,
        ceiling(
          n_current *
          MINIBATCH_FRACTION
        )
      )
    )

  cat(
    "\n============================================================\n"
  )

  cat(
    experiment_name,
    "\n"
  )

  cat(
    "Dataset:",
    dataset_label,
    "\n"
  )

  cat(
    "n =",
    n_current,
    " p =",
    p_current,
    " k =",
    K,
    "\n"
  )

  cat(
    "Mini-batch size =",
    minibatch_size,
    "\n"
  )

  cat(
    "============================================================\n"
  )


  for (
    method_id in seq_along(
      methods_to_run
    )
  ) {

    method <-
      methods_to_run[
        method_id
      ]

    cat(
      "\nMethod:",
      method,
      "\n"
    )


    for (
      run_id in seq_len(
        N_RUNS
      )
    ) {

      cat(
        "  Run",
        run_id,
        "/",
        N_RUNS,
        "\n"
      )


      seed <-
        BASE_SEED +
        n_current * 10 +
        p_current * 100 +
        run_id * 1000 +
        method_id


      result <-
        benchmark_method_safe(

          X =
            X,

          method =
            method,

          k =
            K,

          seed =
            seed,

          minibatch_size =
            minibatch_size,

          max_iter =
            MAX_ITER,

          timing_repeats =
            TIMING_REPEATS

        )


      results <-
        rbind(

          results,

          data.frame(

            Experiment =
              experiment_name,

            Dataset =
              dataset_label,

            Observations =
              n_current,

            Features =
              p_current,

            Clusters =
              K,

            MiniBatch_Size =
              minibatch_size,

            Run =
              run_id,

            Method =
              method,

            Status =
              result$Status,

            Initialization_Mean_ms =
              result$Initialization_Mean_ms,

            Initialization_Median_ms =
              result$Initialization_Median_ms,

            Initialization_SD_ms =
              result$Initialization_SD_ms,

            Total_Mean_ms =
              result$Total_Mean_ms,

            Total_Median_ms =
              result$Total_Median_ms,

            Total_SD_ms =
              result$Total_SD_ms,

            Peak_Memory_MiB =
              result$Peak_Memory_MiB,

            Iterations =
              result$Iterations,

            Error =
              result$Error,

            stringsAsFactors =
              FALSE

          )

        )

    }
  }

  results
}


# ======================================================================
# 21. OBSERVATION SCALING
# ======================================================================

observation_results <- data.frame()

for (
  n_current in OBSERVATION_SIZES_COMMON
) {

  X <- generate_scalability_data(

    n =
      n_current,

    p =
      OBS_P,

    k =
      K,

    seed =
      BASE_SEED +
      n_current

  )


  current_results <-
    run_condition(

      X =
        X,

      experiment_name =
        "Observation Scaling - All 7 Methods",

      dataset_label =
        paste0(
          n_current,
          "x",
          OBS_P
        ),

      methods_to_run =
        METHODS

    )


  observation_results <-
    rbind(

      observation_results,

      current_results

    )


  rm(X)

  gc()
}


# ======================================================================
# 22. SAVE OBSERVATION RAW
# ======================================================================

observation_raw_file <-
  file.path(

    OUTPUT_DIR,

    "SCALABILITY_ALL7_OBSERVATION_RAW.csv"

  )


write.csv(

  observation_results,

  observation_raw_file,

  row.names = FALSE

)


# ======================================================================
# 23. OBSERVATION SUMMARY
# ======================================================================

valid_observation <-
  observation_results[
    observation_results$Status ==
      "OK",
    ,
    drop = FALSE
  ]


observation_summary <-
  if (
    nrow(valid_observation) > 0
  ) {

    aggregate(

      cbind(

        Initialization_Mean_ms,

        Initialization_Median_ms,

        Total_Mean_ms,

        Total_Median_ms,

        Peak_Memory_MiB,

        Iterations

      )

      ~

        Observations +
        Features +
        Clusters +
        Method,

      data =
        valid_observation,

      FUN =
        function(x) {

          mean(
            x,
            na.rm = TRUE
          )

        }

    )

  } else {

    data.frame()
  }


# ======================================================================
# 24. FEATURE SCALING
# ======================================================================

feature_results <- data.frame()


for (
  p_current in FEATURE_SIZES
) {

  X <- generate_scalability_data(

    n =
      FEATURE_N,

    p =
      p_current,

    k =
      K,

    seed =
      BASE_SEED +
      500000 +
      p_current

  )


  current_results <-
    run_condition(

      X =
        X,

      experiment_name =
        "Feature Scaling - All 7 Methods",

      dataset_label =
        paste0(
          FEATURE_N,
          "x",
          p_current
        ),

      methods_to_run =
        METHODS

    )


  feature_results <-
    rbind(

      feature_results,

      current_results

    )


  rm(X)

  gc()
}


# ======================================================================
# 25. SAVE FEATURE RAW
# ======================================================================

feature_raw_file <-
  file.path(

    OUTPUT_DIR,

    "SCALABILITY_ALL7_FEATURE_RAW.csv"

  )


write.csv(

  feature_results,

  feature_raw_file,

  row.names = FALSE

)


# ======================================================================
# 26. FEATURE SUMMARY
# ======================================================================

valid_feature <-
  feature_results[
    feature_results$Status ==
      "OK",
    ,
    drop = FALSE
  ]


feature_summary <-
  if (
    nrow(valid_feature) > 0
  ) {

    aggregate(

      cbind(

        Initialization_Mean_ms,

        Initialization_Median_ms,

        Total_Mean_ms,

        Total_Median_ms,

        Peak_Memory_MiB,

        Iterations

      )

      ~

        Observations +
        Features +
        Clusters +
        Method,

      data =
        valid_feature,

      FUN =
        function(x) {

          mean(
            x,
            na.rm = TRUE
          )

        }

    )

  } else {

    data.frame()
  }


# ======================================================================
# 27. EXTENDED OBSERVATION SCALING
# ======================================================================
#
# Canopy is excluded ONLY here because the exact implementation uses
# a full n x n distance matrix.
#
# This is clearly labeled; no approximate Canopy replacement is used.
#
# ======================================================================

extended_results <- data.frame()

extended_summary <- data.frame()

if (
  RUN_EXTENDED
) {

  EXTENDED_METHODS <- c(

    "Farthest First",

    "K-Means++",

    "Forgy",

    "ECKM",

    "Mini-Batch K-Means++",

    "COSOHUC"

  )


  for (
    n_current in
    OBSERVATION_SIZES_EXTENDED
  ) {

    X <- generate_scalability_data(

      n =
        n_current,

      p =
        OBS_P,

      k =
        K,

      seed =
        BASE_SEED +
        900000 +
        n_current

    )


    current_results <-
      run_condition(

        X =
          X,

        experiment_name =
          "Extended Observation Scaling - 6 Methods",

        dataset_label =
          paste0(
            n_current,
            "x",
            OBS_P
          ),

        methods_to_run =
          EXTENDED_METHODS

      )


    extended_results <-
      rbind(

        extended_results,

        current_results

      )


    rm(X)

    gc()
  }


  extended_raw_file <-
    file.path(

      OUTPUT_DIR,

      "SCALABILITY_EXTENDED_OBSERVATION_RAW.csv"

    )


  write.csv(

    extended_results,

    extended_raw_file,

    row.names = FALSE

  )


  valid_extended <-
    extended_results[
      extended_results$Status ==
        "OK",
      ,
      drop = FALSE
    ]


  if (
    nrow(valid_extended) > 0
  ) {

    extended_summary <-
      aggregate(

        cbind(

          Initialization_Mean_ms,

          Initialization_Median_ms,

          Total_Mean_ms,

          Total_Median_ms,

          Peak_Memory_MiB,

          Iterations

        )

        ~

          Observations +
          Features +
          Clusters +
          Method,

        data =
          valid_extended,

        FUN =
          function(x) {

            mean(
              x,
              na.rm = TRUE
            )

          }

      )

  }
}


# ======================================================================
# 28. FIX METHOD ORDER
# ======================================================================

method_order <- c(

  "Farthest First",

  "Canopy",

  "K-Means++",

  "Forgy",

  "ECKM",

  "Mini-Batch K-Means++",

  "COSOHUC"

)


order_summary <- function(
    x
) {

  if (
    nrow(x) == 0
  ) {

    return(
      x
    )
  }

  x$Method <- factor(
    x$Method,
    levels =
      method_order
  )

  x <-
    x[
      order(
        x$Observations,
        x$Features,
        x$Method
      ),
      ,
      drop = FALSE
    ]

  x$Method <-
    as.character(
      x$Method
    )

  x
}


observation_summary <-
  order_summary(
    observation_summary
  )


feature_summary <-
  order_summary(
    feature_summary
  )


extended_summary <-
  order_summary(
    extended_summary
  )


# ======================================================================
# 29. SAVE SUMMARIES
# ======================================================================

observation_summary_file <-
  file.path(

    OUTPUT_DIR,

    "SCALABILITY_ALL7_OBSERVATION_SUMMARY.csv"

  )


feature_summary_file <-
  file.path(

    OUTPUT_DIR,

    "SCALABILITY_ALL7_FEATURE_SUMMARY.csv"

  )


extended_summary_file <-
  file.path(

    OUTPUT_DIR,

    "SCALABILITY_EXTENDED_OBSERVATION_SUMMARY.csv"

  )


write.csv(

  observation_summary,

  observation_summary_file,

  row.names = FALSE

)


write.csv(

  feature_summary,

  feature_summary_file,

  row.names = FALSE

)


if (
  RUN_EXTENDED
) {

  write.csv(

    extended_summary,

    extended_summary_file,

    row.names = FALSE

  )
}


# ======================================================================
# 30. PLOTTING
# ======================================================================

METHOD_LTY <- c(

  "Farthest First" = 1,

  "Canopy" = 2,

  "K-Means++" = 3,

  "Forgy" = 4,

  "ECKM" = 5,

  "Mini-Batch K-Means++" = 6,

  "COSOHUC" = 1

)


METHOD_PCH <- c(

  "Farthest First" = 16,

  "Canopy" = 17,

  "K-Means++" = 15,

  "Forgy" = 18,

  "ECKM" = 8,

  "Mini-Batch K-Means++" = 7,

  "COSOHUC" = 0

)


make_plot <- function(

    summary_data,

    x_column,

    y_column,

    x_values,

    x_label,

    y_label,

    title,

    output_file,

    log_x = FALSE

) {

  if (
    nrow(summary_data) == 0
  ) {
    return(
      invisible(NULL)
    )
  }


  max_y <-
    max(
      summary_data[
        ,
        y_column
      ],
      na.rm = TRUE
    )


  if (
    !is.finite(max_y) ||
    max_y <= 0
  ) {
    max_y <- 1
  }


  png(

    output_file,

    width =
      1600,

    height =
      1000

  )


  plot(

    x_values,

    rep(
      NA_real_,
      length(x_values)
    ),

    type =
      "n",

    log =
      if (log_x) "x" else "",

    xlim =
      range(x_values),

    ylim =
      c(
        0,
        max_y * 1.15
      ),

    xlab =
      x_label,

    ylab =
      y_label,

    main =
      title

  )


  methods_present <-
    unique(
      summary_data$Method
    )


  for (
    method in method_order
  ) {

    if (
      !(method %in%
        methods_present)
    ) {
      next
    }


    temp <-
      summary_data[
        summary_data$Method ==
          method,
        ,
        drop = FALSE
      ]


    temp <-
      temp[
        order(
          temp[
            ,
            x_column
          ]
        ),
        ,
        drop = FALSE
      ]


    valid <- is.finite(
      temp[
        ,
        y_column
      ]
    )


    temp <-
      temp[
        valid,
        ,
        drop = FALSE
      ]


    if (
      nrow(temp) == 0
    ) {
      next
    }


    lines(

      temp[
        ,
        x_column
      ],

      temp[
        ,
        y_column
      ],

      type =
        "b",

      lty =
        METHOD_LTY[
          method
        ],

      pch =
        METHOD_PCH[
          method
        ]

    )
  }


  legend(

    "topleft",

    legend =
      method_order[
        method_order %in%
          methods_present
      ],

    lty =
      METHOD_LTY[
        method_order[
          method_order %in%
            methods_present
        ]
      ],

    pch =
      METHOD_PCH[
        method_order[
          method_order %in%
            methods_present
        ]
      ],

    bty =
      "o",

    cex =
      0.85

  )


  dev.off()


  invisible(NULL)
}


# ======================================================================
# 31. OBSERVATION PLOTS
# ======================================================================

make_plot(

  summary_data =
    observation_summary,

  x_column =
    "Observations",

  y_column =
    "Initialization_Median_ms",

  x_values =
    OBSERVATION_SIZES_COMMON,

  x_label =
    "Number of Observations",

  y_label =
    "Median Initialization Time (ms)",

  title =
    "Initialization Time vs Number of Observations",

  output_file =
    file.path(

      OUTPUT_DIR,

      "FIG_01_Observation_Initialization.png"

    ),

  log_x =
    TRUE

)


make_plot(

  summary_data =
    observation_summary,

  x_column =
    "Observations",

  y_column =
    "Total_Median_ms",

  x_values =
    OBSERVATION_SIZES_COMMON,

  x_label =
    "Number of Observations",

  y_label =
    "Median Total Pipeline Time (ms)",

  title =
    "Total K-means Pipeline Time vs Number of Observations",

  output_file =
    file.path(

      OUTPUT_DIR,

      "FIG_02_Observation_Total_Time.png"

    ),

  log_x =
    TRUE

)


make_plot(

  summary_data =
    observation_summary,

  x_column =
    "Observations",

  y_column =
    "Peak_Memory_MiB",

  x_values =
    OBSERVATION_SIZES_COMMON,

  x_label =
    "Number of Observations",

  y_label =
    "Peak Memory (MiB)",

  title =
    "Peak Memory vs Number of Observations",

  output_file =
    file.path(

      OUTPUT_DIR,

      "FIG_03_Observation_Peak_Memory.png"

    ),

  log_x =
    TRUE

)


# ======================================================================
# 32. FEATURE PLOTS
# ======================================================================

make_plot(

  summary_data =
    feature_summary,

  x_column =
    "Features",

  y_column =
    "Initialization_Median_ms",

  x_values =
    FEATURE_SIZES,

  x_label =
    "Number of Features",

  y_label =
    "Median Initialization Time (ms)",

  title =
    "Initialization Time vs Number of Features",

  output_file =
    file.path(

      OUTPUT_DIR,

      "FIG_04_Feature_Initialization.png"

    ),

  log_x =
    FALSE

)


make_plot(

  summary_data =
    feature_summary,

  x_column =
    "Features",

  y_column =
    "Total_Median_ms",

  x_values =
    FEATURE_SIZES,

  x_label =
    "Number of Features",

  y_label =
    "Median Total Pipeline Time (ms)",

  title =
    "Total K-means Pipeline Time vs Number of Features",

  output_file =
    file.path(

      OUTPUT_DIR,

      "FIG_05_Feature_Total_Time.png"

    ),

  log_x =
    FALSE

)


make_plot(

  summary_data =
    feature_summary,

  x_column =
    "Features",

  y_column =
    "Peak_Memory_MiB",

  x_values =
    FEATURE_SIZES,

  x_label =
    "Number of Features",

  y_label =
    "Peak Memory (MiB)",

  title =
    "Peak Memory vs Number of Features",

  output_file =
    file.path(

      OUTPUT_DIR,

      "FIG_06_Feature_Peak_Memory.png"

    ),

  log_x =
    FALSE

)


# ======================================================================
# 33. EXPERIMENT SETTINGS FILE
# ======================================================================

settings_file <-
  file.path(

    OUTPUT_DIR,

    "SCALABILITY_EXPERIMENT_SETTINGS.txt"

  )


settings_text <- c(

  "COSOHUC SCALABILITY EXPERIMENT",

  "",

  paste(
    "Independent runs:",
    N_RUNS
  ),

  paste(
    "Timing repetitions:",
    TIMING_REPEATS
  ),

  paste(
    "k:",
    K
  ),

  paste(
    "Maximum Lloyd iterations:",
    MAX_ITER
  ),

  "",

  "Common observation scaling:",

  paste(
    OBSERVATION_SIZES_COMMON,
    collapse = ", "
  ),

  paste(
    "Fixed p:",
    OBS_P
  ),

  "",

  "Feature scaling:",

  paste(
    FEATURE_SIZES,
    collapse = ", "
  ),

  paste(
    "Fixed n:",
    FEATURE_N
  ),

  "",

  paste(
    "Mini-Batch fraction:",
    MINIBATCH_FRACTION
  ),

  "",

  "Methods:",

  paste(
    METHODS,
    collapse = ", "
  ),

  "",

  "Common refinement: Lloyd",

  "Timing units: milliseconds",

  "Peak memory: MiB",

  "Quality indices were not calculated during profiling.",

  "Parallel execution was not evaluated.",

  "Extended n experiment excludes exact Canopy because it creates an n x n distance matrix."

)


writeLines(

  settings_text,

  settings_file

)


# ======================================================================
# 34. PRINT SUMMARY
# ======================================================================

cat(
  "\n\n============================================================\n"
)

cat(
  "SCALABILITY EXPERIMENT COMPLETED\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nOutput directory:\n",
  OUTPUT_DIR,
  "\n"
)


cat(
  "\nObservation summary:\n\n"
)


print(
  observation_summary
)


cat(
  "\nFeature summary:\n\n"
)


print(
  feature_summary
)


if (
  RUN_EXTENDED &&
  nrow(extended_summary) > 0
) {

  cat(
    "\nExtended observation summary:\n\n"
  )

  print(
    extended_summary
  )
}


cat(
  "\nFiles generated:\n\n"
)


print(
  list.files(
    OUTPUT_DIR,
    full.names = FALSE
  )
)


cat(
  "\n============================================================\n"
)

cat(
  "EXPERIMENT DESIGN CHECK\n"
)

cat(
  "============================================================\n"

)

cat(
  "Farthest First             : YES\n"
)

cat(
  "Canopy                     : YES (common experiment)\n"
)

cat(
  "K-Means++                  : YES\n"
)

cat(
  "Forgy                      : YES\n"
)

cat(
  "ECKM                       : YES\n"
)

cat(
  "Mini-Batch K-Means++      : YES\n"
)

cat(
  "COSOHUC                    : YES\n"
)

cat(
  "Same X within condition    : YES\n"
)

cat(
  "Same k                     : YES\n"
)

cat(
  "Same Lloyd refinement      : YES\n"
)

cat(
  "Initialization timing      : YES\n"
)

cat(
  "Total timing               : YES\n"
)

cat(
  "Peak memory                : YES\n"
)

cat(
  "Parallel experiment        : NO\n"
)

cat(
  "============================================================\n"
)
