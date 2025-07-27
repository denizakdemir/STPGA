#' @title Performance Optimization Functions
#' @description Functions to improve computational performance for large datasets
#' @author Deniz Akdemir

# =============================================================================
# MEMORY AND PERFORMANCE UTILITIES
# =============================================================================

#' Add progress bar for long-running computations
#' 
#' @param total Total number of iterations
#' @param title Title for the progress bar
#' @param verbose Whether to show progress bar
#' @return Progress bar object or NULL
#' @export
create_progress_bar <- function(total, title = "Processing", verbose = TRUE) {
  if (!verbose || total <= 1) {
    return(NULL)
  }
  
  if (requireNamespace("progress", quietly = TRUE)) {
    progress::progress_bar$new(
      format = paste0(title, " [:bar] :percent in :elapsed"),
      total = total,
      clear = FALSE,
      width = 60
    )
  } else {
    # Fallback to base R
    structure(
      list(
        total = total,
        current = 0,
        start_time = Sys.time(),
        title = title
      ),
      class = "simple_progress"
    )
  }
}

#' Update progress bar
#' 
#' @param pb Progress bar object
#' @param increment How much to increment (default: 1)
#' @export
update_progress <- function(pb, increment = 1) {
  if (is.null(pb)) return(invisible(NULL))
  
  if (inherits(pb, "Progress")) {
    pb$tick(increment)
  } else if (inherits(pb, "simple_progress")) {
    pb$current <- pb$current + increment
    if (pb$current %% max(1, floor(pb$total / 10)) == 0) {
      pct <- round(100 * pb$current / pb$total)
      elapsed <- difftime(Sys.time(), pb$start_time, units = "secs")
      cat(sprintf("\r%s: %d%% (%d/%d) in %.1fs", 
                  pb$title, pct, pb$current, pb$total, elapsed))
      if (pb$current >= pb$total) cat("\n")
      flush.console()
    }
  }
  
  invisible(NULL)
}

#' Finish progress bar
#' 
#' @param pb Progress bar object
#' @export
finish_progress <- function(pb) {
  if (is.null(pb)) return(invisible(NULL))
  
  if (inherits(pb, "simple_progress") && pb$current < pb$total) {
    elapsed <- difftime(Sys.time(), pb$start_time, units = "secs")
    cat(sprintf("\r%s: 100%% (%d/%d) in %.1fs\n", 
                pb$title, pb$total, pb$total, elapsed))
    flush.console()
  }
  
  invisible(NULL)
}

# =============================================================================
# CHUNKED MATRIX OPERATIONS
# =============================================================================

#' Compute matrix operations in chunks for memory efficiency
#' 
#' @param X Matrix to process
#' @param FUN Function to apply to each chunk
#' @param chunk_size Maximum chunk size (default: 1000)
#' @param ... Additional arguments to FUN
#' @return Combined result from all chunks
#' @export
chunked_matrix_operation <- function(X, FUN, chunk_size = 1000, ...) {
  n_rows <- nrow(X)
  
  if (n_rows <= chunk_size) {
    return(FUN(X, ...))
  }
  
  # Calculate chunks
  n_chunks <- ceiling(n_rows / chunk_size)
  chunk_starts <- seq(1, n_rows, by = chunk_size)
  chunk_ends <- pmin(chunk_starts + chunk_size - 1, n_rows)
  
  results <- vector("list", n_chunks)
  
  # Process chunks with progress bar
  pb <- create_progress_bar(n_chunks, "Processing chunks")
  
  for (i in 1:n_chunks) {
    chunk_data <- X[chunk_starts[i]:chunk_ends[i], , drop = FALSE]
    results[[i]] <- FUN(chunk_data, ...)
    update_progress(pb)
  }
  
  finish_progress(pb)
  
  # Combine results
  if (is.matrix(results[[1]])) {
    do.call(rbind, results)
  } else if (is.numeric(results[[1]])) {
    unlist(results)
  } else {
    results
  }
}

#' Efficient block matrix multiplication
#' 
#' @param X First matrix
#' @param Y Second matrix
#' @param block_size Block size for computation (default: 500)
#' @return Matrix product X %*% Y
#' @export
block_matrix_mult <- function(X, Y, block_size = 500) {
  if (ncol(X) != nrow(Y)) {
    stop("Incompatible matrix dimensions for multiplication")
  }
  
  n <- nrow(X)
  p <- ncol(Y)
  k <- ncol(X)
  
  # If matrices are small, use regular multiplication
  if (n <= block_size && p <= block_size && k <= block_size) {
    return(X %*% Y)
  }
  
  # Initialize result matrix
  result <- matrix(0, nrow = n, ncol = p)
  
  # Block multiplication
  n_blocks_i <- ceiling(n / block_size)
  n_blocks_j <- ceiling(p / block_size)
  n_blocks_k <- ceiling(k / block_size)
  
  total_blocks <- n_blocks_i * n_blocks_j * n_blocks_k
  pb <- create_progress_bar(total_blocks, "Block multiplication")
  
  for (i in 1:n_blocks_i) {
    i_start <- (i - 1) * block_size + 1
    i_end <- min(i * block_size, n)
    
    for (j in 1:n_blocks_j) {
      j_start <- (j - 1) * block_size + 1
      j_end <- min(j * block_size, p)
      
      for (l in 1:n_blocks_k) {
        l_start <- (l - 1) * block_size + 1
        l_end <- min(l * block_size, k)
        
        # Compute block product
        X_block <- X[i_start:i_end, l_start:l_end, drop = FALSE]
        Y_block <- Y[l_start:l_end, j_start:j_end, drop = FALSE]
        
        result[i_start:i_end, j_start:j_end] <- 
          result[i_start:i_end, j_start:j_end] + X_block %*% Y_block
        
        update_progress(pb)
      }
    }
  }
  
  finish_progress(pb)
  return(result)
}

#' Memory-efficient crossprod computation
#' 
#' @param X Matrix
#' @param Y Optional second matrix (default: X)
#' @param chunk_size Chunk size for computation (default: 1000)
#' @return t(X) %*% Y
#' @export
chunked_crossprod <- function(X, Y = NULL, chunk_size = 1000) {
  if (is.null(Y)) {
    Y <- X
  }
  
  if (nrow(X) != nrow(Y)) {
    stop("X and Y must have the same number of rows")
  }
  
  n <- nrow(X)
  
  if (n <= chunk_size) {
    return(crossprod(X, Y))
  }
  
  # Initialize result
  result <- matrix(0, nrow = ncol(X), ncol = ncol(Y))
  
  # Process in chunks
  n_chunks <- ceiling(n / chunk_size)
  chunk_starts <- seq(1, n, by = chunk_size)
  chunk_ends <- pmin(chunk_starts + chunk_size - 1, n)
  
  pb <- create_progress_bar(n_chunks, "Chunked crossprod")
  
  for (i in 1:n_chunks) {
    idx <- chunk_starts[i]:chunk_ends[i]
    X_chunk <- X[idx, , drop = FALSE]
    Y_chunk <- Y[idx, , drop = FALSE]
    
    result <- result + crossprod(X_chunk, Y_chunk)
    update_progress(pb)
  }
  
  finish_progress(pb)
  return(result)
}

# =============================================================================
# PARALLEL PROCESSING UTILITIES
# =============================================================================

#' Parallel matrix operations with proper error handling
#' 
#' @param X List of matrices or data to process
#' @param FUN Function to apply
#' @param mc.cores Number of cores to use
#' @param ... Additional arguments to FUN
#' @return List of results
#' @export
parallel_apply <- function(X, FUN, mc.cores = parallel::detectCores() - 1, ...) {
  if (length(X) <= 1 || mc.cores <= 1) {
    # Sequential processing
    pb <- create_progress_bar(length(X), "Sequential processing")
    results <- lapply(X, function(x) {
      result <- FUN(x, ...)
      update_progress(pb)
      result
    })
    finish_progress(pb)
    return(results)
  }
  
  # Parallel processing
  if (.Platform$OS.type == "windows") {
    # Use parLapply on Windows
    cl <- parallel::makeCluster(mc.cores)
    on.exit(parallel::stopCluster(cl))
    
    # Export necessary objects
    parallel::clusterEvalQ(cl, library(STPGA))
    
    results <- parallel::parLapply(cl, X, FUN, ...)
  } else {
    # Use mclapply on Unix-like systems
    results <- parallel::mclapply(X, FUN, ..., 
                                  mc.cores = mc.cores, 
                                  mc.cleanup = TRUE,
                                  mc.silent = FALSE)
  }
  
  # Check for errors
  errors <- sapply(results, function(x) inherits(x, "try-error"))
  if (any(errors)) {
    warning("Some parallel computations failed: ", sum(errors), " errors")
  }
  
  return(results)
}

#' Get optimal number of cores for computation
#' 
#' @param data_size Approximate size of data
#' @param operation_type Type of operation ("cpu", "memory", "mixed")
#' @return Recommended number of cores
#' @export
get_optimal_cores <- function(data_size, operation_type = "mixed") {
  max_cores <- parallel::detectCores()
  
  # Conservative approach: leave one core for system
  available_cores <- max(1, max_cores - 1)
  
  # Adjust based on operation type and data size
  if (operation_type == "memory") {
    # Memory-bound operations benefit less from many cores
    recommended <- min(available_cores, 4)
  } else if (operation_type == "cpu") {
    # CPU-bound operations can use more cores
    recommended <- available_cores
  } else {  # mixed
    # Balance between CPU and memory
    if (data_size < 1e6) {
      recommended <- min(available_cores, 2)
    } else if (data_size < 1e8) {
      recommended <- min(available_cores, 4)
    } else {
      recommended <- min(available_cores, 8)
    }
  }
  
  return(recommended)
}

# =============================================================================
# CACHING UTILITIES
# =============================================================================

#' Simple hash-based cache for computation results
#' 
#' @param key Cache key (will be hashed)
#' @param value Value to cache (only if not retrieving)
#' @param cache_env Environment to store cache (default: private environment)
#' @return Cached value if exists, NULL otherwise
#' @export
simple_cache <- function(key, value = NULL, cache_env = NULL) {
  if (is.null(cache_env)) {
    if (!exists(".stpga_cache", envir = .GlobalEnv)) {
      assign(".stpga_cache", new.env(hash = TRUE), envir = .GlobalEnv)
    }
    cache_env <- get(".stpga_cache", envir = .GlobalEnv)
  }
  
  # Create hash key
  hash_key <- digest::digest(key, algo = "md5")
  
  if (is.null(value)) {
    # Retrieve from cache
    if (exists(hash_key, envir = cache_env)) {
      return(get(hash_key, envir = cache_env))
    } else {
      return(NULL)
    }
  } else {
    # Store in cache
    assign(hash_key, value, envir = cache_env)
    return(value)
  }
}

#' Clear computation cache
#' 
#' @param cache_env Cache environment to clear
#' @export
clear_cache <- function(cache_env = NULL) {
  if (is.null(cache_env)) {
    if (exists(".stpga_cache", envir = .GlobalEnv)) {
      rm(list = ls(get(".stpga_cache", envir = .GlobalEnv)), 
         envir = get(".stpga_cache", envir = .GlobalEnv))
    }
  } else {
    rm(list = ls(cache_env), envir = cache_env)
  }
  
  invisible(NULL)
}

#' Get cache statistics
#' 
#' @param cache_env Cache environment
#' @return List with cache statistics
#' @export
cache_stats <- function(cache_env = NULL) {
  if (is.null(cache_env)) {
    if (!exists(".stpga_cache", envir = .GlobalEnv)) {
      return(list(size = 0, memory_usage = 0))
    }
    cache_env <- get(".stpga_cache", envir = .GlobalEnv)
  }
  
  cache_objects <- ls(cache_env)
  cache_size <- length(cache_objects)
  
  # Estimate memory usage
  memory_usage <- sum(sapply(cache_objects, function(obj) {
    tryCatch({
      object.size(get(obj, envir = cache_env))
    }, error = function(e) 0)
  }))
  
  list(
    size = cache_size,
    memory_usage = as.numeric(memory_usage),
    memory_mb = round(as.numeric(memory_usage) / 1024^2, 2)
  )
}