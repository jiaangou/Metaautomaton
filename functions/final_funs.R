# Required packages
library(dplyr)
library(tidyr)
library(infotheo)

#############################################
# Model components ------------------------
#############################################

# 1. Competition: update densities with stochasticity
competition <- function(N, alpha, r) {
  no.spp <- length(r)
  comp <- aperm(apply(N, c(1,2), function(x) x %*% alpha), c(2,3,1))
  N_new <- r * N / (1 + comp)
  N_out <- aperm(apply(N_new, c(1,2), function(x) rpois(no.spp, x)), c(2,3,1))
  return(N_out)
}

# 2. Emigration: individuals leave patches
emigration <- function(N, disp) {
  out <- aperm(apply(N, c(1,2), function(x) rbinom(length(x), size = x, prob = disp)), c(2,3,1))
  return(out)
}

# 3. Viability: survival of emigrants
viability <- function(N, cost) {
  out <- aperm(apply(N, c(1,2), function(x) rbinom(length(x), size = x, prob = 1 - cost)), c(2,3,1))
  return(out)
}

# 4. Invasion fitness: used to weight immigrants if desired
invasion_fitness <- function(N, r, A) {
  W <- r / (1 + N %*% A)
  return(W)
}

# 5. Neighborhood array: for computing local weights and rearranging outputs
neighborhood_array <- function(mat, L = NULL, rearrange = FALSE) {
  if (is.null(L)) L <- dim(mat)[1]
  minus_one <- c(tail(1:L, 1), head(1:L, -1))
  plus_one  <- c(tail(1:L, -1), head(1:L, 1))
  if (!rearrange) {
    out <- array(c(mat[, minus_one],
                   mat[, plus_one],
                   mat[minus_one, ],
                   mat[plus_one, ]),
                 dim = c(L, L, 4))
  } else {
    out <- array(c(mat[, plus_one, 1],
                   mat[, minus_one, 2],
                   mat[plus_one, , 3],
                   mat[minus_one, , 4]),
                 dim = c(L, L, 4))
  }
  return(out)
}

# 6. Weighted frequency: calculates neighbor weights based on density
weighted_freq <- function(N, delta) {
  w <- (N + 1)^delta
  w_freq <- w / sum(w)
  return(w_freq)
}

# 7. Immigration iterator: sample immigrants for each patch
immi_iterator <- function(w, surv) {
  out <- array(0, dim = dim(w))
  rows <- dim(surv)[1]
  cols <- dim(surv)[2]
  for (r in 1:rows) {
    for (c in 1:cols) {
      out[r, c, ] <- rmultinom(n = 1, size = surv[r, c], prob = w[r, c, ])
    }
  }
  return(out)
}

# 8. Immigration function: computes arrivals based on neighborhood weights
immi_function <- function(N_array, survivors, delta, no.spp = 2) {
  weights <- lapply(1:no.spp, function(s) {
    nh <- neighborhood_array(N_array[,,s])
    apply(nh, c(1,2), function(x) weighted_freq(x, delta = delta)) %>%
      aperm(c(2,3,1))
  })
  out <- lapply(1:no.spp, function(s) {
    immi_iterator(w = weights[[s]], surv = survivors[,,s]) %>%
      neighborhood_array(rearrange = TRUE) %>%
      apply(c(1,2), sum)
  })
  return(out)
}

#############################################
# Functions for computing metrics ----------------------
#############################################

# 9. Coarse graining: reduces densities to state codes ("11", "10", "01", or "00")
coarse_graining <- function(N_array, states = c(4,3)){

  #4 states: extinction, sp1 alone, co-occurrence, sp2 alone
  if(states == 4){
    out <- N_array%>%
      apply(c(1,2), function(x)(x>0)%>%
              as.integer%>%
              paste(collapse = ''))

  }
  #3 states: extinction, sp1 dominant, sp2 dominant
  else if(states == 3){
    out <- N_array%>%
      apply(c(1,2), function(x)
        ifelse(sum(x) == 0, 0, which.max(x)))

  }else(stop("select either 3 or 4 for the number of coarse grained states"))

  return(out)
}

# 10. Nearest neighbor coordinates
neighborhood_coord <- function(matrix = NULL, Lx = NULL, Ly = NULL, r = 1) {
  # Determine dimensions
  if (is.null(Lx) && is.null(Ly) && !is.null(matrix)) {
    Lx <- dim(matrix)[1]
    Ly <- dim(matrix)[2]
  } else if (is.null(Lx) || is.null(Ly)) {
    stop("Either 'matrix' or both 'Lx' and 'Ly' must be provided.")
  }

  #output_type <- match.arg(output_type) # Ensure valid output_type

  # Helper function
  wrap <- function(x, L) ((x - 1) %% L) + 1

  # Generate offsets for von Neumann neighborhood of radius r
  get_von_neumann_offsets <- function(r) {
    offsets <- expand.grid(dx = -r:r, dy = -r:r)
    offsets <- offsets[abs(offsets$dx) + abs(offsets$dy) <= r, ]
    offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0), ] # Exclude the center
    return(offsets)
  }

  offsets <- get_von_neumann_offsets(r)

  lattice <- matrix(1:(Lx * Ly), nrow = Lx)
  num_neighbors <- nrow(offsets)
  out <- data.frame(local_ID = integer(Lx * Ly * num_neighbors), neighborID = integer(Lx * Ly * num_neighbors))
  row_index <- 1
  for (i in 1:Lx) {
    for (j in 1:Ly) {
      current <- lattice[i, j]

      for (k in 1:num_neighbors) {
        dx <- offsets$dx[k]
        dy <- offsets$dy[k]
        neighbor_i <- wrap(i + dx, Lx)
        neighbor_j <- wrap(j + dy, Ly)
        neighbor <- lattice[neighbor_i, neighbor_j]

        out[row_index, ] <- c(current, neighbor)
        row_index <- row_index + 1
      }
    }
  }
  return(out)
}

# 11. Entropy function
entropy <- function(N, log_b = 2){
  p <- N/sum(N)
  if(any(p == 0)){
    hx <- 0
  }else{
    hx <- -sum(p*log(p, base = log_b))
  }

  return(hx)
}

# 12. Niche overlap: a measure of similarity between two species’ "coefficient-weighted" density distributions
niche <- function(N_array, A = NULL, weight = FALSE){

  #if weight = TRUE, weight densities by competition coefficients
  if(!is.null(A) & weight == TRUE){
    N <- lapply(1:nrow(A), function(x)t(A[x,] * t(apply(N_array, 3, function(y)y)))) #weight densities by coefficients
  }else{
    N <- lapply(1:dim(N_array)[3], function(x)N_array%>%
                  apply(3, function(x)x)) #collapse into a Spp x Patch matrix
  }

  #Nornmalize densities
  norm_N <- lapply(N, function(x)x/sum(x))

  #Macroscopic (realized) coefficients
  #Intra
  intra <- lapply(1:length(norm_N), function(x)norm_N[[x]][,x]^2%>%
                    sum)%>%
    unlist()

  #Inter
  inter <- lapply(norm_N, function(x)apply(x, 1, prod)%>%sum)%>%
    unlist

  #Ratio
  ratios <- inter/intra

  #Geometric mean (niche overlap, rho)
  rho <- log(ratios)%>%
    mean()%>%
    exp()

  #return(rho)
  #Return
  out <- c(ratios, rho)
  names(out) <- c(paste0('a', 1:length(ratios)), 'rho')
  return(out)

}


# 13. Compute metrics: regional N, Hx, MI, overlap
compute_metrics <- function(N_array, A = NULL, weight = FALSE, neighbors = NULL){

  #Check N is correct dimensions (3: M x N x Spp)
  #dims <- dim(N)
  if (length(dim(N_array)) != 3) stop("N should have 2 spatial and 1 species dimensions")

  #If neighbors are not supplied, compute it from N (default with r = 1)
  if(is.null(neighbors)){
    neighbors <- neighborhood_coord(N_array[,,1])
  }

  #Compute coarse-grained states to be used to compute metrics
  coarse_grainings <- N_array%>%
    coarse_graining(states = 4)

  states <- data.frame(local = coarse_grainings[neighbors$local_ID],
                       neighbor = coarse_grainings[neighbors$neighborID])%>%
    mutate(joint = paste0(local, neighbor))


  #1. Regional densities
  rN <- apply(N_array, 3, sum)


  #2. Compute individual, joint, and mutual entropies
  entropies <- apply(states, 2, function(x)x%>%
                       table%>%
                       entropy(log_b = 2))%>%
    t()%>%
    as.data.frame()%>%
    setNames(c('Hx','Hy','Hxy'))%>%
    mutate(MI = Hx + Hy - Hxy)%>%
    mutate(NMI = ifelse(Hx >0, MI/Hx, 0))%>%
    select(MI, NMI)

  entropies[entropies<10^(-5)] <- 0 #anything less than 5 decimal places is set to 0


  #3. Co-occurrence probability
  cooccur_p <- sum(states$local == '11') / nrow(states)


  #4. Niche overlap
  overlap <- niche(N_array = N_array, A = A, weight = weight)


  #Combine and return
  out <- c('rN' = rN, entropies,
           overlap,
           'cooccurrence' = cooccur_p)%>%unlist

  return(out)

}



#############################################
# Main simulator function: metaautomaton
#############################################
metaautomaton <- function(N, r = 1.5, aii = 1, aij, K = 50,
                       d_disp, delta, cost_c, ext_p = 0,
                       shuffle_p = 0, time = 1,
                       weight_fitness = FALSE, weight_overlap = FALSE,
                       metrics = FALSE, record_every = NULL) {

  # Setup
  L <- dim(N)[1]               # spatial extent (assumed square)
  no.spp <- length(r)          # number of species
  A <- matrix(aij, nrow = no.spp, ncol = no.spp)
  diag(A) <- aii
  A <- A / (K / (r - 1))       # scale by carrying capacity

  # If record_every is provided, preallocate recording arrays:
  if (!is.null(record_every)) {
    n_stored <- floor(time / record_every)
    j <- 1
    if (metrics) {
      metrics_out <- matrix(NA, nrow = n_stored, ncol = no.spp + 7)
      colnames(metrics_out) <- c(paste0('r_N', 1:no.spp),
                                 'MI', 'NMI', 'a1', 'a2', 'rho', 'cooccurrence','time')
    } else {
      N_t <- array(0, dim = c(L, L, no.spp, n_stored))
    }
  }

  # Pre-compute cost (constant over time) and neighborhood coordinates (for metrics)
  cost <- 1 - exp(-cost_c * delta)
  coords <- neighborhood_coord(Lx = L, Ly = L)

  # Main simulation loop
  for (i in 1:time) {
    N0 <- N  # store state before updating

    # 1. Competition
    N_prime <- competition(N, alpha = A, r = r)
    # 2. Emigration
    E <- emigration(N_prime, disp = d_disp)
    # 3. Survival
    S <- viability(E, cost = cost)
    # 4a. Compute weights if using fitness weighting; otherwise, use N_prime
    if (weight_fitness) {
      W <- aperm(apply(N_prime, c(1,2), function(x)
        invasion_fitness(x, r = r, A = A) - 1),
        c(2,3,1))
    } else {
      W <- N_prime
    }
    # 4b. Immigration (using vectorized immi_function)
    I_mig <- immi_function(W, survivors = S, delta = delta, no.spp = no.spp)
    I_array <- array(unlist(I_mig), dim = c(L, L, no.spp))
    # 5. Disturbance: random extinction events
    D_mat <- array(rep(matrix(rbinom(L^2, 1, prob = 1 - ext_p), nrow = L, ncol = L),
                       no.spp),
                   dim = c(L, L, no.spp))
    # Update state
    N <- (N_prime - E + I_array) * D_mat
    # 6. Optional spatial shuffle
    if (rbinom(1, 1, prob = shuffle_p) == 1) {
      shuffled_indices <- sample(L^2)
      N <- simplify2array(lapply(1:no.spp, function(s)
        matrix(N[,,s][shuffled_indices], nrow = L, ncol = L)))
    }

    # Check for extinction (if any species has zero total abundance)
    extinction <- any(apply(N, 3, sum) == 0)

    # Record state/metrics at intervals if record_every is specified
    if (!is.null(record_every) && (i %% record_every == 0)) {
      if (metrics) {
        if (!extinction) {
          metrics_out[j, ] <- c(compute_metrics(N_array = N, A = A, weight = weight_overlap, neighbors = coords), i)
          j <- j + 1
        }
      } else {
        N_t[,,, j] <- N
        j <- j + 1
      }
    }

    # If extinction occurs, record the previous state/metrics (if not already recorded) and break
    if (extinction) {
      last_recorded_time <- i - 1
      if (!is.null(record_every) && (last_recorded_time %% record_every != 0)) {
        if (metrics) {
          metrics_out[j, ] <- c(compute_metrics(N_array = N0, A = A, weight = TRUE, neighbors = coords), last_recorded_time)

          #metrics_out <- metrics_out[!is.na(metrics_out[, 1]), ] #remove empty row

        } else {
          N_t[,,, j] <- N0
        }
      }
      break
    }
  } # end main loop

  # If record_every is NULL, then record only the final state/metrics
  if (is.null(record_every)) {
    if (metrics) {

      metrics_out <- c(compute_metrics(N_array = N, A = A, weight = TRUE, neighbors = coords), i)
      names(metrics_out) <- c(paste0('r_N', 1:no.spp),
                              'MI', 'NMI', 'a1', 'a2', 'rho', 'cooccurrence','time')

    } else {
      N_t <- N
    }
  }

  if (metrics) {
    if(!is.null(record_every)){
      metrics_out <- metrics_out[!is.na(metrics_out[, 1]), ] #remove empty rows
    }
    return(metrics_out)
  } else {
    return(N_t)
  }
}




#Tests
# L <- 10 #Landscape size (64 x 64)
# no.spp <- 2
# r <- rep(1.5, no.spp)
# N0 <- array(50, dim = c(L, L, no.spp)) #initial condition
#
#
# metaautomaton(N = N0, r = r, aij = 2, K = 50,
#                    d_disp = 0.5, delta = 0.1,
#                    cost_c = 1, ext_p = 0.05,
#                    shuffle_p = 0, time = 300,
#                    weight_fitness = FALSE, weight_overlap = FALSE,
#                    metrics = TRUE, record_every = 2)
#
# #Visualize
# lapply(1:dim(N)[4], function(x)coarse_graining(N[,,,x], states = 4)%>%
#          reshape2::melt()%>%
#          setNames(c('row','col', 'State')))%>%
#   bind_rows(.id = 'time')%>%
#   mutate(time = as.integer(time))%>%
#   ggplot(aes(x = col, y = row))+
#   geom_tile(aes(fill = State))+
#   facet_wrap(~time)
#





