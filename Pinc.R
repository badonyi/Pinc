#!/usr/bin/env Rscript
# author: Mihaly Badonyi <badonyi@mpi-cbg.de>

if (!requireNamespace('jsonlite', quietly = TRUE)) {
  stop('The jsonlite package is required but not installed.')
}

CONTACT_RADIUS <- 12 # hard set for Pinc score definition, not user-configurable

# Functions ------------------------------------------------------------------ #
parse_command_args <- function() {
  # Parse and validate command line arguments.
  #
  # Returns:
  #   A named list with elements:
  #     - jsn_file: path to AlphaFold PAE matrix (JSON)
  #     - str_file: path to AlphaFold structure file (PDB/mmCIF)
  #     - usr_opt: character vector of user-defined output options
  #         always includes 'default', may also include 'all' and 'pairlist'
  #     - out_path: output file root name derived from structure file,
  #         prefixed with results folder if provided

  exit <- function() quit(status = 1, save = 'no')
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    cat(
      'Description:\n',
      '   R script to compute the Pinc score from an AlphaFold PAE\n',
      '   matrix (JSON) and the corresponding structure file (PDB/CIF).\n\n',
      'Usage:\n',
      '   Rscript Pinc.R <json_file> <structure_file> [options]\n\n',
      'Options:\n',
      '   None        Pinc score for all unique chain pairs (default)\n',
      '   --all       Full contact probability matrix\n',
      '   --pairlist  Residue pair list of non-zero contact probabilities\n',
      '   --results   Optional path to output folder\n\n',
      'Output format:\n',
      '   - default:  <structure_file_name>_Pinc.csv\n',
      '   - all:      <structure_file_name>_contact_probability.json\n',
      '   - pairlist: <structure_file_name>_pairlist.csv\n\n'
    )
    exit()
  }

  jsn_file <- args[1]
  if (tolower(tools::file_ext(jsn_file)) != 'json') {
    warning(
      'JSON file must have .json extension: ',
      shQuote(jsn_file),
      call. = FALSE
    )
    exit()
  }

  str_file <- args[2]
  if (!tolower(tools::file_ext(str_file)) %in% c('pdb', 'cif')) {
    warning(
      'Structure file must have .pdb or .cif extension: ',
      shQuote(str_file),
      call. = FALSE
    )
    exit()
  }

  missing_files <- c(jsn_file, str_file)[!file.exists(c(jsn_file, str_file))]
  if (length(missing_files)) {
    warning(
      'File(s) not found: ',
      paste(shQuote(missing_files), collapse = ', '),
      call. = FALSE
    )
    exit()
  }

  usr_opt <- 'default'
  results_dir <- getwd()

  if (length(args) >= 3) {
    opt <- args[-c(1, 2)]

    if ('--all' %in% opt) {
      usr_opt <- c(usr_opt, 'all')
    }

    if ('--pairlist' %in% opt) {
      usr_opt <- c(usr_opt, 'pairlist')
    }

    if ('--results' %in% opt) {
      i <- match('--results', opt)

      if (is.na(i) || i == length(opt)) {
        warning(
          'Option --results requires a folder path.',
          call. = FALSE
        )
        exit()
      }

      results_dir <- opt[i + 1]
    }
  }

  if (!dir.exists(results_dir)) {
    warning(
      'Results directory does not exist: ',
      shQuote(results_dir),
      call. = FALSE
    )
    exit()
  }

  if (file.access(results_dir, 2) != 0) {
    warning(
      'Results directory is not writable: ',
      shQuote(results_dir),
      call. = FALSE
    )
    exit()
  }

  out_path <- file.path(
    sub('[\\\\/]+$', '', results_dir),
    tools::file_path_sans_ext(basename(str_file))
  )

  list(
    jsn_file = jsn_file,
    str_file = str_file,
    usr_opt = usr_opt,
    out_path = out_path
  )
}

parse_pae <- function(jsn_file) {
  # Parse the PAE matrix from the AlphaFold JSON file.
  #
  # Args:
  #   jsn_file (character): path to JSON file containing PAE matrix
  #
  # Returns:
  #   An n x n matrix.

  jsn <- jsonlite::fromJSON(jsn_file)

  if (!is.null(jsn$pae)) {
    return(jsn$pae)
  }

  if (!is.null(jsn$predicted_aligned_error)) {
    return(jsn$predicted_aligned_error)
  }

  stop(
    'No PAE matrix found in JSON. Expected field ',
    '\'pae\' or \'predicted_aligned_error\'.'
  )
}

reduce_to_com <- function(coords) {
  # Collapse atomic coordinates to residue centres of mass.
  #
  # Args:
  #   coords (data.frame): structure data with columns:
  #     - c: chain
  #     - r: residue
  #     - e: element
  #     - x, y, z: coordinates
  #
  # Returns:
  #   A data.frame of residue centre of mass coordinates.

  xyz <- c('x', 'y', 'z')
  coords$m <- c(
    'S' = 32.0650,
    'P' = 30.9738,
    'O' = 15.9994,
    'N' = 14.0067
  )[coords$e]
  coords$m[is.na(coords$m)] <- 12.0107 # assume carbon

  # group for aggregation
  grp <- paste0(coords$c, '_', coords$r)
  grp <- factor(grp, levels = unique(grp))

  # compute centres of mass
  com <- aggregate(coords[, xyz] * coords$m, by = list(grp), FUN = sum)
  com[, xyz] <- com[, xyz] / tapply(coords$m, grp, sum)

  com
}

str_to_dist <- function(str_file) {
  # Parse a PDB/mmCIF file as a centre of mass distance matrix.
  #
  # Args:
  #   str_file (character): path to structure file
  #
  # Returns:
  #   An n x n matrix. Column & row names are chain IDs.

  lines <- readLines(str_file, warn = FALSE)
  atom_idx <- grep('^ATOM|^HETATM', lines)

  if (length(atom_idx) < 10) {
    stop('Non-structure file or insufficient number of ATOM records.')
  }

  if (any(grepl('_atom_site.group_PDB', lines))) {
    # CIF route
    start <- grep('_atom_site.group_PDB', lines)
    records <- trimws(lines[start:(atom_idx[1] - 1)])
    header <- c(
      c = '_atom_site.label_asym_id',
      r = '_atom_site.label_seq_id',
      e = '_atom_site.type_symbol',
      x = '_atom_site.Cartn_x',
      y = '_atom_site.Cartn_y',
      z = '_atom_site.Cartn_z'
    )

    keep <- records %in% header
    coords <- read.table(textConnection(lines[atom_idx]), na.strings = NULL)[
      keep
    ]
    colnames(coords) <- names(header)[match(records[keep], header)]

    # consecutive IDs for '.' placeholder in non-polymer entities
    coords$r <- ave(coords$r, coords$c, FUN = function(x) {
      m <- x == '.'
      if (any(m)) {
        x[m] <- 1:sum(m)
      }
      x
    })
  } else {
    # PDB route
    atom <- lines[atom_idx]
    coords <- data.frame(
      c = trimws(substr(atom, 22, 22)),
      r = trimws(substr(atom, 23, 26)),
      e = trimws(substr(atom, 77, 78)),
      x = as.numeric(substr(atom, 31, 38)),
      y = as.numeric(substr(atom, 39, 46)),
      z = as.numeric(substr(atom, 47, 54))
    )
  }

  # centre of mass distances (excluding hydrogens)
  com <- reduce_to_com(coords[coords$e != 'H', ])
  dist_mat <- as.matrix(dist(com[c('x', 'y', 'z')], upper = TRUE))

  # important: chain for columns, chain_tokenID for rows;
  # inherited by cont_asym and used in contact_pairlist()
  colnames(dist_mat) <- sub('_.*$', '', com[, 1])
  rownames(dist_mat) <- com[, 1]

  dist_mat
}

pae_to_prob <- function(pae_mat, dist_mat) {
  # Convert a PAE matrix into contact probabilities.
  #
  # Args:
  #   pae_mat (matrix): pairwise aligned error matrix
  #   dist_mat (matrix): residue distance matrix
  #
  # Returns:
  #   An n x n matrix.

  # sphere volume helper
  vol_sphere <- function(r) (4 / 3) * pi * (r^3)

  # sphere-sphere intersection volume (fixed reference sphere radius Rc,
  # uncertainty sphere radius Ru, centroid distance D)
  vol_intersect <- function(Ru, D) {
    vol_vec <- numeric(length(Ru))
    ok <- is.finite(Ru) & is.finite(D) & (Ru > 0)

    if (!any(ok)) {
      return(vol_vec)
    }

    rc <- CONTACT_RADIUS
    ru <- Ru[ok]
    dd <- D[ok]
    vol <- numeric(length(dd))

    # if disjoint, volume is zero (already initialised to zero)
    disjoint <- dd >= (rc + ru)

    # if one sphere fully contains the other, intersection is the smaller volume
    contain <- dd <= abs(rc - ru)
    if (any(contain)) {
      contain[is.na(contain)] <- FALSE
      vol[contain] <- vol_sphere(pmin(rc, ru[contain]))
    }

    # otherwise compute the standard lens formula analytically
    other <- !(contain | disjoint)
    if (any(other)) {
      other[is.na(other)] <- FALSE
      ddo <- dd[other]
      ruo <- ru[other]
      vol[other] <- pmax(
        0,
        pi *
          (rc + ruo - ddo)^2 *
          (ddo^2 + 2 * ddo * (ruo + rc) - 3 * (ruo - rc)^2) /
          (12 * ddo)
      )
    }

    vol_vec[ok] <- vol
    vol_vec
  }

  # flatten for vectorisation
  n <- nrow(pae_mat)
  unc_ij <- as.numeric(t(pae_mat))

  dist_vec <- as.numeric(t(dist_mat))
  vol_unc_ij <- vol_sphere(unc_ij)

  # raw contact probabilities
  p_ij <- vol_intersect(unc_ij, dist_vec) / vol_unc_ij

  # clamp to [0, 1] and guard against non-finite values
  p_ij[!is.finite(p_ij)] <- 0
  p_ij <- pmax(0, pmin(1, p_ij))

  # asymmetric matrix: row i, col j holds p(i->j); p(j->i) is the transpose
  cont_asym <- t(matrix(p_ij, nrow = n, ncol = n))
  diag(cont_asym) <- 1
  colnames(cont_asym) <- colnames(dist_mat)
  rownames(cont_asym) <- rownames(dist_mat)

  cont_asym
}

compute_pinc <- function(cont_asym, dist_mat) {
  # Compute the Pinc score between all non-redundant chain pairs.
  #
  # Args:
  #   cont_asym (matrix): asymmetric contact probability matrix
  #   dist_mat (matrix): residue distance matrix
  #
  # Returns:
  #   A data frame with columns chain1, chain2, Pinc1, Pinc2, Pinc.

  # list of chain indices for pairwise comparison
  chains <- colnames(dist_mat)
  uchain <- unique(chains)
  idx <- lapply(uchain, function(x) which(chains == x))
  names(idx) <- uchain

  # unique chain pair combinations
  pairs <- t(combn(uchain, 2))
  pinc1 <- numeric(nrow(pairs))
  pinc2 <- numeric(nrow(pairs))

  for (k in seq_len(nrow(pairs))) {
    i <- idx[[pairs[k, 1]]]
    j <- idx[[pairs[k, 2]]]

    dvec <- dist_mat[i, j, drop = FALSE]
    keep <- is.finite(dvec) & (dvec < CONTACT_RADIUS)

    if (any(keep)) {
      # Pinc1: chain1 residues carry the uncertainty sphere -> use cont_asym[j, i]
      # Pinc2: chain2 residues carry the uncertainty sphere -> use cont_asym[i, j]
      pinc1[k] <- mean(cont_asym[j, i][t(keep)])
      pinc2[k] <- mean(cont_asym[i, j][keep])
    }
  }

  pinc_df <- data.frame(
    chain1 = pairs[, 1],
    chain2 = pairs[, 2],
    Pinc1 = sprintf('%.4f', pinc1),
    Pinc2 = sprintf('%.4f', pinc2),
    Pinc = sprintf('%.4f', (pinc1 + pinc2) / 2)
  )

  # sort by Pinc score descending
  pinc_df[order(pinc_df$Pinc, decreasing = TRUE), ]
}

contact_pairlist <- function(cont_asym, dist_mat) {
  # Generate a list of contacting token pairs with their probabilities.
  #
  # Args:
  #   cont_asym (matrix): asymmetric contact probability matrix
  #   dist_mat (matrix): token distance matrix
  #
  # Returns:
  #   A data frame with columns: token1, token2, contact_p, distance,
  #   sorted by contact_p descending.

  # symmetrise
  cont_sym <- (cont_asym + t(cont_asym)) / 2

  # non-zero interchain contacts in the upper triangle
  keep <- upper.tri(cont_sym, diag = FALSE) &
    (cont_sym > 0) &
    (dist_mat < CONTACT_RADIUS) &
    (colnames(cont_sym)[row(cont_sym)] != colnames(cont_sym)[col(cont_sym)])

  # long format pair representation
  ij <- which(keep, arr.ind = TRUE)
  pair_df <- data.frame(
    token1 = rownames(cont_sym)[ij[, 1]],
    token2 = rownames(cont_sym)[ij[, 2]],
    contact_p = sprintf('%.4f', cont_sym[keep]),
    distance = sprintf('%.4f', dist_mat[keep])
  )

  # sort by contact probability descending
  pair_df[order(pair_df$contact_p, decreasing = TRUE), ]
}

# Process user input --------------------------------------------------------- #
args <- parse_command_args()

# parse the PAE matrix from the JSON file
pae_mat <- parse_pae(args$jsn_file)

# parse the distance matrix from the structure file
dist_mat <- str_to_dist(args$str_file)

# check for shape mismatch between PAE and distance matrices
if (!identical(dim(dist_mat), dim(pae_mat))) {
  stop(
    'Dimension mismatch between distance and contact matrices.\n',
    'Distance matrix dimensions: ',
    paste(dim(dist_mat), collapse = ' x '),
    '\n',
    'PAE matrix dimensions: ',
    paste(dim(pae_mat), collapse = ' x '),
    '\n',
    'Check that the structure file and PAE matrix correspond to the same model.\n'
  )
}

# convert PAE to contact probabilities
cont_asym <- pae_to_prob(pae_mat, dist_mat)

# compute and write Pinc scores (default; always produced)
write.csv(
  x = compute_pinc(cont_asym, dist_mat),
  file = paste0(args$out_path, '_Pinc.csv'),
  row.names = FALSE,
  quote = FALSE
)

# write full asymmetric contact probability matrix if requested
if ('all' %in% args$usr_opt) {
  jsonlite::write_json(
    x = list(
      token_chain_ids = colnames(cont_asym),
      contact_probability = cont_asym
    ),
    path = paste0(args$out_path, '_contact_probability.json'),
    digits = 4,
    pretty = TRUE
  )
}

# write residue-pair contact list if requested
if ('pairlist' %in% args$usr_opt) {
  write.csv(
    x = contact_pairlist(cont_asym, dist_mat),
    file = paste0(args$out_path, '_pairlist.csv'),
    row.names = FALSE,
    quote = FALSE
  )
}