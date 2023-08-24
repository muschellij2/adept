library(Matrix)
library(dplyr)
devtools::load_all()
# muschellij2/adept@nhanes
library(adeptdata)
all_wrist_templates = adeptdata::stride_template$left_wrist
template_list = do.call(rbind, all_wrist_templates)
template_list = apply(template_list, 1, identity, simplify = FALSE)

template = template_list[1:2]
# devtools::load_all()
options(digits.secs = 3)
sample_rate = 10L

# read in some test data
data = readr::read_csv(
  "~/Dropbox/Projects/nhanes_80hz/data/csv/pax_h/73557.csv.gz",
  # "https://github.com/martakarass/adept/files/12423699/test_data.csv.gz",
  # n_max = 1e6)
)
xyz = data %>% select(X, Y, Z) %>% as.matrix()
x.fs = 80L
xyz <- as.matrix(xyz)
# walk_out = segmentWalking(xyz,
#                           xyz.fs = x.fs,
#                           template = template)
xyzptr <- as.data.frame(cbind(xyz, pracma::cart2sph(xyz)))
x <- xyzptr[, 6]

pattern.dur.seq = seq(0.5, 4, length.out = 30)
similarity.measure = "cor"
similarity.measure.thresh = -2
x.adept.ma.W = 0.2
finetune = "maxima"
finetune.maxima.ma.W = NULL
finetune.maxima.nbh.W = 0.6
run.parallel = FALSE
run.parallel.cores = 1L
x.cut = TRUE
x.cut.vl = 6000L


template.vl <- pattern.dur.seq * x.fs
template.vl <- sort(unique(round(template.vl)))
template.vl.max <- max(template.vl)
template.vl.min <- min(template.vl)
if (!is.list(template)) template <- list(template)
template.scaled <- scaleTemplate(template, template.vl)


x.smoothed <- adept:::get.x.smoothed(x = x, W = x.adept.ma.W,
                                     x.fs = x.fs)

finetune.maxima.x <- x
finetune.maxima.nbh.vl <- round(finetune.maxima.nbh.W *
                                  x.fs)
if (!x.cut) x.cut.vl <- length(x)
x.cut.margin <- template.vl.max - 1
x.cut.seq <- seq(1, to = length(x), by = x.cut.vl)
template.idx.mat.i <- NULL

##############
# New Code
##############
all_x = xcut_to_data_matrix(
  x.cut.seq = x.cut.seq,
  x.smoothed = x.smoothed,
  x.cut.vl = x.cut.vl,
  template.vl = template.vl,
  x.cut.margin = x.cut.margin)
nc_all_x = ncol(all_x)
not_na_all_x = !is.na(all_x)
all_x[!not_na_all_x] = 0
# all_x = Matrix(all_x, sparse = TRUE)
not_na_all_x = Matrix(not_na_all_x, sparse = TRUE)


# template = lapply(1:10, function(r) template)
x_stats_length_template = function(
    length_template,
    nc_all_x,
    all_x,
    not_na_all_x,
    similarity.measure) {
  max_first_index = nc_all_x - length_template + 1

  ones = c(rep(1, length_template), rep(0, nc_all_x - length_template))
  ones = toeplitz(ones)
  ones[upper.tri(ones)] = 0
  ones = ones[,1:max_first_index]
  ones <- Matrix::Matrix(ones, sparse = TRUE)
  n_x = not_na_all_x %*% ones
  if (similarity.measure == "cor") {
    sum_x = all_x %*% ones
    mean_x = sum_x / n_x

    sum_x2 = all_x^2 %*% ones
  } else {
    mean_x = sum_x2 = NULL
  }
  list(
    # sum_x = sum_x,
    n_x = n_x,
    mean_x = mean_x,
    sum_x2 = sum_x2
  )
}

template_set = template.scaled[[1]]
ti = template_set[[1]]

result = pbapply::pblapply(template.scaled, function(template_set) {

  length_template = length(template_set[[1]])
  out = x_stats_length_template(
    length_template = length_template,
    nc_all_x = nc_all_x,
    all_x = all_x,
    not_na_all_x = not_na_all_x,
    similarity.measure = similarity.measure)

  en_minus_1 = length_template - 1
  if (similarity.measure == "cov") {
    denominator = (out$n_x - 1)
  } else {
    denominator = (sqrt(out$sum_x2 - out$n_x * out$mean_x^2) *
                     sqrt(en_minus_1))
  }
  denominator[out$n_x <= 0] = NA_real_
  denominator[!is.finite(denominator)] = NA_real_

  itemp = 1
  for (itemp in seq_along(template_set)) {
    ti = template_set[[itemp]]
    stopifnot(mean(ti) <= 1e-10)
    stopifnot( abs(sd(ti) - 1) <= 1e-10)

    system.time({
      max_first_index = nc_all_x - length_template + 1

      tt = c(ti, rep(0, nc_all_x - length_template))
      M = toeplitz(tt)
      M[upper.tri(M)] = 0
      M = M[,1:max_first_index]
      M <- Matrix::Matrix(M, sparse = TRUE)
      sum_xy = all_x %*% M
      rm(M)
      output = sum_xy / denominator
      output
    })
  }
})

