library(dplyr)
library(adept)
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
walk_out = segmentWalking(xyz,
                          xyz.fs = x.fs,
                          template = template)
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

system.time({
  vec_result = adept:::similarityMatrix_vectorized(
    x.cut.seq = x.cut.seq,
    x.smoothed = x.smoothed,
    template.scaled = template.scaled,
    similarity.measure = similarity.measure,
    x.cut.vl = x.cut.vl,
    template.vl = template.vl,
    x.cut.margin = x.cut.margin)
})
index = vec_result$index
vec_result = vec_result$similarity

vec_result_reordered = lapply(seq_along(x.cut.seq), function(x) {
  matrix(nrow = length(template.scaled), ncol = ncol(vec_result[[1]]))
})
for (i in seq_along(x.cut.seq)) {
  for (itemp in seq_along(template.scaled)) {
    vec_result_reordered[[i]][itemp,] = vec_result[[itemp]][i,]
  }
}

system.time({
  orig_result = adept:::similarityMatrix_original(
    x.cut.seq = x.cut.seq,
    x.smoothed = x.smoothed,
    template.scaled = template.scaled,
    similarity.measure = similarity.measure,
    x.cut.vl = x.cut.vl,
    template.vl = template.vl,
    x.cut.margin = x.cut.margin)
})
all.equal(orig_result[[1]], vec_result_reordered[[1]])
all.equal(orig_result[1:16], vec_result_reordered[1:16])

r = vec_result_reordered[[17]]
all_bad = apply(is.na(r), 2, all)
head(which(all_bad))
vec_r = r[, !all_bad]

r = orig_result[[17]]
all_bad = apply(is.na(r), 2, all)
head(which(all_bad))
orig_r = r[, !all_bad]



