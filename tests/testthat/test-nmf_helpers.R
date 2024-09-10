example(transfer_labels, echo = TRUE)
test_that("outputted NMF factors are all nonnegative", {
  source_factors <- t(res$nmf$h)
  expect_equal(sum(source_factors >= 0), length(source_factors))
})

test_that("predicted probabilities are in [0,1]",{


  vis_anno_target <- res$targets
  multinom_mod <- res$multinom
  projections <- reducedDim(vis_anno_target, "nmf_projections")

  probs <- predict(multinom_mod, newdata=projections, type='probs',
                                    na.action=na.exclude)
  expect_equal(sum(probs >= 0 & probs <= 1), length(probs))
  expect_equal(as.vector(rowSums(probs)), rep(1, nrow(probs)))
})
