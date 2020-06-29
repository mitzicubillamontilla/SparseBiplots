context("Ridge HJ-Biplot")

# Generate Some Data in R
V1 <- rnorm(500, 0, 1)
V2 <- rnorm(500, 0, 1)
V3 <- rnorm(500,0,1)
V4 <-
  -0.1*V1 +
  0.1*V2 +
  rnorm(500,0,1)

# Create imput matrix
X <- cbind(V1, V2, V3, V1, V2, V3, V1, V2, V3, V4, V4, V4)

X <- X +
  matrix(rnorm(length(X),0,1),
         ncol = ncol(X),
         nrow = nrow(X)
         )

#--------------------------------------------------------------------
# Test HJBiplot output
#--------------------------------------------------------------------

out <- Ridge_HJBiplot(X, Lambda = 0.3) # HJ-Biplot

test_that("1. Test output", {
  expect_identical(
    names(out),
    c("eigenvalues", "explvar", "loadings", "coord_ind", "coord_var")
    )
  })

test_that("2. Test for loadings", {
  expect_identical(
    typeof(out$loadings),
    "double"
  )
  # expect_identical(
  #   class(out$loadings),
  #   "matrix"
  # )
  expect_identical(
    nrow(out$loadings),
    ncol((X))
    )
  expect_identical(
    ncol(out$loadings),
    ncol(X)
    )
  })

test_that("3. Test for row coordinates", {
  expect_identical(
    typeof(out$coord_ind),
    "double"
  )
  # expect_identical(
  #   class(out$coord_ind),
  #   "matrix"
  # )
  expect_identical(
    nrow(out$coord_ind),
    nrow((X))
  )
  expect_identical(
    ncol(out$coord_ind),
    ncol(X)
    )
  })

test_that("4. Test for column coordinates", {
  expect_identical(
    typeof(out$coord_var),
    "double"
  )
  # expect_identical(
  #   class(out$coord_var),
  #   "matrix"
  # )
  expect_identical(
    nrow(out$coord_var),
    ncol((X))
  )
  expect_identical(
    ncol(out$coord_var),
    ncol(X)
    )
  })

test_that("5. Test for eigenvalues", {
  expect_identical(
    typeof(out$eigenvalues),
    "double"
  )
  expect_identical(
    class(out$eigenvalues),
    "numeric"
  )
  expect_identical(
    length(out$eigenvalues),
    ncol(X)
    )
  })

test_that("6. Test for explained variance", {
  expect_identical(
    typeof(out$explvar),
    "double"
  )
  expect_identical(
    class(out$explvar),
    "numeric"
  )
  expect_identical(
    length(out$explvar),
    ncol(X)
    )
  })
