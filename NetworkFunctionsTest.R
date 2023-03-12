source("NetworkFunctions.R") # loading the functions from the NetworkFunctions.R file

########## - Testing the algorithms - ##########

## Testing function
test_fn <- function(mat, exp_autoregsize = 0, exp_numFFL = 0, exp_pathlength = 0, print_results = FALSE) {

  autoreg <- find_autoreg(mat)
  num_autoreg <- length(autoreg)
  FFL <- find_feed_forw(mat)
  num_FFL <- nrow(FFL)
  longestpath <- find_longest_path(mat)
  
  if (print_results) {
    print("Testing:")
    print(mat)
    print("Autoregulations in the matrix")
    print(autoreg)
    print(num_autoreg)
    print("Feed forward loops in the matrix")
    print(FFL)
    print(num_FFL)
    print("Longest path in the matrix")
    print(longestpath)
  }
  
  if (num_autoreg != exp_autoregsize) {
    print(paste0("ERROR autoreg: expected ", exp_autoregsize, " but got ", num_autoreg))
  }
  
  if (num_FFL != exp_numFFL) {
    print(paste0("ERROR FFL: expected ", exp_numFFL, " but got ", num_FFL))
  }
  
  if (longestpath != exp_pathlength) {
    print(paste0("ERROR FFL: expected ", exp_pathlength, " but got ", longestpath))
  }
  
  if (num_autoreg == exp_autoregsize & num_FFL == exp_numFFL & longestpath == exp_pathlength) {
    print("Successful test.")
  }
}

## Testing with some matrices

# 4x4 matrix, all 1
print("Test: 4x4 matrix, all 1")
mat1 <- matrix(rep(1,16), nrow=4, ncol=4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
test_fn(mat1, exp_autoregsize = 4, exp_pathlength = 3)

# 4x4 matrix, all 0
print("Test: 4x4 matrix, all 0")
mat0 <- matrix(rep(0,16), nrow=4, ncol=4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
test_fn(mat0)

# 4x4 matrix, indentity matrix
print("Test: 4x4 matrix, indentity matrix")
mat4i <- diag(4)
rownames(mat4i) <- LETTERS[1:4]
colnames(mat4i) <- LETTERS[1:4]
test_fn(mat4i, exp_autoregsize = 4)

# 3x3 matrix, FFL with autoregs
print("Test: 3x3 matrix, FFL with autoregs")
fflar3 <- rbind(c(1, 1, 1), c(0, 1, 1), c(0, 0, 1))
rownames(fflar3) <- LETTERS[1:3]
colnames(fflar3) <- LETTERS[1:3]
test_fn(fflar3, exp_autoregsize = 3, exp_numFFL = 1, exp_pathlength = 2)

# 3x3 matrix, FFL without autoregs
print("Test: 3x3 matrix, FFL without autoregs")
ffl3 <- rbind(c(1, 1, 1), c(0, 1, 1), c(0, 0, 1)) - diag(3)
rownames(ffl3) <- LETTERS[1:3]
colnames(ffl3) <- LETTERS[1:3]
test_fn(ffl3, exp_numFFL = 1, exp_pathlength = 2)

# 3x3 matrix, almost FFL but also B->A
print("Test: 3x3 matrix, almost FFL but also B->A")
fflba3 <- rbind(c(1, 1, 1), c(1, 1, 1), c(0, 0, 1)) - diag(3)
rownames(fflba3) <- LETTERS[1:3]
colnames(fflba3) <- LETTERS[1:3]
test_fn(fflba3, exp_pathlength = 2)

# 3x3 matrix, almost FFL but also C->A
print("Test: 3x3 matrix, almost FFL but also C->A")
fflca3 <- rbind(c(1, 1, 1), c(0, 1, 1), c(1, 0, 1)) - diag(3)
rownames(fflca3) <- LETTERS[1:3]
colnames(fflca3) <- LETTERS[1:3]
test_fn(fflba3, exp_pathlength = 2)

# 3x3 matrix, almost FFL but also C->B
print("Test: 3x3 matrix, almost FFL but also C->B")
fflcb3 <- rbind(c(1, 1, 1), c(0, 1, 1), c(0, 1, 1)) - diag(3)
rownames(fflcb3) <- LETTERS[1:3]
colnames(fflcb3) <- LETTERS[1:3]
test_fn(fflba3, exp_pathlength = 2)

# 3x3 matrix, A->B->C->A
print("Test: 3x3 matrix, A->B->C->A")
abca3 <- rbind(c(0, 1, 0), c(0, 0, 1), c(1, 0, 0))
rownames(abca3) <- LETTERS[1:3]
colnames(abca3) <- LETTERS[1:3]
test_fn(abca3, exp_pathlength = 2)

# 6x6 matrix, A->B->C->D->E->F, F->all
# FFL: FAB, FBC, FCD
print("Test: 6x6 matrix, A->B->C->D->E->F, F->all")
mat6 <- rbind(c(0,1,0,0,0,0),
              c(0,0,1,0,0,0),
              c(0,0,0,1,0,0),
              c(0,0,0,0,1,0),
              c(0,0,0,0,0,1),
              c(1,1,1,1,1,1))
rownames(mat6) <- LETTERS[1:6]
colnames(mat6) <- LETTERS[1:6]
test_fn(mat6, exp_autoregsize = 1, exp_pathlength = 5, exp_numFFL = 3)

# 5x5 matrix, 1 FFL (A, B, C), D<->E and autoreg D, E
print("Test: 5x5 matrix, 1 FFL (A, B, C), D<->E and autoreg D, E")
mat5 <- rbind(c(0,1,1,0,0),
              c(0,0,1,0,0),
              c(0,0,0,0,0),
              c(0,0,0,1,1),
              c(0,0,0,1,1))
rownames(mat5) <- LETTERS[1:5]
colnames(mat5) <- LETTERS[1:5]
test_fn(mat5, exp_autoregsize = 2, exp_numFFL = 1, exp_pathlength = 2)

# 8x8 matrix, A->A, B->B, C->all, all->C, D->E->F->G, H->H
print("Test: 8x8 matrix, A->A, B->B, C->all, all->C, D->E->F->G, H->H")
mat8 <- rbind(c(1,0,1,0,0,0,0,0),
              c(0,1,1,0,0,0,0,0),
              c(1,1,1,1,1,1,1,1),
              c(0,0,1,0,1,0,0,0),
              c(0,0,1,0,0,1,0,0),
              c(0,0,1,0,0,0,1,0),
              c(0,0,1,0,0,0,0,0),
              c(0,0,1,0,0,0,0,1))
rownames(mat8) <- LETTERS[1:8]
colnames(mat8) <- LETTERS[1:8]
test_fn(mat8, exp_autoregsize = 4, exp_pathlength = 5)

# 6x6 matrix, A->B, C->D, E->F
print("Test: 6x6 matrix, A->B, C->D, E->F")
pair6 <- rbind(c(0,1,0,0,0,0),
               c(0,0,0,0,0,0),
               c(0,0,0,1,0,0),
               c(0,0,0,0,0,0),
               c(0,0,0,0,0,1),
               c(0,0,0,0,0,0))
rownames(pair6) <- LETTERS[1:6]
colnames(pair6) <- LETTERS[1:6]
test_fn(pair6, exp_pathlength = 1)

# Test that should fail, to see the errors:
# 3x3 matrix, FFL with autoregs
fflar3 <- rbind(c(1, 1, 1), c(0, 1, 1), c(0, 0, 1))
rownames(fflar3) <- LETTERS[1:3]
colnames(fflar3) <- LETTERS[1:3]
# Right answer:
# test_fn(fflar3, exp_autoregsize = 3, exp_numFFL = 1, exp_pathlength = 2)
# Wrong answer:
print("Test: The following test is designed to fail, to check that the three errors work correctly")
test_fn(fflar3, exp_autoregsize = 0, exp_numFFL = 0, exp_pathlength = 0)
