## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)
names = c("White Li","Great Meng","Pure-lake Tao","Easy-life Bai")
born = c(701,NA,365,772)
DATASET = data.frame(name = names, born = born)

usethis::use_data(DATASET, overwrite = TRUE)
