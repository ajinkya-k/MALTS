

d1 <- function(Xc, Mc, disc=FALSE) {
 Dc <- make_dist_array(Xc, disc)

 Dct<- colSums((Dc * Mc) ^ 2)
 return(as.matrix(Dct))
}


n <- 1000
p <- 20

Mc <- sample(1:10, size = 20, replace = T)


nt <- 2000
XX <- matrix(rnorm(n*p), nrow = n)


dd <- cross_dist(XX, XX, Mc)

dd2 <- d1(XX, Mc)

dd3 <- cross_dist_rfast(XX, XX, Mc)

dd4 <- cross_dist_rfast2(XX, XX, Mc)
all.equal(dd, dd2)

mbench <- microbenchmark::microbenchmark(
 "MALTS" = d1(XX, Mc),
 # "proxy" = cross_dist(XX, XX, Mc),
 "rfast" = cross_dist_rfast(XX, XX, Mc),
 "rfast2" = cross_dist_rfast2(XX, XX, Mc)
)


ggplot() +
 geom_boxplot(aes(x = mbench$expr, y = mbench$time, fill = mbench$expr)) + 
 scale_y_log10() +
 theme_bw() +
 khroma::scale_fill_discreterainbow()


mbench2 <- 

bench::mark(
 "MALTS" = d1(XX, Mc),
 "proxy" = cross_dist(XX, XX, Mc),
 "rfast" = cross_dist_rfast(XX, XX, Mc),
 "rfast2" = cross_dist_rfast2(XX, XX, Mc),
 iterations = 1000,
 check = all.equal
)

times <- mbench2$time

names(mbench2) <- c("MALTS", "proxy", )

library(bench)
library(tidyr)

mbench2 %>%
 tidyr::unnest(c(time, memory_alloc)) %>%
 filter(gc == "none") %>%
 mutate(expression = as.character(expression)) %>%
 ggplot(aes(x = mem_alloc, y = time, color = expression)) +
 geom_point() +
 scale_color_bench_expr(scales::brewer_pal(type = "qual", palette = 3))



mbench2$time |> str()
mbench2$gc |> str()
mbench2$mem_alloc |> str()


mbench2 %>% 
 unnest(c(time, mem_alloc))




Am <- rno

library(microbenchmark)

mbenc_dmult <- 
 microbenchmark(
  "swp" = sweep(XX, 2L, Mc, `*`),
  "mat" = XX %*% diag(Mc),
  times = 10000L
 )


boxplot(log10(time)~expr, data = mbenc_dmult)

all.equal(
 sweep(XX, 2L, Mc, `*`),
 XX %*% diag(Mc)
)


run_check <- function() {
 Mc <- sample(1:10, size = 20, replace = T)
 XX <- matrix(rnorm(n*p), nrow = n)


 ddm <- d1(XX, Mc)
 ddf <- cross_dist_rfast(XX, XX, Mc)
 return(all.equal(ddm, ddf))
}


all(replicate(1000, run_check()))





cross_dist_rfast(model.matrix(~., Xd_C)[, -1], model.matrix(~., Xd_C)[, -1], c(1/sqrt(2),1/sqrt(2),sqrt(2),sqrt(2)))


d1(Xd_C, c(1,2), disc=TRUE )


Xd5 <- data.frame(
 f1 = as.factor(sample(letters[1:20], 1000, replace = T)),
 f2 = as.factor(sample(letters[2:25], 1000, replace = T)),
 f3 = as.factor(sample(letters[8:16], 1000, replace = T))
)


#TODO: Potential issue: if too many factor levels this could be really slow
hamming <- function(Xd, w) {
 n_levels <- unname(sapply(Xd, nlevels))

 if (any(n_levels < 2)) stop(paste0(names(Xd)[which(n_levels == 1)], " has only one level"))

 w_exp <- rep(w/sqrt(2), n_levels)

 # print(w_exp)
 Xm <- model.matrix(~., Xd, contrasts.arg = lapply(Xd, contrasts, contrasts=FALSE))[,-1]


 #TODO: Use tcrossprod(Xm * w_, Xm)
 # where w_sq_by_2 = rep(w/sqrt(2), n_levels)
 d_mat <- cross_dist_rfast(Xm, Xm, w_exp)

 return(d_mat)
}


hamming2 <- function(Xd, w) {
 n_levels <- unname(sapply(Xd, nlevels))

 if (any(n_levels < 2)) stop(paste0(names(Xd)[which(n_levels == 1)], " has only one level"))

 w_exp <- rep(w/sqrt(2), n_levels)

 # print(w_exp)
 Xm <- model.matrix(~., Xd, contrasts.arg = lapply(Xd, contrasts, contrasts=FALSE))[,-1]

 Xmw <- reweight_mat(Xm, w_exp)
 d_mat <- tcrossprod(Xmw, Xmw)

 return(d_mat)
}

hamming3 <- function(Xd, w) {
 n_levels <- unname(sapply(Xd, nlevels))

 if (any(n_levels < 2)) stop(paste0(names(Xd)[which(n_levels == 1)], " has only one level"))

 w_exp <- rep(w^2/2, n_levels)

 # print(w_exp)
 Xm <- model.matrix(~., Xd, contrasts.arg = lapply(Xd, contrasts, contrasts=FALSE))[,-1]


 Xmp <- 1-Xm
 Xmw <- reweight_mat(Xm, w_exp)
 dhalf <- tcrossprod(Xmw, Xmp)
 d_mat <- setop(dhalf, "+", t(dhalf), rowwise = TRUE)

 return(d_mat)
}



all.equal(d1(Xd5, c(1,2,4), disc=TRUE ),
hamming(Xd5, c(1,2,4))
)


 bench::mark(
  "malts" = d1(Xd5, c(1,2,4), disc=TRUE ),
  "mine" = hamming(Xd5, c(1,2,4)),
  iterations = 1000L
 )


plot(
 microbenchmark(
  "malts" = d1(Xd5, c(1,2,4), disc=TRUE ),
  "mine" = hamming(Xd5, c(1,2,4)),
  times = 100L
 )
)



run_check2 <- function() {
  Xd5 <- data.frame(
    f1 = as.factor(sample(letters[1:20], 1000, replace = T)),
    f2 = as.factor(sample(letters[2:25], 1000, replace = T)),
    f3 = as.factor(sample(letters[8:16], 1000, replace = T))
  )

  w = runif(3, 1, 5)
  all.equal(
    d1(Xd5, c(1,2,4), disc=TRUE ),
    hamming(Xd5, c(1,2,4))
  )

}

all(replicate(100, run_check2()))





Xd4 <- data.frame(
 as.factor(sample(letters[1:20], 1000, replace = T)),
 as.factor(sample(letters[2:25], 1000, replace = T)),
 as.factor(sample(letters[8:16], 1000, replace = T)),
 as.factor(sample(letters[1:20], 1000, replace = T)),
 as.factor(sample(letters[2:25], 1000, replace = T)),
 as.factor(sample(letters[8:16], 1000, replace = T)),
 as.factor(sample(letters[1:20], 1000, replace = T)),
 as.factor(sample(letters[2:25], 1000, replace = T)),
 as.factor(sample(letters[8:16], 1000, replace = T)),
 as.factor(sample(letters[1:20], 1000, replace = T)),
 as.factor(sample(letters[2:25], 1000, replace = T)),
 as.factor(sample(letters[8:16], 1000, replace = T)),
 as.factor(sample(letters[1:20], 1000, replace = T)),
 as.factor(sample(letters[2:25], 1000, replace = T)),
 as.factor(sample(letters[8:16], 1000, replace = T))
)


outer(Xd5[,1], Xd5[,1],function(x,y) as.numeric(x == y))



Rfast::Outer(Xd5[,1], Xd5[,1],function(x,y) as.numeric(x == y))


Xd4_small <- data.frame(
  f1 = as.factor(sample(letters[1:2], 5, replace = T)), 
  f2 = as.factor(sample(letters[1:2], 5, replace = T))
)

Xd4b <- data.frame(
  f1 = sample(0:1, 3, replace = T), 
  f2 = sample(0:1, 3, replace = T)
)

Reduce(
  `+`,
lapply(
  Xd4,
  function(x) {
    outer(x,x, `==`)
  }
)
)



w_ = runif(15, 1, 5)

 
bench::mark(
  "mr" = purrr::reduce(
  lapply(
    1:ncol(Xd5),
    FUN = function(i, x,y, w) {
      outer(x[,i],y[,i], `!=`) * (w[i])^2
    },
    x = Xd5, y = Xd5, w = w_
  ),
  `+`
  ),
  "hm" = hamming(Xd5, w_),
  iterations = 100
)
Reduce(
  `+`,
  lapply(
    1:ncol(Xd4),
    FUN = function(i, x,y, w) {
      outer(x[,i],y[,i], `!=`) * (w[i])^2
    },
    x = Xd4, y = Xd4, w = c(1,2)
  )
) 


purrr::reduce(
  lapply(
    1:ncol(Xd4),
    FUN = function(i, x,y, w) {
      outer(x[,i],y[,i], `!=`) * (w[i])^2
    },
    x = Xd4, y = Xd4, w = c(1,2)
  ),
  `+`
) 




hamming(Xd4, c(1,2))



cmark <- 
bench::mark(
  "d1" = d1(Xd4, w_, TRUE),
  "hm" = hamming(Xd4, w_),
  "hm-cr" = unname(hamming3(Xd4, w_)) ,
  iterations = 1000
)

cmark

plot(cmark, type = "boxplot")
hamming(Xd4, c(1,2)) %>% str
hamming2(Xd4, c(1,2))
hamming3(Xd4, c(1,2)) %>% unname


d1(Xd4, c(1,2), TRUE)



all.equal(hamming(Xd4, c(1,2)), hamming3(Xd4, c(1,2)) %>% unname)







rand_mat <- matrix(sample(0:1, 30, replace = T), ncol = 5)


cp <- tcrossprod(1 - rand_mat, rand_mat)


make_model_mat <- function(Xd) {
  model.matrix(~., Xd, contrasts.arg = lapply(Xd, contrasts, contrasts=FALSE))[,-1]
}


bench::press(
  rows = c(1,2)
  {
    bench::mark(
      "mm" = make_model_mat(df),
      "fd" = fastDummies::dummy_cols(df, remove_selected_columns=TRUE)
    )
  }
)



bench::mark(
      "mm" = unname(make_model_mat(Xd4)),
      "fd" = unname(data.matrix(fastDummies::dummy_cols(Xd4, remove_selected_columns=TRUE)))
    ) %>% plot(type = "boxplot")

make_model_mat(Xd4_small) %>% unname

fastDummies::dummy_columns(Xd4_small, remove_selected_columns=TRUE) %>%  data.matrix(., ) %>% unname 


unname(make_model_mat(Xd4_small))
unname(data.matrix(fastDummies::dummy_cols(Xd4, remove_selected_columns=TRUE)))





 model.matrix(~., Xd4_small, contrasts.arg = lapply(Xd, contrasts, contrasts=FALSE))[,-1]

fm <-
as.formula(
  paste0(c("~0", names(Xd4_small)), collapse = "+")
)


 model.matrix(~fm1+fm2, data = Xd_c)


Xd_c <- lapply(Xd4_small, as.character)
