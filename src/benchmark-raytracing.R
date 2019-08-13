reprex::reprex({
  library(armacmp)
  # code from https://nextjournal.com/wolfv/how-fast-is-r-with-fastr-pythran
  # which in turn comes in part from http://www.tylermw.com/throwing-shade/
  # Author: Tyler Morgan-Wall

  # first the R version
  faster_bilinear <- function (Z, x0, y0){
    i = floor(x0)
    j = floor(y0)
    XT = (x0 - i)
    YT = (y0 - j)
    result = (1 - YT) * (1 - XT) * Z[i, j]
    nx = nrow(Z)
    ny = ncol(Z)
    if(i + 1 <= nx){
      result = result + (1-YT) * XT * Z[i + 1, j]
    }
    if(j + 1 <= ny){
      result = result + YT * (1-XT) * Z[i, j + 1]
    }
    if(i + 1 <= nx && j + 1 <= ny){
      result = result + YT * XT * Z[i + 1, j + 1]
    }
    result
  }

  bench_rays <- function() {
    volcanoshadow = matrix(1, ncol = ncol(volcano), nrow = nrow(volcano))

    volc = list(x=1:nrow(volcano), y=1:ncol(volcano), z=volcano)
    sunangle = 45 / 180*pi
    angle = -90 / 180 * pi
    diffangle = 90 / 180 * pi
    numberangles = 25
    anglebreaks = tan(seq(angle, diffangle, length.out = numberangles))
    maxdistance = floor(sqrt(ncol(volcano)^2 + nrow(volcano)^2))
    sinsun = sin(sunangle)
    cossun = cos(sunangle)
    for (i in 1:nrow(volcano)) {
      for (j in 1:ncol(volcano)) {
        vij = volcano[i, j]
        for (anglei in anglebreaks) {
          for (k in 1:maxdistance) {
            xcoord = (i + sinsun*k)
            ycoord = (j + cossun*k)
            if(xcoord > nrow(volcano) ||
               ycoord > ncol(volcano) ||
               xcoord < 0 || ycoord < 0) {
              break
            } else {
              tanangheight = vij + anglei * k
              if (all(c(volcano[ceiling(xcoord), ceiling(ycoord)],
                        volcano[floor(xcoord), ceiling(ycoord)],
                        volcano[ceiling(xcoord), floor(ycoord)],
                        volcano[floor(xcoord), floor(ycoord)]) < tanangheight)) next
              if (tanangheight < faster_bilinear(volcano, xcoord, ycoord)) {
                volcanoshadow[i, j] =  volcanoshadow[i, j] - 1 / length(anglebreaks)
                break
              }
            }
          }
        }
      }
    }
    return(volcanoshadow)
  }

  # Let's try to compile it to C++
  # We have to pass an initial matrix and the volcano matrix
  # because 1) we cannot yet init matrices and 2) the volcano matrix only exists in R
  # but apart from that everything can be translated with minimal changes (e.g. <- assignments and there is no `all` function)
  bench_rays_cpp <- compile(function(volcanoshadow_init, volcano) {
    volcanoshadow <- volcanoshadow_init
    faster_bilinear <- function(Z, x0 = type_scalar_numeric(), y0 = type_scalar_numeric()) {
      i <- floor(x0)
      j <- floor(y0)
      XT <- (x0 - i)
      YT <- (y0 - j)
      result <- (1 - YT) * (1 - XT) * Z[i, j]
      nx <- nrow(Z)
      ny <- ncol(Z)
      if (i + 1 <= nx) {
        result <- result + (1 - YT) * XT * Z[i + 1, j]
      }
      if (j + 1 <= ny) {
        result <- result + YT * (1 - XT) * Z[i, j + 1]
      }
      if (i + 1 <= nx && j + 1 <= ny) {
        result <- result + YT * XT * Z[i + 1, j + 1]
      }
      return(result, type = type_scalar_numeric())
    }

    sunangle <- 45 / 180 * pi
    angle <- -90 / 180 * pi
    diffangle <- 90 / 180 * pi
    numberangles <- 25
    anglebreaks <- tan(seq(angle, diffangle, numberangles))
    maxdistance <- floor(sqrt(ncol(volcano)^2 + nrow(volcano)^2))
    sinsun <- sin(sunangle)
    cossun <- cos(sunangle)
    for (i in seq_len(nrow(volcano))) {
      for (j in seq_len(ncol(volcano))) {
        vij <- volcano[i, j]
        for (anglei in anglebreaks) {
          for (k in seq_len(maxdistance)) {
            xcoord <- (i + sinsun * k)
            ycoord <- (j + cossun * k)
            if (xcoord > nrow(volcano) ||
                ycoord > ncol(volcano) ||
                xcoord < 0 || ycoord < 0) {
              break
            } else {
              tanangheight <- vij + anglei * k
              if (volcano[ceiling(xcoord), ceiling(ycoord)] < tanangheight &&
                  volcano[floor(xcoord), ceiling(ycoord)] < tanangheight &&
                  volcano[ceiling(xcoord), floor(ycoord)]< tanangheight &&
                  volcano[floor(xcoord), floor(ycoord)] < tanangheight) {
                next
              }
              if (tanangheight < faster_bilinear(volcano, xcoord, ycoord)) {
                volcanoshadow[i, j] <- volcanoshadow[i, j] - 1 / length(anglebreaks)
                break
              }
            }
          }
        }
      }
    }
    return(volcanoshadow)
  })


  all.equal(
    bench_rays(),
    bench_rays_cpp(matrix(1, ncol = ncol(volcano), nrow = nrow(volcano)), volcano)
  )

  system.time(bench_rays()) # takes too long for the microbenchmark

  microbenchmark::microbenchmark(
    cpp = bench_rays_cpp(matrix(1, ncol = ncol(volcano), nrow = nrow(volcano)), volcano)
  )
}, outfile = "benchmark-raytracing", venue = "R", show = FALSE)
