####################################################################################################
# pseudo.presences.sampling
# Michel J.M. Bieleveld
# July 2017
####################################################################################################

#' Creates pseudo presences based on other presences to increase
#'
#' @param sp coordiantes of reponse points (2 column matrix)
#' @param resp Response Variable (monospecific) as vector, 1=pres, 0=true_abs, NA=no_info
#' @param env Explanatory Variable as rasterStack
#' @param strategy Strategy to be used to generate points, one of c("random", "sre", "similar", "dissimilar")
#' @param cell.distance Maximum cell.distance from presence point to be used for pseudo presence generation
#' @param np.points Number of pseudo presence points to generate
#' @param quant.SRE Quant parameter of sre function.
#' @return a SpatialPointsDataFrame object that will be given to BIOMOD_FormatingData function
'pseudo.presences.sampling' <- function(sp,resp, env, strategy='random', cell.distance=0, np.points=NULL, quant.SRE = 0){

  args <- .check.params.pseudo.presences.sampling(sp, resp, env, strategy, cell.distance, np.points, quant.SRE)

  sp <- args$sp
  env <- args$env
  strategy <- args$strategy
  cell.distance <- args$cell.distance
  np.points <- args$np.points
  quant.SRE <- args$quant.SRE

  mask = which(sp@data > 0)
  sp.presenceonly = as(sp[mask,],'SpatialPoints')
  sp.remainder = sp[-mask,]

  if((np.points <= 0 | cell.distance <= 0)){
    out <- sp[NULL,]
  } else {
    out <- switch(strategy,
                  random = .random.presences.sampling( sp.presenceonly, env, cell.distance, np.points ),
                  sre = .sre.presences.sampling(sp.presenceonly, env, np.points, quant.SRE),
                  similar = .fun.presences.sampling(sp.presenceonly, env, cell.distance, np.points, which.min),
                  dissimilar = .fun.presences.sampling(sp.presenceonly, env, cell.distance, np.points, which.max))

    out <- SpatialPointsDataFrame(out,data.frame(resp=rep(1,length(out))),proj4string=crs(env),match.ID = FALSE)
    if (length(out) < np.points) {
      warning(paste("Requested ",np.points," but only returned ",length(out),", try increasing cell.distance",sep=""))
    }
  }
  return(out)
}

.check.params.pseudo.presences.sampling <- function(sp, resp, env, strategy, cell.distance, np.points, quant.SRE){
  cat("\n   > Pseudo Presences Selection checkings...")

  # define here the implemented strategies
  availableStrategies <- c("random", "sre", "similar", "dissimilar")

  if (ncol(as.matrix(sp)) != 2)
  {
    stop("SP input must contain exacly two columns")
  }

  if(is.matrix(resp) | is.data.frame(resp)){
    if(ncol(resp) > 1){
      stop("You must give a monospecific response variable (1D object)")
    } else{
      resp <- as.numeric(resp[,1])
    }
  }


  # 1. sp input checking
  sp <- SpatialPointsDataFrame(as.matrix(sp), data.frame(resp),proj4string=crs(env))

  if(!(inherits(sp, 'SpatialPoints'))){
    stop("species input must be a SpatialPointsDataFrame object")
  }

  # 2. env input checking
  if(!inherits(env, 'Raster')){
    stop("Explanatory variables input must be a RasterStack object")
  }

  # 3. Strategy checking
  if( ! (strategy %in% availableStrategies) ){
    stop("You must select a valid strategy")
  }

  # 4. Np points checking
  if(is.null(np.points)){
    stop("You must give the number of pseudo absences you want")
  }

  # 5. Distances checking
  if(!is.null(cell.distance)){
    if(cell.distance <= 0){
      stop("Cell distance needs to be higher than 0")
    } else {
      dist <- area(env)@data@max * cell.distance
      cat("\n   > Approximate maximum distance",dist,"km\n")
    }
  }

  # 6. SRE quantil checking
  if(strategy == 'sre'){
    if( quant.SRE >= 0.5 | quant.SRE <0 ){
      stop("\n    ! SRE Quant should be a value between 0 and 0.5 ")
    }
  }

  # 7. return checked params
  return(list(sp = sp,
              env = env,
              strategy = strategy,
              cell.distance = cell.distance,
              np.points = np.points,
              quant.SRE = quant.SRE))
}

setGeneric(".get_near_col_row",
           function(sp,env,cell.distance){
             standardGeneric(".get_near_col_row")
           })

setMethod(".get_near_col_row",
          signature(sp="SpatialPoints", env = "BasicRaster",cell.distance="numeric"),
          function(sp,env,cell.distance)
          {
            env.mask <- calc(env, fun = sum)
            x = colFromX(env.mask, sp$x)
            y = rowFromY(env.mask, sp$y)
            xy.orig <- data.frame(x,y)
            x <- rep(x,each=(2*cell.distance+1)^2) + rep(seq(-cell.distance,cell.distance,1),length(x)*(2*cell.distance+1))
            y <- rep(rep(y,each=(2*cell.distance+1)) + seq(-cell.distance,cell.distance,1),each=(2*cell.distance+1))
            xy <- unique(data.frame(x,y))
            xy$cell <- env.mask[cellFromRowCol(env.mask,xy$y,xy$x)]
            xy <- xy[which(!is.na(xy$cell)),c("x","y")]
            all <- dplyr::anti_join(xy,xy.orig,by = c("x", "y"))
            rownames(all) <- NULL
            all
          })


setGeneric(".random.presences.sampling",
           function(sp,env,cell.distance,np.points){
             standardGeneric(".random.presences.sampling")
           })

setMethod(".random.presences.sampling",
          signature(sp="SpatialPoints", env = "BasicRaster",
                    cell.distance="numeric", np.points="numeric"),
          function(sp,env,cell.distance,np.points)
          {
            colrows <- .get_near_col_row(sp,env,cell.distance)
            colrows <- colrows[sample(nrow(colrows),min(np.points,nrow(colrows))),]
            points <- data.frame(x = xFromCol(env,colrows$x),
                                 y= yFromRow(env,colrows$y))
            coordinates(points) <- ~ x + y
            crs(points) <- crs(env)
            points
          })

setGeneric(".sre.presences.sampling",
           function(sp,env,np.points,...){
             standardGeneric(".sre.presences.sampling")
           })

setMethod(".sre.presences.sampling",
          signature(sp="SpatialPoints", env = "BasicRaster",
                    np.points="numeric"),
          function(sp,env,np.points,quant=0.025)
          {
            myResp <- reclassify(subset(env,1,drop=TRUE), c(-Inf,Inf,0))
            myResp[sp] <- 1
            mask.in <- sre(Response = myResp, Explanatory = env, NewData = env, quant)
            mask.in[sp] <- 0
            points <- as.data.frame(rasterToPoints(mask.in, function(x) x ==1)[,c("x","y")])
            points <- sample(points,replace=FALSE,size=min(length(points),np.points))
            coordinates(points) <- ~ x + y
            crs(points) <- crs(env)
            points
          })


setGeneric(".fun.presences.sampling",
           function(sp,env,cell.distance,np.points,...){
             standardGeneric(".fun.presences.sampling")
           })

setMethod(".fun.presences.sampling",
          signature(sp="SpatialPoints", env = "BasicRaster",
                    cell.distance="numeric", np.points="numeric"),
          function(sp,env,cell.distance,np.points, f=which.max)
          {
            expls <- scale(env)
            colrows <- .get_near_col_row(sp,expls,cell.distance)
            points <- data.frame(x = xFromCol(env,colrows$x),
                                 y = yFromRow(env,colrows$y))
            coordinates(points) <- ~ x + y
            crs(points) <- crs(env)
            values.presences <- raster::extract(expls,sp)
            values.candidates <- raster::extract(expls,points)
            selected <- NULL
            for (i in seq(1,min(np.points,nrow(values.candidates)-1)))
            {
              d <- f(apply(values.candidates,1,function(x) {sum(colSums(values.presences - x)^2)} ))
              values.presences <- rbind(values.presences,values.candidates[d,])
              values.candidates <- values.candidates[-c(d),]
              selected <- append(selected, d)
            }
            points[selected,]
          })






