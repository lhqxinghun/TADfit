genRectFromTADCoor <- function(TADsCoor,binstart,binend){
  a <- matrix(ncol=4)
  a_add <- matrix(ncol=4,nrow=1)
  #print(a[1,2])
  
  length_heatmap <-binend-binstart+1
  ####-----------------select TAD——————————
  TADsCoor_useful <- NULL
  for(i in 1:nrow(TADsCoor))
  {
    if((TADsCoor[i,1])>=binstart | (TADsCoor[i,2]>binstart))
    {
      if(((TADsCoor[i,1]<binend) | (TADsCoor[i,2])<=binend)){
        current <- TADsCoor[i,2] : TADsCoor[i,1]
        TADsCoor_useful <- rbind(TADsCoor_useful,range(current))
      }
    }
  }
  TADsCoor_useful <- as.data.frame(TADsCoor_useful)
  #-------------------------------------\
  print(TADsCoor_useful)
  
  
  if((TADsCoor_useful[1,1]<binstart)&&(TADsCoor[1,2]>binstart)){
    a[1,1] <- 1
    a[1,2] <- TADsCoor_useful[1,2] -binstart +1
    a[1,3] <- length_heatmap
    a[1,4] <- length_heatmap-a[1,2]+1
  }
  if((TADsCoor_useful[1,1]>=binstart)&&(TADsCoor_useful[1,2]<=binend)){
    a[1,1] <- TADsCoor_useful[1,1]-binstart+1
    a[1,2] <- TADsCoor_useful[1,2]-binstart+1
    a[1,3] <- length_heatmap+1-a[1,1]
    a[1,4] <- length_heatmap+1-a[1,2]
  }
  if((TADsCoor_useful[1,1]<binend) && (TADsCoor_useful[1,2]>binend)){
    a[1,1] <- TADsCoor_useful[1,1]-binstart+1
    a[1,2] <- length_heatmap
    a[1,3] <- length_heatmap-a[1,1]+1
    a[1,4] <- length_heatmap-a[1,2]+1
  }
  
  #a[1,1] <- 1
  #a[1,2] <- TADsCoor_useful[1,2]-binstart+1
  #a[1,3] <- length_heatmap
  #a[1,4] <- length_heatmap+binstart-TADsCoor_useful[1,2]
  
  if(nrow(TADsCoor_useful >=2 )){
    for(j in 2:nrow(TADsCoor_useful))
    {
      
      if((TADsCoor_useful[j,1]<binstart)&&(TADsCoor[j,2]>binstart)){
        a_add[1,1] <- 1
        a_add[1,2] <- TADsCoor_useful[j,2] -binstart +1
        a_add[1,3] <- length_heatmap
        a_add[1,4] <- length_heatmap-a_add[1,2]+1
        a <- rbind(a,a_add)
        
      }
      if((TADsCoor_useful[j,1]>=binstart)&&(TADsCoor_useful[j,2]<=binend)){
        a_add[1,1] <- TADsCoor_useful[j,1]-binstart+1
        a_add[1,2] <- TADsCoor_useful[j,2]-binstart+1
        a_add[1,3] <- length_heatmap+1-a_add[1,1]
        a_add[1,4] <- length_heatmap+1-a_add[1,2]
        a <- rbind(a,a_add)
        
      }
      if((TADsCoor_useful[j,1]<binend) && (TADsCoor_useful[j,2]>binend)){
        a_add[1,1] <- TADsCoor_useful[j,1]-binstart+1
        a_add[1,2] <- length_heatmap
        a_add[1,3] <- length_heatmap-a_add[1,1]+1
        a_add[1,4] <- length_heatmap-a_add[1,2]+1
        a <- rbind(a,a_add)
        
      }
      
      
    }
  }
  
  #__________________________________
  
  return(a)
  
}


#-------------------------------
TADs_dup_rm <- function(TADsCoor){
  result2temp <- TADsCoor[!duplicated(TADsCoor),]
  #TADsCoor.Out <- result2temp[order(result2temp[,1]),]
  #return(TADsCoor.Out)
  return(result2temp)
}

## Function from original "comparator.R" script by Aaron Lun
##
altmatch <- function(x1, x2, y1, y2) {
  # Quick dirty matcher.
  pair.x <- paste(x1, x2, sep=".")
  pair.y <- paste(y1, y2, sep=".")
  matched <- match(pair.x, pair.y)
  keep <- !is.na(matched)
  return(list(unique1=which(!keep),
              unique2=(1:length(pair.y))[-matched[keep]],
              match1=which(keep),
              match2=matched[keep]))
}




#####################################################
#####################################################


spawnTADs_withNesting_dualOutput<- function(ntads, min.tad=10, max.tad=30, nestedLevels=3, nestedProb=0.3, nestedSizeCap=NULL, overlapProbs=NULL) 
{
  total_ref<-0
  if (!is.null(nestedSizeCap)) {
    if (nestedSizeCap<=max.tad) {
      stop("nestedSizeCap should be larger than max.tad")
    }
  }
  ntads <- as.integer(ntads)
  mysize <- as.integer(runif(ntads, min.tad, max.tad))
  tads <- list()
  tads_intervals<-NULL
  tads_intervals_one <-NULL
  mysize_ref<-mysize
  
  #--------------------- first layer ----------------
  last.launch <- 0L
  pre_existing_tads<-length(tads)
  for (i in 1:length(mysize_ref)) { 
    current <- last.launch + mysize_ref[i]:1   ## %%% catch max min of this to return the TADs start end positions
    tads_intervals<-rbind(tads_intervals, range(current))
    tads[[i+pre_existing_tads]] <- t(combn(current, 2L))
    colnames(tads[[i+pre_existing_tads]]) <- c("anchor", "target")
    last.launch <- last.launch + mysize_ref[i]
  }
  print("first layer")
  print(length(tads))
  
  
  #-----------------------------------second layer--------------------------
  size_one <- length(tads_intervals)/2
  #print(size_two)
  borderSkippingProbs <- overlapProbs
  #remove_size_two <- floor(borderSkippingProbs*(size_one))
  remove_size_two <- round(borderSkippingProbs*(size_one))
  #print(remove_size_two)
  every_space_two <- round(size_one/remove_size_two)
  #print(every_space_two)
  #current_BordersRemoval_two <- array(1:(remove_size_two))
  current_BordersRemoval_two <- rep(0, remove_size_two)
  

  for (i in 1:(remove_size_two)){
    top_space <- sample(2:(every_space_two-1),1)
    current_BordersRemoval_two[i] <- (every_space_two*(i-1) + top_space)
  }
  #print(current_BordersRemoval_two)
  for(imerge in 1:(size_one-1)){
    if (imerge %in% current_BordersRemoval_two) {
      random_overlap_number_left <- as.integer(runif(1,min=1,max=1.1))
      random_overlap_number_right <- as.integer(runif(1,min=1,max=1.1))
      current <-  (tads_intervals[length(mysize_ref)+imerge]):tads_intervals[imerge - random_overlap_number_left]
      tads_intervals_one<-rbind(tads_intervals_one, range(current))
      current <- (tads_intervals[length(mysize_ref)+imerge+random_overlap_number_right]):tads_intervals[imerge]
      tads_intervals_one<-rbind(tads_intervals_one, range(current))
    }
  }
  print(tads_intervals_one)
  j <- 1
  
  
  for(i in 1:(size_one)){
    if((tads_intervals[i+length(tads_intervals)/2] <= tads_intervals_one[length(tads_intervals_one)]))
    {
      if(!((i) %in% current_BordersRemoval_two) && !((i-1) %in% current_BordersRemoval_two) && !((i+1) %in% current_BordersRemoval_two)){
        current <- tads_intervals[length(tads_intervals)/2+i]:tads_intervals[i]
        tads_intervals<-rbind(tads_intervals, range(current))
      }
      if(i %in% current_BordersRemoval_two){
        current <- tads_intervals_one[length(tads_intervals_one)/2+j]:tads_intervals_one[j]
        tads_intervals<-rbind(tads_intervals, range(current))
        current <- tads_intervals_one[length(tads_intervals_one)/2+j+1]:tads_intervals_one[j+1]
        tads_intervals<-rbind(tads_intervals, range(current))
        j <- j+2
      }
    }
    if((tads_intervals[i+length(tads_intervals)/2] > tads_intervals_one[length(tads_intervals_one)])){
      current <- tads_intervals[length(tads_intervals)/2+i]:tads_intervals[i]
      tads_intervals<-rbind(tads_intervals, range(current))
    }
  }
  
  
  #----------------------- third layer ------------------------
  
  current_mysize <-NULL
  
  tads_intervals_three <- NULL
  # print(length(tads_intervals_three))
  borderSkippingProbs <- nestedProb
  size_two <- length(tads_intervals)/2-length(mysize)
  remove_size <- (round(borderSkippingProbs*(size_two))-1)
  every_space <- round(size_two/remove_size)
  #current_BordersRemoval_one <- array(1:remove_size)
  current_BordersRemoval_one <- rep(0, remove_size)
  for (i in 1:(remove_size-1)){
    top_space <- sample(1:every_space, 1)
    current_BordersRemoval_one[i] <- every_space*i + top_space
  }
  
  #current_BordersRemoval_one<-sort(sample(4:((size_two)-3), size=round(borderSkippingProbs*(size_two))))
  
  for(imerge in 1:(size_two-5)){
    if (imerge %in% current_BordersRemoval_one) {
      j <- sample(1:4,size = 1,prob = c(8,2,3,4))
      if((tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]-tads_intervals[length(mysize)+imerge]) > nestedSizeCap ){
        j <- j-1
      }
      if((tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]-tads_intervals[length(mysize)+imerge]) > nestedSizeCap ){
        j <- j-1
      }
      if((tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]-tads_intervals[length(mysize)+imerge]) < nestedSizeCap ){
        current <- tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+j]:tads_intervals[length(mysize)+imerge]
        if(length(tads_intervals_three) == 0 ){
          tads_intervals_three <- rbind(tads_intervals_three, range(current))
        }
        if(!(max(current) == tads_intervals_three[length(tads_intervals_three)])){
          tads_intervals_three <- rbind(tads_intervals_three, range(current))
        }
      }
    }
    
  }

  print(tads_intervals_three)
  
  j <- 1
  for(i in 1:(size_two)){
    if(j <= (length(tads_intervals_three)/2)){
      if(tads_intervals[i+length(mysize)] < tads_intervals_three[j]){
        current <- tads_intervals[length(mysize)+length(tads_intervals)/2+i]:tads_intervals[length(mysize)+i]
        tads_intervals <- rbind(tads_intervals, range(current))
      }
      if(tads_intervals[i+length(mysize)] == tads_intervals_three[j]){
        current <- tads_intervals_three[length(tads_intervals_three)/2+j]:tads_intervals_three[j]
        tads_intervals <- rbind(tads_intervals, range(current))
      }
      if(tads_intervals[i+length(mysize)+length(tads_intervals)/2] == tads_intervals_three[length(tads_intervals_three)/2+j]){
        j <- j+1
      }
    }
    if(tads_intervals[i+length(mysize)+length(tads_intervals)/2] > tads_intervals_three[length(tads_intervals_three)]){
      if(j == (length(tads_intervals_three)/2+1)){
        current <- tads_intervals[length(mysize)+length(tads_intervals)/2+i]:tads_intervals[length(mysize)+i]
        tads_intervals <- rbind(tads_intervals, range(current))
      }
      
    }
  }  
  #tads_intervals_copy <- tads_intervals
  #tads_intervals <-  TADs_dup_rm(tads_intervals)
  
  #-----------------------write tads————————————————————————
  #print(length(tads_intervals))
  #print(tads_intervals)
  pre_existing_tads<-length(tads)
  for (i in ((size_one+1):(length(tads_intervals)/2))){
    current <- tads_intervals[length(tads_intervals)/2+i] :tads_intervals[i]
    tads[[i]] <- t(combn(current, 2L))
    colnames(tads[[i]]) <- c("anchor", "target")
  }
  
  colnames(tads_intervals)<-c("startbin", "endbin")
  return(list(tads_list=tads, tads_intervals=tads_intervals))
}
  


spawnTADs_withNesting_dualOutput_6_3<- function(ntads, min.tad=10, max.tad=30, nestedLevels=3, nestedProb=0.3, nestedSizeCap=NULL, overlapProbs=NULL) 
{
  total_ref<-0
  if (!is.null(nestedSizeCap)) {
    if (nestedSizeCap<=max.tad) {
      stop("nestedSizeCap should be larger than max.tad")
    }
  }
  ntads <- as.integer(ntads)
  mysize <- as.integer(runif(ntads, min.tad, max.tad))
  tads <- list()
  tads_intervals<-NULL
  tads_intervals_one <-NULL
  mysize_ref<-mysize
  
  
  #----------------------first layer----------------
  last.launch <- 0L
  pre_existing_tads<-length(tads)
  for (i in 1:length(mysize_ref)) { 
    current <- last.launch + mysize_ref[i]:1   ## %%% catch max min of this to return the TADs start end positions
    tads_intervals<-rbind(tads_intervals, range(current))
    tads[[i+pre_existing_tads]] <- t(combn(current, 2L))
    colnames(tads[[i+pre_existing_tads]]) <- c("anchor", "target")
    last.launch <- last.launch + mysize_ref[i]
  }
  print("第一层")
  print(length(tads))
  
  
  #-----------------------------------second layer--------------------------
  size_one <- length(tads_intervals)/2
  borderSkippingProbs <- overlapProbs
  remove_size_two <- round(borderSkippingProbs*(size_one))
  every_space_two <- floor(size_one/remove_size_two)
  current_BordersRemoval_two <- rep(0, remove_size_two)
  
  for (i in 1:(remove_size_two)){
    top_space <- sample(2:(every_space_two-1))
    current_BordersRemoval_two[i] <- every_space_two*(i-1) + top_space
  }
  #print(current_BordersRemoval_two)
  for(imerge in 1:(size_one)){
    if (imerge %in% current_BordersRemoval_two) {
      random_overlap_number_left <- as.integer(runif(1,min=1,max=1.1))
      random_overlap_number_right <- as.integer(runif(1,min=1,max=1.1))
      current <-  (tads_intervals[length(mysize_ref)+imerge]):tads_intervals[imerge - random_overlap_number_left]
      tads_intervals_one<-rbind(tads_intervals_one, range(current))
      current <- (tads_intervals[length(mysize_ref)+imerge+random_overlap_number_right]):tads_intervals[imerge]
      tads_intervals_one<-rbind(tads_intervals_one, range(current))
    }
  }
  print(tads_intervals_one)
  j <- 1
  
  for(i in 1:(size_one)){
    if((tads_intervals[i+length(tads_intervals)/2] <= tads_intervals_one[length(tads_intervals_one)]))
    {
      if(!((i) %in% current_BordersRemoval_two) && !((i-1) %in% current_BordersRemoval_two) && !((i+1) %in% current_BordersRemoval_two)){
        current <- tads_intervals[length(tads_intervals)/2+i]:tads_intervals[i]
        tads_intervals<-rbind(tads_intervals, range(current))
      }
      if(i %in% current_BordersRemoval_two){
        current <- tads_intervals_one[length(tads_intervals_one)/2+j]:tads_intervals_one[j]
        tads_intervals<-rbind(tads_intervals, range(current))
        current <- tads_intervals_one[length(tads_intervals_one)/2+j+1]:tads_intervals_one[j+1]
        tads_intervals<-rbind(tads_intervals, range(current))
        j <- j+2
      }
    }
    if((tads_intervals[i+length(tads_intervals)/2] > tads_intervals_one[length(tads_intervals_one)])){
      current <- tads_intervals[length(tads_intervals)/2+i]:tads_intervals[i]
      tads_intervals<-rbind(tads_intervals, range(current))
    }
  }
  
  
  #-----------------------third layer------------------------
  
  current_mysize <-NULL
  
  tads_intervals_three <- NULL
  borderSkippingProbs <- nestedProb
  size_two <- length(tads_intervals)/2-length(mysize)
  remove_size <- round(borderSkippingProbs*(size_two))-1
  every_space <- round(size_two/remove_size)
  current_BordersRemoval_one <- array(1:remove_size)
  
  for (i in 1:(remove_size-1)){
    top_space <- sample(1:every_space)
    current_BordersRemoval_one[i] <- every_space*i + top_space
  }
  
  #current_BordersRemoval_one<-sort(sample(4:((size_two)-3), size=round(borderSkippingProbs*(size_two))))
  
  
  for(imerge in 1:size_two){
    if (imerge %in% current_BordersRemoval_one) {
      j <- sample(1:4,size = 1)
      if((tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]-tads_intervals[length(mysize)+imerge]) > nestedSizeCap){
        j <- j-1
      }
      if((tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]-tads_intervals[length(mysize)+imerge]) > nestedSizeCap ){
        j <- j-1
      }
      if((tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]-tads_intervals[length(mysize)+imerge]) < nestedSizeCap ){
        current <- tads_intervals[length(mysize)+length(tads_intervals)/2+imerge+1+j]:tads_intervals[length(mysize)+imerge]
        tads_intervals_three <- rbind(tads_intervals_three, range(current))
      }
    }
    
  }
  #print(tads_intervals_three)
  
  j <- 1
  for(i in 1:size_two){
    if(j <= (length(tads_intervals_three)/2)){
      if(tads_intervals[i+length(mysize)] < tads_intervals_three[j]){
        current <- tads_intervals[length(mysize)+length(tads_intervals)/2+i]:tads_intervals[length(mysize)+i]
        tads_intervals <- rbind(tads_intervals, range(current))
      }
      if(tads_intervals[i+length(mysize)] == tads_intervals_three[j]){
        current <- tads_intervals_three[length(tads_intervals_three)/2+j]:tads_intervals_three[j]
        tads_intervals <- rbind(tads_intervals, range(current))
      }
      if(tads_intervals[i+length(mysize)+length(tads_intervals)/2] == tads_intervals_three[length(tads_intervals_three)/2+j]){
        j <- j+1
      }
    }
    if(tads_intervals[i+length(mysize)+length(tads_intervals)/2] > tads_intervals_three[length(tads_intervals_three)]){
      if(j == (length(tads_intervals_three)/2+1)){
        current <- tads_intervals[length(mysize)+length(tads_intervals)/2+i]:tads_intervals[length(mysize)+i]
        tads_intervals <- rbind(tads_intervals, range(current))
      }
      
    }
  }  
  
  
  #-----------------------writetads————————————————————————
  #print(length(tads_intervals))
  #print(tads_intervals)
  pre_existing_tads<-length(tads)
  for (i in ((size_one+1):(length(tads_intervals)/2))){
    current <- tads_intervals[length(tads_intervals)/2+i] :tads_intervals[i]
    tads[[i]] <- t(combn(current, 2L))
    colnames(tads[[i]]) <- c("anchor", "target")
  }
  
  colnames(tads_intervals)<-c("startbin", "endbin")
  return(list(tads_list=tads, tads_intervals=tads_intervals))
}

#---------------------------------------------------
#---------------------------------------------------


#####################################################
#####################################################

spawnTADs_withNesting_dualOutput_origin <- function(ntads, min.tad=10, max.tad=30, nestedLevels=3, nestedProb=0.3, nestedSizeCap=NULL) 
{
  
  #check parameters
  if (!is.null(nestedSizeCap)) {
    if (nestedSizeCap<=max.tad) {
      stop("nestedSizeCap should be larger than max.tad")
    }
  }
  ####as.integer()
  ntads <- as.integer(ntads)
  ###runif()
  mysize <- as.integer(runif(ntads, min.tad, max.tad))
  tads <- list()
  tads_intervals<-NULL
  
  # copy object for subsequent manipulations
  mysize_ref<-mysize
  
  for (nestingLevel in 1:nestedLevels) {
    if (nestingLevel==1) {
      borderSkippingProbs=0
    } else {
      borderSkippingProbs=nestedProb
    }
    
    
    # sampling TADs toc 
    
    
    #  be removed with probability nestedProb  
    # the selected TADs will be merged to the following one, as such we sample up to (mysize length -1)
    current_BordersRemoval<-sort(sample(1:(length(mysize_ref)-1), size=round(borderSkippingProbs*(length(mysize_ref)-1))))
    
    current_mysize<-NULL
    for (imerge in 1:length(mysize_ref)) {
      if (imerge %in% current_BordersRemoval) {
        # cap the maximum size of nested TADs
        # cap the maximum size of nested TADs 
        if ( !is.null(nestedSizeCap)) {
          if ((mysize_ref[imerge]+mysize_ref[(imerge+1)]) > nestedSizeCap) {
            # skip the merging it the merged TAD would be larger than the cap size
            current_mysize<-c(current_mysize, mysize_ref[imerge])
          } else {
            mysize_ref[(imerge+1)]<-(mysize_ref[imerge]+mysize_ref[(imerge+1)])
          }
        } else {
          mysize_ref[(imerge+1)]<-(mysize_ref[imerge]+mysize_ref[(imerge+1)])
        }
      } else {
        current_mysize<-c(current_mysize, mysize_ref[imerge])
      }
    }
    mysize_ref<-current_mysize
    
    
    # > current_BordersRemoval
    #  [1]  49  44  51  33 148 134  13 111 169  99  35  39  72 108  65  90   2  78 121
    # [20]  31   4 125 124 150  16 145 113  26  54  95  83 105  58 117 137 163 155  56
    # [39]  25  38  85 154  43 106  73 118 153  70  74  14 133
    # > sort(current_BordersRemoval)
    #  [1]   2   4  13  14  16  25  26  31  33  35  38  39  43  44  49  51  54  56  58
    # [20]  65  70  72  73  74  78  83  85  90  95  99 105 106 108 111 113 117 118 121
    # [39] 124 125 133 134 137 145 148 150 153 154 155 163 169
    # 
    
    last.launch <- 0L
    pre_existing_tads<-length(tads)
    for (i in 1:length(mysize_ref)) { 
      current <- last.launch + mysize_ref[i]:1   ## %%% catch max min of this to return the TADs start end positions
      tads_intervals<-rbind(tads_intervals, range(current))
      tads[[i+pre_existing_tads]] <- t(combn(current, 2L))
      colnames(tads[[i+pre_existing_tads]]) <- c("anchor", "target")
      last.launch <- last.launch + mysize_ref[i]
    }
  }	
  
  colnames(tads_intervals)<-c("startbin", "endbin")
  return(list(tads_list=tads, tads_intervals=tads_intervals))
}
#------------------------------------------------------------------------
#------------------------------------------------------------------------

spawnTADs_dualOutput <- function(ntads, min.tad=10, max.tad=30) 
{
  ntads <- as.integer(ntads)
  mysize <- as.integer(runif(ntads, min.tad, max.tad))
  tads <- list()
  last.launch <- 0L
  tads_intervals<-NULL
  for (i in 1:ntads) {
    current <- last.launch + mysize[i]:1   ## %%% catch max min of this to return the TADs start end positions
    tads_intervals<-rbind(tads_intervals, range(current))
    tads[[i]] <- t(combn(current, 2L))
    colnames(tads[[i]]) <- c("anchor", "target")
    last.launch <- last.launch + mysize[i]
  }
  colnames(tads_intervals)<-c("startbin", "endbin")
  return(list(tads_list=tads, tads_intervals=tads_intervals))
}


### This is the same as original Aaron Lun function
getMeans <- function(mat, base=10, decay=0.5, prior=1, rate=1) 
  # This makes mean values for the requested matrices that are coming in,
  # based on the log-log relationship observed between mean and distance.
  
{
  distance <- mat[,1]-mat[,2] + prior
  true.mean <- base*distance^(-decay)
  if (is.na(rate)) { return(true.mean) } else {
    return(rgamma(nrow(mat), shape=true.mean*rate, rate=rate)) } ## rate is set to NA below thus it is returning the treu mean
}


### This is the same as original Aaron Lun function
spawnDiagon <- function(len) 
  # Adds coordinates for the diagonals. 
 
{
  len <- as.integer(len)
  cbind(anchor=1:len, target=1:len)#cbind
}

### This is the same as original Aaron Lun function
# but with more comments added
spawnLigations2 <- function(nspec, len, diag=FALSE, prior=1, decay=0.5) 
  # This adds ligation events where the weighting of the selected ligation events is
  # skewed towards the diagonal. This avoids putting too much work into the rest of the
  # contact matrix that is unlikely to be interesting. 
  
{ 
  len <- as.integer(len) ## %%% number of rows (diagonals) in the matrix
  all.dists <- (len-1):ifelse(diag, 0, 1) ## %%%  define whether to include or not the diagonal x=y (no in the script below)
  n.entries <- len - all.dists ## %%%  define how many data points are available at each considered diagonal (distance)
  chosen.dists <- sample(all.dists, nspec, p=n.entries*(all.dists+prior)^(-decay), replace=TRUE)  ## %%% sample the data points with replacement and probability proportional to how many data points are available at each diagonal and decreasing with power law (x-y+p)^c
  
  chosen <- split(1:nspec, chosen.dists) ## %%% group the selected data points by diagonals and format as list
  anchor <- target <- integer(nspec) ## %%% define empty vector
  for (dist in names(chosen)) { ## %%% for each diagonal
    relevants <- chosen[[dist]] ## %%% how many data points needs to be samples (result from sample above)
    allowable.anchors <- (as.integer(dist)+1L):len ## %%%  how many possible data points in that diagonal?
    if(length(allowable.anchors)>1){
      cur.anchors <- sample(allowable.anchors, length(relevants), replace=TRUE) ## %%% sample with replacement among possible data points of that diagonals, sample length(relevants) elements which is the result of the previous sample function call
    }else{
      cur.anchors <- rep(allowable.anchors, length(relevants))
    }
    cur.targets <- cur.anchors - as.integer(dist)  ## %%% define the y coordinates by subtracting the diagonals
    anchor[relevants] <- cur.anchors ## %%%  save value sin vector
    target[relevants] <- cur.targets ## %%%  save value sin vector 
  }
  
  o <- order(anchor, target) ## %%% sort and return
  cbind(anchor=anchor[o], target=target[o])
}


### This is the same as original Aaron Lun function
spawnLigations <- function(nspec, max.range, len, diag=TRUE) 
  # This adds ligation events within a certain range of each other. This avoids
  # simulating the far-out interactions which are mostly irrelevant.
{ 
  max.range <- as.integer(max.range)
  len <- as.integer(len)
  if (max.range > len) { stop("max.range must be less than the total length of the chromosome") }
  
  per.height <- pmin(max.range + 1L, len-1:len+1L)
  if (!diag) { per.height <- per.height[-len] - 1L }
  total.possibles <- sum(per.height)
  if (total.possibles <= nspec) { 
    warning("reporting everything as requested") 
    chosen <- 1:total.possibles
  } else {
    chosen <- sample(total.possibles, nspec) 
  }
  
  ratchet <- c(1L, cumsum(per.height)+1L)
  target <- findInterval(chosen, ratchet) 
  anchor <- chosen - ratchet[target] + target
  if (!diag) { anchor <- anchor + 1L }
  cbind(anchor=anchor, target=target)
}

### This is the same as original Aaron Lun function
dummyFragments <- function(len, name="chrA") {
  names(len) <- name
  GRanges(name, IRanges(1:len, 1:len), seqlengths=len)
}


### This is the same as original Aaron Lun function
aggregateAll <- function(...) 
  # This takes a bunch of (coordinate matrix, count vector) lists, and
  # merges them to compute the aggregate mean.
{
  everything <- list(...)
  all.anchors <- all.targets <- all.counts <- list()
  counter <- 1L
  for (x in everything) { 
    all.anchors[[counter]] <- x[[1]][,1]
    all.targets[[counter]] <- x[[1]][,2]
    all.counts[[counter]] <- x[[2]]
    counter <- counter + 1L
  }
  
  # Reordering.
  all.anchors <- unlist(all.anchors)
  all.targets <- unlist(all.targets)
  all.counts <- unlist(all.counts)
  o <- order(all.anchors, all.targets)
  all.anchors <- all.anchors[o]
  all.targets <- all.targets[o]
  all.counts <- all.counts[o]
  
  # Finding the sum.
  keep <- which(diff(all.anchors)!=0L | diff(all.targets)!=0L)
  first <- c(1L, keep + 1L)
  last <- c(keep, length(o))
  cum.count <- cumsum(all.counts)
  collected <- cum.count[last] - cum.count[first] + all.counts[first]
  
  # Returning the aggregated output.
  list(cbind(anchor=all.anchors[first], target=all.targets[first]), collected)
}



#####################################################
#####
##### FROM NOW ON THE CODE IS MORE DISTAND FROM THE ORIGINAL CODE BY Aaron Lun
##### still the structure of the procedure (in term of steps and sequence of calls 
##### to the accessory functions defined above) is conserved as the underlying rationale
##### of the procedure is used here as well
#####



#####################################################
#####################################################
#####
##### ADDING VANILLA COVERAGE FUNCTION FOR NORMALIZATION
#####


# x is the squared symmetric contact matrix
VanillaCoverage<-function(OBS) {
  
  # just make sure it is a matrix
  if (!is.matrix(OBS)) {
    stop("OBS must be a matrix")
  }
  
  # grand total
  T<-sum(OBS)
  
  # rows and column sums
  Cs<-colSums(OBS)
  Rs<-rowSums(OBS)
  
  # expected counts given the total by column/row and the grandtotal of read counts
  EXP<-((Cs) %*% t(Rs))/T
  
  # returning the normalized matrix
  return(OBS/EXP)
  
}



#####################################################
#####################################################


#### SAMPLE PARAMETERS FOR DEBUGGING
# simChromName<-"ChrS"
# BINSIZE<-40000
# max.tad <- 20  
# min.tad<-3
# TARGET.chrsize<-180.92*1e6 # same as chr5 in hg19
# base.tad.signal<-80
# interaction.decay.rate<-0.8
# percent_of_sampling_points_nonSpecificsignal<-0.04
# n.interactions=NULL
# base.diagonal.signal=NULL
# max.range=NULL
# base.interaction.signal=NULL
# dispersion <- 0.01 
# PNGplot=FALSE
# saveVanilla=FALSE
# saveDiffHiC=FALSE
# saveHOMER=FALSE
# addIntConstant=FALSE
# base.nonspecific.signal=NULL
# nestedTADs=TRUE
# nestedLevels=3
# nestedProb=0.25
# nestedSizeCap=75
#--------------------------------------------------------------------------------
#################################################################################
#--------------------------------------------------------------------------------

runHiCSimulation<-function(BINSIZE, simChromName="ChrS", max.tad=50, min.tad=3,
                           TARGET.chrsize=180.92*1e6, base.tad.signal=80, interaction.decay.rate=0.8,
                           percent_of_sampling_points_nonSpecificsignal=0.04, n.interactions=NULL,
                           base.diagonal.signal=NULL, max.range=NULL, base.interaction.signal=NULL, base.nonspecific.signal=NULL,
                           dispersion=0.01, PNGplot=FALSE, saveVanilla=TRUE, saveDiffHiC=FALSE, saveHOMER=FALSE, addIntConstant=FALSE,
                           nestedTADs=FALSE, nestedLevels=3, nestedProb=0.3, nestedSizeCap=NULL,overlapProbs=NULL,output=NULL, samplenumber, nreplicate) {
  
  
  
  require(diffHic)
  require(edgeR)
  
  #debugSource('drawHeatmap.R')
  
  
  TIME<-(Sys.time()) # can be used to set seed
  # in this version it is just used to define the output filename
  
  ## make this linked to the target chromosome size
  n.tads <- round((TARGET.chrsize/BINSIZE)/mean(min.tad:max.tad))
  
  ## round()
  #  signal in diagonal
  if (is.null(base.diagonal.signal)) {
    # proportional to original Aron Lun choices (28*100/80)
    base.diagonal.signal<-(base.tad.signal*100/80)
  }
  
  estimated_binsperchromsome<-round((TARGET.chrsize/BINSIZE))
  estimated_matrix_size_upperTriangle<-(((estimated_binsperchromsome)^2-(estimated_binsperchromsome))/2)
  
  ## %%% changed to make this proportional to changing parameters
  SampledPointNonSpecificSignal<-round(percent_of_sampling_points_nonSpecificsignal * estimated_matrix_size_upperTriangle)
  
  if (is.null(base.nonspecific.signal)) {
    ## %%% changed to make this proportional to changing parameters
    base.nonspecific.signal<-round((5/80)*base.tad.signal) 
  }
  
  # ## original was 100 --> over matrix size 100/4840716 =  2.06581e-05 (i.e. 2 every 100k)
  if (is.null(n.interactions)) {
    ## %%% changed to make this proportional to changing parameters
    n.interactions <-  round(2e-5*estimated_matrix_size_upperTriangle)  
  }
  
  
  
  if (is.null(max.range)) {
    # max distance from diagonal for the interactions (no more than 10 bins beyond TAD max size as there are few interactions very far from diagonal in real data)
    ## %%% changed to make this proportional to changing parameters 
    ### %%% not going beyond diagonal #40... we may consider also changing this parameter
    max.range <- max.tad + round(max.tad/3) 
  }
  
  if (is.null(base.interaction.signal)) {
    # signal at interaction
    ## %%% changed to make this proportional to changing parameters
    base.interaction.signal<-base.tad.signal*2
  }
  
  
  print("Simulate data...")
  
  ###generate tads
  if (nestedTADs) {
    if(overlapProbs==0){
       TadsSimulated <- spawnTADs_withNesting_dualOutput_origin(ntads=n.tads, min.tad=min.tad, max.tad=max.tad, nestedLevels=nestedLevels, nestedProb=nestedProb, nestedSizeCap=nestedSizeCap)
    tad.co<-TadsSimulated$tads_list
    tad.co_intervals<-TadsSimulated$tads_intervals
    }else{
      TadsSimulated <- spawnTADs_withNesting_dualOutput(ntads=n.tads, min.tad=min.tad, max.tad=max.tad, nestedLevels=nestedLevels, nestedProb=nestedProb, nestedSizeCap=nestedSizeCap,overlapProbs=overlapProbs) 
      tad.co<-TadsSimulated$tads_list
      tad.co_intervals<-TadsSimulated$tads_intervals
    }
  } else {
    TadsSimulated <- spawnTADs_dualOutput(ntads=n.tads, min.tad=min.tad, max.tad=max.tad) 
    tad.co<-TadsSimulated$tads_list
    tad.co_intervals<-TadsSimulated$tads_intervals
  }
  
  tad.co_intervals <- tad.co_intervals[!duplicated(tad.co_intervals), ]
  tad.co_intervals[, 2] <- tad.co_intervals[, 2]+1
 
  tad.co <- do.call(rbind, tad.co)
  tad.mu <- getMeans(tad.co, base=base.tad.signal, decay=interaction.decay.rate, rate=NA)
  

  #generate noise, nb distribution, interaction, self
  for (r in 1:nreplicate) {
    
    nlen <- max(tad.co)
    self.co <- spawnDiagon(nlen) ## 
    self.mu <- rep(base.diagonal.signal, nlen)   ### diagonal signal
    # Adding some random noise. 
    nnspec <- SampledPointNonSpecificSignal
    
    if(nnspec == 0){
      random.co <- NULL
      random.mu <- NULL
    }else{
      random.co <- spawnLigations2(nnspec, nlen, decay=interaction.decay.rate) ## %%% not defined
      random.mu <- rep(base.nonspecific.signal, nnspec)
      #randomNoiseTotals<-tapply(random.mu, INDEX=paste("anch", random.co[,1], "tar", random.co[,2], "bin", sep="_"), FUN=sum)
    }
    if(n.interactions == 0){
      inter.co <- NULL
      inter.mu <- NULL
    }else{
      inter.co <- spawnLigations(n.interactions, max.range, nlen, diag=FALSE)  ## %%% not defined
      inter.mu <- getMeans(inter.co, base=base.interaction.signal, decay=interaction.decay.rate, rate=NA) ## %%% not defined
    }
    
    if (addIntConstant) {
      inter.mu <- (inter.mu + base.nonspecific.signal)
    }
    blocks <- dummyFragments(nlen, name=simChromName) ## %%% not defined
    xparam <- pairParam(blocks) 
    # 		for (grouping in 1:2) {
    cur.mus <- inter.mu
    # 		if (grouping==1L) {
    cur.inters <- inter.co # keeping everything as we simulate just one group
    collated <- aggregateAll( list(tad.co, tad.mu), list(self.co, self.mu), list(random.co, random.mu), list(cur.inters, cur.mus) )
    
    counts <- rnbinom(length(collated[[2]]), mu=collated[[2]], size=1/dispersion)
    
    
    
    COUNT.MATRIX<-matrix(0, nrow=nlen, ncol=nlen)
    COUNT.MATRIX[(collated[[1]])]<-counts
    ## we want a symmetric matrix
    COUNT.MATRIX[((collated[[1]])[,c(2,1)])]<-counts
    Bin_labels<-paste(simChromName, 1:nlen, sep="_")
    colnames(COUNT.MATRIX)<-Bin_labels
    OUT.COUNT.MATRIX<-data.frame(COUNT.MATRIX, stringsAsFactors=FALSE)
    OUTFILES_suffix<-paste(simChromName, "Bin", BINSIZE, "TIME", format(TIME, "%Y-%m-%d_%Hh%Mm%Ss"), sep=".")
    
    ########################write tads#######################
    print("Save TADs...")
    if(nreplicate == 1){
      outfilename_tads<-paste(output,"/ChrS_TADs", "_noise", percent_of_sampling_points_nonSpecificsignal,"_POP",overlapProbs, ".txt", sep="")
    }else if(nreplicate >1){
      outfilename_tads<-paste(output,"/ChrS_TADs", "_noise", percent_of_sampling_points_nonSpecificsignal, "_POP",overlapProbs, "_rep", r, ".txt", sep="")
    }else{
      stop("replicate number error")
    }
    write.table(tad.co_intervals, file=outfilename_tads, sep="\t", row.names=FALSE, quote=FALSE)
    
    
      ########################write matrix#######################
    print("Save counts...")
    if(nreplicate == 1){
      chr_filename<-paste(output,"/ChrS_MAT", "_noise", percent_of_sampling_points_nonSpecificsignal,"_POP",overlapProbs, ".txt", sep="")
    }else if(nreplicate > 1){
      chr_filename<-paste(output,"/ChrS_MAT", "_noise", percent_of_sampling_points_nonSpecificsignal, "_POP",overlapProbs, "_rep", r, ".txt", sep="")
    }else{
      stop("replicate number error")
    }
    write.table(OUT.COUNT.MATRIX, file=chr_filename, sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)
    
    
    ########################draw heatmap#######################
    print("Drawing heatmap")
    binStart=1
    binEnd=nrow(COUNT.MATRIX)
    if(nreplicate == 1){
      chr_figurename<-paste(output,"/ChrS_MAT_Figure", "_noise", percent_of_sampling_points_nonSpecificsignal,"_POP",overlapProbs,"_range",binStart,"_", binEnd, ".tiff", sep="")
      chr_figurename2 <- paste(output,"/ChrS_MAT_Figure_Mark", "_noise", percent_of_sampling_points_nonSpecificsignal,"_POP",overlapProbs,"_range",binStart,"_", binEnd, ".tiff", sep="")
    }else if(nreplicate > 1){
      chr_figurename<-paste(output,"/ChrS_MAT_Figure", "_noise", percent_of_sampling_points_nonSpecificsignal,"_POP",overlapProbs,"_range",binStart,"_", binEnd, "_rep", r, ".tiff", sep="")
      
      chr_figurename2 <- paste(output,"/ChrS_MAT_Figure_Mark", "_noise", percent_of_sampling_points_nonSpecificsignal,"_POP",overlapProbs,"_range",binStart,"_", binEnd, "_rep", r, ".tiff", sep="")
    }else{
      stop("replicate number error")
    }
    
    COUNT.PLOT <- COUNT.MATRIX[binStart:binEnd,binStart:binEnd]
    COUNT.PLOT <- log10(1+COUNT.PLOT)
    TADsCoor.PLOT <<- genRectFromTADCoor(tad.co_intervals, binStart, binEnd)
    
    library(gplots)
    tiff(chr_figurename, width = 18, height = 18, units = "cm", pointsize=1, res=350)
    heatmap.2(COUNT.PLOT, Rowv=FALSE,col=colorpanel(128,"lightyellow","red"),symm=TRUE,dendrogram="none",trace = "none",labRow = FALSE,labCol = FALSE,density.info="none",key=FALSE,lmat=rbind(c(1,4),c(3,2)),lhei=c(7,1),lwid=c(7,1), add.expr=for (i in 1:nrow(TADsCoor.PLOT)){rect(xleft = TADsCoor.PLOT[i,1],ybottom = TADsCoor.PLOT[i,3],xright = TADsCoor.PLOT[i,2],ytop = TADsCoor.PLOT[i,4],angle = 0, border = "blue", col = NULL,lwd = 1, lty = 1)})
    dev.off()
    
    tiff(chr_figurename2, width = 18, height = 18, units = "cm", pointsize=1, res=350)
    heatmap.2(COUNT.PLOT, Rowv=FALSE,col=colorpanel(128,"lightyellow","red"),symm=TRUE,dendrogram="none",trace = "none",labRow = FALSE,labCol = FALSE,density.info="none",key=FALSE,lmat=rbind(c(1,4),c(3,2)),lhei=c(7,1),lwid=c(7,1))
    dev.off()
  }
  return(OUTFILES_suffix)
}


KtBaselineTADs<-28
nestedLevels<-3
nsample <- 1
nrep <- 5


###create corresponding folders before running the scriptes 
for (overlap in c(0.15, 0)) {
  for (noise in c(0.20, 0.04, 0.08, 0.12, 0.16)) {
    outdir <- paste0("/media/disk1/liuerhu/hicda_liuerhu/simdata/overlap/output_final/overlap",overlap,"/noise",noise*100)
    for (i in 1:nsample) {
      check<-runHiCSimulation(BINSIZE=25000, simChromName="ChrS", max.tad=20, min.tad=5,TARGET.chrsize=110*1e5, base.tad.signal=round(KtBaselineTADs/nestedLevels), interaction.decay.rate=0.45, percent_of_sampling_points_nonSpecificsignal=noise, n.interactions=10,base.diagonal.signal=(KtBaselineTADs*100/80), max.range=(50+round(50*1/3)), base.interaction.signal=(KtBaselineTADs*2),base.nonspecific.signal=round((9/80)*KtBaselineTADs), dispersion=0.001, PNGplot=TRUE, saveVanilla=TRUE,saveDiffHiC=FALSE, saveHOMER=FALSE, addIntConstant=TRUE,nestedTADs=TRUE,nestedLevels=nestedLevels, nestedProb=0.25, nestedSizeCap=75,overlapProbs = overlap, output = outdir, samplenumber = i, nreplicate = nrep)
      print(paste("completed", check))
    }
  }
}
