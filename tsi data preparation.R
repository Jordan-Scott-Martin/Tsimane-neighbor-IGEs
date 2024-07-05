#prep workspace############################################
# setwd("...")
library(reshape2)
library(MCMCglmm)
library(kinship2)
library(dplyr)
library(Matrix)

#load data
ped = read.csv("pedigrees.csv")
birth = read.csv("births.csv")
comm = read.csv("communities.csv")
mom = read.csv("mothers.csv")
neigh = read.csv("neighbor_dyads.csv")

#fix pedigree errors############################################
pedg = ped[,c("pid","mother_pid","father_pid","year_of_birth")]
colnames(pedg) = c("id","dam","sire","time_born")

#remove uncertain pedigree entries based on conflicting pids
error_id = as.vector(na.omit(intersect(pedg$dam, pedg$sire))) #same pid for mothers and fathers 
pedg = pedg[!(pedg$id %in% error_id | pedg$dam %in% error_id | pedg$sire %in% error_id  ), ]

#remove rows without information on parents
pedg = na.omit(pedg)

orderPed <-function(ped, time_born=NULL){
  
  if(length(time_born)!=0){
    reorder<-order(time_born)
  }else{
    reorder<-order(kindepth(ped[,1],ped[,2],ped[,3]), decreasing=FALSE)
  }

  ped[reorder,]
}

#order
pedg = orderPed(pedg)

#add phantom IDs for missing parents of mothers and fathers
phnt.id1 = unique(pedg[(pedg[, "dam"] %in% pedg[, "id"])==FALSE,"dam"])
phnt1 = data.frame(id = phnt.id1, dam = NA, sire = NA, time_born = NA) #NA birth dates to preserve order
phnt.id2 = unique(pedg[(pedg[, "sire"] %in% pedg[, "id"])==FALSE,"sire"])
phnt2 = data.frame(id = phnt.id2, dam = NA, sire = NA, time_born = NA)
pedg2 = rbind(phnt1, phnt2, pedg)

#prune to only individuals with phenotypic data
mom.comm = na.omit(unique(birth[,c("mother_pid","father_pid","community_id")]))
keep_id = unique(c(mom.comm$mother_pid,mom.comm$father_pid))
pedg3 = prunePed(pedg2[,1:3], keep = keep_id, make.base = TRUE)

#create initial pedigree############################################
invA = inverseA(pedg3[,c("id","dam","sire")])
A = solve(invA$Ainv)
A = cov2cor(A)
rownames(A) = rownames(invA$Ainv)
colnames(A) = rownames(invA$Ainv)

#prepare phenotypic data for analysis##############################

#remove births without community IDs
b = birth[,c("pid","community_id","year_of_birth","mother_pid","father_pid")]
b = b[order(b$mother_pid),]
b = b[!is.na(b$community_id),]

#sampling window for focal fertility/community
years = aggregate(year_of_birth ~ mother_pid * community_id, data = b, 
              FUN = function(x) max(x) - min(x) )
colnames(years)[3] = "years"
b = left_join(b, years)

#year range
year_min = aggregate(year_of_birth ~ mother_pid * community_id, data = b, 
              FUN = function(x) min(x)  )
colnames(year_min)[3] = "year_min"
b = left_join(b, year_min)

year_max = aggregate(year_of_birth ~ mother_pid * community_id, data = b, 
              FUN = function(x) max(x)  )
colnames(year_max)[3] = "year_max"
b = left_join(b, year_max)

#create fertility/community outcome variable
birth.tbl = melt(table((b[,c("mother_pid","community_id")])), value.name = "fert")
birth.tbl = birth.tbl[birth.tbl$fert>0,]
birth.df = left_join(birth.tbl, 
             unique(b[,c("mother_pid","community_id","years","year_min","year_max")]))

#add community location
birth.df$lat = comm[match(birth.df$community_id,comm$community_id),"latitude"]
birth.df$long = comm[match(birth.df$community_id,comm$community_id),"longitude"]

#add mother age (relative to current year)
birth.df$mom_age = mom[match(birth.df$mother_pid, mom$pid),"year_of_birth"]
birth.df$mom_age = birth.df$year_min - birth.df$mom_age + 3 
#3 = + year prior, current year, subsequent year

#create neighbor matrix for Stan
neighbs = list()
fids = unique(unlist(neigh[,1:2]))
fids = na.omit(years$mother_pid[match(fids, years$mother_pid)])
for(i in 1:length(fids)){
  fid = fids[i]
  neighbh = unique(unlist(neigh[neigh$pid1==fid | neigh$pid2==fid, 1:2]))
  neighbs[[i]] = neighbh[neighbh!=fid]
}
neigh_id = matrix(NA, nrow = length(fids), ncol = max(mapply(neighbs, FUN = length)))
for(i in 1:nrow(neigh_id)){
  focal = fids[i]
  partners = unlist(neighbs[i])
  neigh_id[i,1:length(partners)] = partners
}

#only neighbors with fertility data
neigh.id = neigh_id[match(birth.df$mother_pid, fids),]

#duplicate birth dataframe for editing
data.df = birth.df

#remove mothers with <2 births/community
#unable to distinguish fertility from sampling window
remove.id = which(data.df$fert==1)
data.df = data.df[-remove.id,]
neigh.id = neigh.id[-remove.id,]

#remove focal mothers without data in pedigree
remove.ped = which(data.df$mother_pid %in% 
                   setdiff(unique(data.df$mother_pid), dimnames(A)[[1]]))
data.df = data.df[-remove.ped,]
neigh.id = neigh.id[-remove.ped,]

#remove communities without GPS
remove.com = which(is.na(data.df$lat))
data.df = data.df[-remove.com,]
neigh.id = neigh.id[-remove.com,]

#remove mothers without birth years
#unknown age-adjusted fertility
remove.id2 = which(is.na(data.df$mom_age))
data.df = data.df[-remove.id2,]
neigh.id = neigh.id[-remove.id2,]

#remove women with impossible birthdays
too.old = na.omit(birth.df[birth.df$mom_age>55 & birth.df$fert>1,c("mother_pid","mom_age")])
too.young = na.omit(birth.df[birth.df$mom_age<15 & birth.df$fert>1,c("mother_pid","mom_age")])
rma = which(data.df$mother_pid %in% c(too.old$mother_pid, too.young$mother_pid))
data.df = data.df[-rma,]
neigh.id = neigh.id[-rma,]

#add maternal ID (mother of the focal woman)
data.df$maternal = pedg3[match(data.df$mother_pid,pedg3$id),"dam"]
data.df[is.na(data.df$maternal),"maternal"] = "unknown" #offset cluster

#add birth father IDs (father of children, generally the spouse of focal woman)
fpid = list()
for(i in 1:nrow(data.df)){
    smtch = as.numeric(data.df[i, c("mother_pid", "community_id")])
    fathers = unique(as.numeric(na.omit(b[b$mother_pid == smtch[1] & b$community_id == smtch[2], "father_pid"])))
    if(length(fathers)==0){fpid[[i]] = NA}
    else(fpid[[i]] = fathers)
  }
father.id = matrix(NA, nrow = nrow(data.df), ncol = max(mapply(fpid, FUN = length)))
for(i in 1:nrow(father.id)){
  fathers = unlist(fpid[i])
  father.id[i,1:length(fathers)] = fathers
  }

#remove neighbors without pedigree information
rmn = setdiff(neigh.id, rownames(A))[-1] #remove NA
neigh.id[which(neigh.id %in% rmn)] = NA

#adjust IDs for indexing in Stan
key.id = unique(na.omit(c(data.df$mother_pid)))
key.fth = unique(na.omit(as.vector(father.id)))
key.com = unique(data.df$community_id)
key.mom = unique(data.df$maternal)
key.year = seq(min(data.df$year_min), max(data.df$year_max))
new.id = seq(1:length(key.id))
new.fth = seq(1:length(key.fth))
new.com = seq(1:length(key.com))
new.mom = seq(1:length(key.mom))
new.year = seq(1:length(key.year))

data.df$mom_id = new.id[match(data.df$mother_pid, key.id)]
data.df$com_id = new.com[match(data.df$community_id, key.com)]
data.df$maternal_id = new.mom[match(data.df$maternal, key.mom)]
data.df$yearf_id = new.year[match(data.df$year_min, key.year)]
data.df$yearl_id = new.year[match(data.df$year_max, key.year)]

neigh.id = matrix(new.id[match(neigh.id, key.id)], nrow = nrow(neigh.id), ncol = ncol(neigh.id))
neigh.id = t(apply(neigh.id, 1, function(x) `length<-`(na.omit(x), length(x))))
neigh.n = apply(neigh.id, 1, FUN = function(x) sum(!is.na(x)) )
neigh.id[is.na(neigh.id)] = 0

father.id = matrix(new.id[match(father.id, key.fth)], nrow = nrow(father.id), ncol = ncol(father.id))
father.id = t(apply(father.id, 1, function(x) `length<-`(na.omit(x), length(x))))
father.n = apply(father.id, 1, FUN = function(x) sum(!is.na(x)) )
father.id[is.na(father.id)] = 0

#prep pedigree
A_moms = A[rownames(A) %in% as.character(key.id), colnames(A) %in% as.character(key.id)]
dimnames(A_moms)[[1]] = new.id[match(dimnames(A_moms)[[1]], key.id)]
dimnames(A_moms)[[2]] = new.id[match(dimnames(A_moms)[[2]], key.id)]
A_moms = as.matrix(A_moms[order(as.numeric(row.names(A_moms))), order(as.numeric(colnames(A_moms)))])

#calculate average pedigree relatedness among neighbors
mean_r = as.vector(NA)
for(n in 1:nrow(data.df)){
    if(neigh.n[n]>0){ #0 = not observed, 1 = observed
      ns = neigh.id[n, 1:neigh.n[n]]
      dims = as.character(c(data.df$mom_id[n], ns))
      A_ns = as.matrix(A_moms[which(rownames(A_moms)%in% dims), 
                                 which(colnames(A_moms) %in% dims)])
      mean_r[[n]] = mean(A_ns[lower.tri(A_ns)]) }
    else{mean_r[[n]] = NA
    }}
data.df$mean_r = round(mean_r,3)

#create spatial autocorrelation matrix among communities
comm2 = comm
comm2$com_id = new.com[match(comm$community_id, key.com)]
comm2 = comm2[!is.na(comm2$com_id),]
comm2 = comm2[order(comm2$com_id),]
S = as.matrix(dist(comm2[,c("latitude","longitude")], diag=T, upper=T))

#create list of data for Stan####################################

#offset for years of sampling
#+3 = +1 year of birth, +1 year pre-birth, +1 year post-birth
data.df$offset = data.df$years + 3

#data list
stan_data = list(
  N = nrow(data.df),
  I = nrow(A_moms),
  F = length(key.fth), 
  C = length(key.com),
  M = length(key.mom),
  Y = max(data.df$yearl_id),
  focal_id = data.df$mom_id,
  com_id = data.df$com_id,
  maternal_id = data.df$maternal_id,
  
  father_id = father.id,
  father_n = father.n,
  social_id = neigh.id,
  social_n = neigh.n,
  social_dim = ncol(neigh.id),
  father_dim = ncol(father.id),
  
  yearf_id = data.df$yearf_id,
  yearl_id = data.df$yearl_id,
  yearobs = data.df$offset,
  age = data.df$mom_age,
  mean_r = data.df$mean_r,
  
  A = A_moms,
  S = S,
  fert = data.df$fert
  )
 
#save data
saveRDS(stan_data, "tsi_fert_data.RDS")

#calculate matrix similarity####################################

#relatedness matrix
Am = stan_data$A[unique(stan_data$focal_id),unique(stan_data$focal_id)]

#expand community distance to individual matrix
cl = data.df[,c("mom_id", "com_id")]#remove repeats
cl = cl[!duplicated(cl$mom_id), ]
N = nrow(cl)
D_cid = matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  for (j in 1:N) {
      D_cid[i, j] = stan_data$S[cl$com_id[i], cl$com_id[j]]
    }
}
#distance to similarity
D_cid = 1/(1 + D_cid)

#generate neighborhood adjacency matrix
library(igraph)
neigh = data.frame(cbind(stan_data$focal_id, stan_data$social_id))
colnames(neigh)[1] = "id"
neigh[neigh==0] = NA
dfg = graph_from_data_frame(neigh)
dfga = as_adjacency_matrix(dfg, sparse = FALSE)
dfga = dfga[-nrow(dfga),-ncol(dfga)]
diag(dfga) = 1

#compute similarity values (0 to 1)
library(MatrixCorrelation)
RVadj(Am, D_cid)
RVadj(Am, dfga)
RVadj(D_cid, dfga)

#plot community map####################################

#prep data
csid = unique(stan_data$com_id[stan_data$social_n>0])
cm = comm2[,c("com_id","latitude","longitude")]
cm$soc = as.factor(ifelse(cm$com_id %in% csid, 1, 0))
cml = as.vector(NA)
cmr = as.vector(NA)
for(n in 1:max(stan_data$com_id)){
  cml[[n]] = ifelse(length(unique(as.vector(stan_data$social_id[stan_data$com_id==n,])))==1,0,
                    length(unique(as.vector(cbind(stan_data$focal_id[stan_data$com_id==n &
                                    stan_data$social_n>0],
                      stan_data$social_id[stan_data$com_id==n & 
                                    stan_data$social_n>0,])))) - 1) #remove 0
  cmr[[n]] = mean(stan_data$mean_r[stan_data$com_id==n & stan_data$mean_r!=-99])
}
cm$cml = cml
cm$cmr = ifelse(is.na(cmr),0,cmr)

library(leaflet)
# Create map
lcm = 
  leaflet(cm, options = leafletOptions(attributionControl=FALSE)) %>% 
    addProviderTiles(providers$Esri.NatGeoWorldMap) %>% 
    addCircleMarkers(
        lng = ~longitude,
        lat = ~latitude,
        radius = ~7 + cml/7,
        fillColor = ~colorFactor(c("darkorange", "gold"), soc)(soc),
        color = "black", 
        stroke = T,
        weight = 2,
        fillOpacity = 0.7) %>% 
  addScaleBar(position = "bottomleft",
              options=scaleBarOptions(imperial = F, maxWidth = 200))

#save image
library(mapview)
mapshot(lcm, file = "tsi_map.png", zoom = 2, remove_controls =c("zoomControl"),
        vwidth = 1000, vheight = 900)
