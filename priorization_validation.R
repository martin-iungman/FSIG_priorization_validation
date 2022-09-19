#package installation and loading
packages<-c( "tidyverse", "pedprobr","pedtools", "pedmut", "forrel", "poibin","pracma")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

#functions
#Remove genotype
comb_familias<-function(vector,len.max=NULL,len.min=NULL){
  
  repet<-abs(pascal(length(vector)+1,1))[length(vector)+1,]
  lista<-list()
  for (i in seq(1, length(vector)+1)) {
    lista1<-list()
    
    while (length(lista1)<repet[i]) {
      v1<-sort(sample(vector, i-1, replace=F))
      lista1[[length(lista1)+1]]<-v1
      lista1<-unique(lista1)
    }
    lista<-c(lista,lista1)
    
  }
  if(!is.null(len.max)){
    lista<-keep(lista,~length(.x)<=len.max)
  }
  if(!is.null(len.min)){
    lista<-keep(lista,~length(.x)>=len.min)
  }
  return(lista)
}
remove_genotype<-function(ped, id){
  for(marker in name(ped)){
    genotype(ped,id=id, marker=marker)<-"-/-"
  }
  return(ped)
}
remove_genotype_mult<-function(ped,ids){
  ped2=ped
  for (id in ids){
    ped2<-remove_genotype(ped2,id)
  }
  return(ped2)
}
modifyID<-function(ped){
  ped$ID<-ped$ID%>%gsub("-.*?$","",.)
  return(ped)
}
remove_partial_marker<-function(ped){
  gtp<-getGenotypes(ped,typedMembers(ped))
  rm_mark<-colnames(gtp)[ceiling(which(gtp=="-/-")/nrow(gtp))]
  ped=removeMarkers(ped, rm_mark)
  if(is_empty(rm_mark)){rm_mark="none"}
  print(paste("Removing markers:", paste(unique(rm_mark), collapse = " ")))
  return(ped)
}
get_IP_table<-function(sim,nProfiles,sum_table=T,lim_inf,lim_sup){
  a<-flatten(sim)%>%flatten()
  a<-a[names(a)=="ip"]
  b<-flatten(a)%>%map(., ~.x$meanLogLR)
  b<-data.frame(LR_indiv=unlist(b), group=ceiling(seq(1,length(b))/nProfiles))
  if(sum_table){
  c<-b%>%group_by(group)%>%summarise(mediaLogLR=mean(LR_indiv), logLR_ICinf=sort(LR_indiv)[lim_inf],logLR_ICsup=sort(LR_indiv)[lim_sup])
  return(c)}else(return(b))
  
}

get_ExpMM_table<-function(sim,nProfiles,sum_table=T,lim_inf,lim_sup){
  a<-flatten(sim)%>%flatten()
  a<-a[names(a)=="ep"]
  b<-flatten(a)%>%map(., ~.x$expectedMismatch)
  b<-data.frame(ExpMM_indiv=unlist(b), group=ceiling(seq(1,length(b))/nProfiles))
  if(sum_table){
  c<-b%>%group_by(group)%>%summarise(mediaExpMM=mean(ExpMM_indiv), ExpMM_ICinf=sort(ExpMM_indiv)[lim_inf],ExpMM_ICsup=sort(ExpMM_indiv)[lim_sup])
  return(c)}else(return(b))
}

fam=readFam("priorization_families.fam", useDVI = T, verbose=F)
xref<-map(fam, ~.x[['Reference pedigree']])

#families to test
nProfiles=50
lrsims=1000
sim_flia_list<-list()
df_flia<-list()
names(df_flia)=names(sim_flia_list)<-names(xref)
list_simID<-list() #each element must be a character vector with the IDs to be simulated for each family
list_tipID<-list() #each element must be a character vector with the IDs of the genotyped members of the baseline family

for (i in families){
  flia<-xref[[i]]
  flia<-remove_partial_marker(flia)
  flia_simID<-list_simID[[i]]
  flia_tipID<-list_tipID[[i]]
  flia_rm<-flia%>%remove_genotype_mult(typedMembers(flia)[!typedMembers(flia)%in%c(flia_simID, flia_tipID)])
  comb_flia<-comb_familias(flia_simID)%>%keep(~length(.x)>=length(flia_simID)-1)
  flia_peds<-map(comb_flia,~remove_genotype_mult(flia_rm, .x))
  sims<-map(flia_peds,~if(length(typedMembers(.x))==(length(flia_tipID)+1)){list(typedMembers(.x))} else {append(list(flia_tipID),lapply(flia_simID, function(x){c(x,flia_tipID)}))})
  sim_flia_list[[i]]<-map2(flia_peds, sims, ~MPPsims(.x,missing = "Missing person", selections = .y, nProfiles = nProfiles, lrSims = 1000, ep=T,ip=T,addBaseline = F, numCores =20, verbose = T))
  names(sims)<-map_chr(flia_peds,  ~typedMembers(.x)%>%str_trim()%>%sort()%>%paste(collapse = " "))
  for(j in 1:length(sims)){if(is_empty(sims[[j]])){sims[[j]]=""}}
  tip<-rep(map_chr(flia_peds, ~typedMembers(.x)%>%str_trim()%>%sort()%>%paste(collapse=" ")), times=map_dbl(sim_flia_list[[i]],~length(.x)))
  tip_n<-rep(map_dbl(flia_peds, ~length(typedMembers(.x))), times=map_dbl(sim_flia_list[[i]],~length(.x)))
  simulados_<-map2(tip,flatten(sims)%>%map(str_trim),~.y[!.y%in%unlist(strsplit(.x,split=" "))])
  simulados<-simulados_%>%sapply(function(x){if(is_empty(x)){""}else{x}})
  simulados_n<-simulados_%>%sapply(function(x){if(is_empty(x)){0}else{length(x)}})
  flia_ip<-get_IP_table(sim_flia_list[[i]],nProfiles = nProfiles,lim_inf = 2,lim_sup=49)
  flia_ep<-get_ExpMM_table(sim_flia_list[[i]],nProfiles = nProfiles,lim_inf = 2,lim_sup=49)
  df_flia[[i]]<-data.frame(tip,tip_n,simulados,simulados_n)%>%cbind(full_join(flia_ip,flia_ep, by="group"))
  df_flia[[i]]<-df_flia[[i]]%>%mutate(tot_indiv=ifelse(simulados_n>0,paste(tip,simulados)%>%strsplit(split=" ")%>%map_chr(~sort(.x)%>%paste(collapse=" ")),tip))%>%
    mutate(is_sim=ifelse(simulados_n>0, "Simulated", "Observed"))
  tiff(paste0(i,".tiff"),width=1024,height=1024)
  ggplot(df_flia[[i]], aes(x=mediaExpMM, y=mediaLogLR, width=ifelse(is_sim=="Simulated",0.20,0.1),col=tot_indiv, shape=as_factor(is_sim),stroke=2,xmin=ExpMM_ICinf, xmax=ExpMM_ICsup, ymin=logLR_ICinf,ymax=logLR_ICsup))+
    geom_point(size=3)+geom_errorbarh()+geom_errorbar()+scale_shape_manual(values=c("Simulated"=21,"Observed"=24))+
    ylab("log10(LR)")+xlab("Expected Mismatch")+theme_bw()+
    theme(legend.title = element_blank(),text = element_text(size=25))+ggtitle(i)
  
  dev.off()     
  
}

