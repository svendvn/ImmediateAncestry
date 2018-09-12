args <- commandArgs(TRUE)
hidden_states_file=args[1]
parameters=args[2]
possible_names=args[3]
individual_name=args[4]
plot_filename=args[5]
#generations=as.numeric(args[4])

short_names=strsplit(possible_names,split='')[[1]]

df=read.csv(hidden_states_file, header=F, stringsAsFactors = F)
#print(df)
colnames(df) <- c('State','Chromosome')

params=strsplit(parameters,split='')[[1]]
m=length(params)/2



#standard values
#long_names=c('Pte','Pts','Ptt','Ptv')
#short_names=c('e','s','t','v')
#df=data.frame(State=c(0,0,0,1,1,3,3,3,3),Chromosome=c(1,1,1,1,1,2,2,2,2))
#parameters='eetv'
#individual_name='test'
#plot_filename='test.pdf'

df$id=0
for(i in unique(df$Chromosome)){
  indices=which(df$Chromosome==i)
  df[indices, 'id']=1:length(indices)
}

get_state=function(state){
  father=state%%m
  mother=floor(state/m)
  return(c(father+1, mother+1, params[m+father+1],params[(mother+1)]))
}

pop_columns=sapply(df$State, get_state, simplify = T)
pop_columns_df=data.frame(t(pop_columns), stringsAsFactors = F)
colnames(pop_columns_df) <- c('Ancestor_1','Ancestor_2','Species_1','Species_2')

N=1000
rectangles=data.frame(xmin=rep(0,N),
                      xmax=rep(0,N),
                      ymin=rep(0,N),
                      ymax=rep(0,N),
                      species=rep('',N),
                      chrom=rep('',N),
                      stringsAsFactors = F) #4 corners of the rectangle, plus species plus chromosome

save_row=function(tier, xmin, xmax,species, chrom, rectangles){
  if(tier=='bottom'){
    row=c(xmin,xmax,0,0.5,species,chrom)
  }
  else if(tier=='top'){
    row=c(xmin,xmax,0.5,1,species,chrom)
 }
  return(rbind(rectangles, row))
}

xmin_top=0
xmin_bottom=0
chrom=df[1,'Chromosome']
top_species=pop_columns_df[1,'Species_1']
bottom_species=pop_columns_df[1,'Species_2']
top_ancestor=pop_columns_df[1,'Ancestor_1']
bottom_ancestor=pop_columns_df[1,'Ancestor_2']
j=1
for(i in 2:(nrow(pop_columns_df))){
  
  if(df[i,'Chromosome']!=chrom){
    xmax_bottom=df$id[i-1]
    xmax_top=df$id[i-1]
    rectangles[j,]=c(xmin_top,xmax_top,0,0.5,top_species,chrom)
    j=j+1
    rectangles[j,]=c(xmin_bottom,xmax_bottom,0.5,1,bottom_species,chrom)
    j=j+1
    xmin_top=0
    xmin_bottom=0
    chrom=df[i,'Chromosome']
    top_species=pop_columns_df[i,'Species_1']
    bottom_species=pop_columns_df[i,'Species_2']
    top_ancestor=pop_columns_df[i,'Ancestor_1']
    bottom_ancestor=pop_columns_df[i,'Ancestor_2']
  }
  else if(pop_columns_df[i,'Ancestor_1']!=top_ancestor){
    xmax_top=df$id[i-1]
    
    rectangles[j,]=c(xmin_top,xmax_top,0,0.5,top_species,chrom)
    j=j+1
    top_species=pop_columns_df[i,'Species_1']
    top_ancestor=pop_columns_df[i,'Ancestor_1']
    xmin_top=df$id[i-1]
  }
  else if(pop_columns_df[i,'Ancestor_2']!=bottom_ancestor){
    xmax_bottom=df$id[i-1]
    
    rectangles[j,]=c(xmin_bottom,xmax_bottom,0.5,1,bottom_species,chrom)
    j=j+1
    bottom_species=pop_columns_df[i,'Species_2']
    bottom_ancestor=pop_columns_df[i,'Ancestor_2']
    xmin_bottom=df$id[i-1]
  }
}
xmax_bottom=df$id[nrow(pop_columns_df)]
xmax_top=df$id[nrow(pop_columns_df)]
rectangles[j,]=c(xmin_top,xmax_top,0,0.5,top_species,chrom)
j=j+1
rectangles[j,]=c(xmin_bottom,xmax_bottom,0.5,1,bottom_species,chrom)
j=j+1
rectangles=subset(rectangles, ymax>0)
rectangles[,1:4]=apply(rectangles[,1:4], c(1,2), as.numeric)

library(RColorBrewer)
#display.brewer.pal(n = length(short_names), name = 'Set2')
cols=brewer.pal(n = length(short_names), name = 'Set2')
names(cols) <- sort(short_names)

#print(rectangles)

library(ggplot2)

ggplot(data=rectangles, aes(xmin=xmin, 
                            xmax=xmax, 
                            ymin=ymin, 
                            ymax=ymax, 
                            group=species))+
  facet_grid(chrom ~.)+theme_bw()+
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        #panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=16))+
  geom_rect(aes(fill=species), color='black', size=0.3)+
  scale_fill_manual(values = cols)+
xlab('locus number')+
  ggtitle(paste(individual_name, parameters))

ggsave(plot_filename,width=8,height=6)

