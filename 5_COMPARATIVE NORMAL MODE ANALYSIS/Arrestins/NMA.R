# Loading libraries -------------------------------------------------------

library(bio3d)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)

# Functions to extract information from PDB files and build the graphs ----

process_pdb <- function(pdb){
 print(pdb)
 pdb_file <- read.pdb(pdb)
 chains <- unique(chain.pdb(pdb_file))[1:2]
 name <- substr(x = pdb,start = 1,stop = str_length(pdb) - 4)
 name_vec <- strsplit(x = name,split = '_')[[1]]
 output <- c(pdb,name_vec,chains)
 return(output)
}

plot_nma <- function(monomer1,unique_mon,all_pdbs,chains,monomer2) {
 warnings_current <- c()
 for (mon in unique_mon) {
  monomer2 <- monomer2[monomer1 == mon]
  mon_df <- all_pdbs[monomer1 == mon,] %>%
   data.frame()
  mon_df$file <- mon_df$file %>%
   as.character()
  if(length(mon_df$file) > 1){
   chains <- chains %>%
    as.character
   pdb_list <- lapply(mon_df$file,FUN = read.pdb)
   trimmed_list <- list()
   interface_list <- list()
   interface_resno <- list()
   for (i in 1:length(pdb_list)) {
    trimmed_list[[i]] <- trim.pdb(pdb_list[[i]], 
                                  atom.select(pdb_list[[i]], chain = chains[i]))
    interface_list[[i]] <- binding.site(pdb_list[[i]], 
                                        a.inds = atom.select(pdb_list[[i]], 
                                                             chain = chains[i]),
                                        b.inds = atom.select(pdb_list[[i]], 
                                                             chain = chains[i],
                                                             inverse = TRUE))
    interface_resno[[i]] <- interface_list[[i]][3]
    to_write <- interface_list[[i]][2]
    write.csv(x = to_write,
              file = paste('interface_',mon,i,'.csv',sep=''))
   }
   
   interface_all <- unlist(interface_resno[[1]])
   seq_list <- lapply(trimmed_list, FUN = pdbseq)
   
   all_seq <- do.call(seqbind,seq_list)
   aligned_seq <- seqaln(all_seq,exefile = 'muscle.exe',id = mon_df$file)
   aligned_fasta <- read.fasta.pdb(aligned_seq)
   all_modes <- nma(aligned_fasta,full=TRUE)
   
   fluctuations <- all_modes$fluctuations %>%
    t() %>%
    data.frame()
   
   interface_all <- interface_all[interface_all < nrow(fluctuations)]
   
   colnames(fluctuations) <- monomer2
   fluctuations$num <- seq(nrow(fluctuations))
   write.csv(x = fluctuations,file = "fluctuations.csv")
   fluctuations_long <- melt(fluctuations,id.vars = "num")
   fluctuations_long$Monomer <- fluctuations_long$variable
   
   name <- paste(mon,'.png')
   png(filename = paste('default_',name,sep=''))
   plot(all_modes)
   dev.off()
   nma_graph <- ggplot(data = fluctuations_long,
                       mapping = aes(x = num, y = value,
                                     colour = Monomer)) +
    geom_line(size = 0.6,alpha=0.8) + 
    geom_vline(xintercept = interface_all,alpha = 0.5,size = 0.25,
               show.legend = TRUE) + 
    theme_classic() + 
    xlab("Residue Number") + 
    ylab("Fluctuation") + 
    scale_x_continuous(breaks = seq(0,1000,by = 50),
                       minor_breaks = seq(0,1000,by = 10),
                       limits = c(0,nrow(fluctuations))) + 
    coord_cartesian(xlim = c(0,nrow(fluctuations)))
   
   ggsave(filename = name,
          plot = nma_graph,height = 4,width = 11.69)
  }
  else {
   warning <- paste(mon, 'has only one partner. No comparative analysis can be done')
   warnings <- c(warnings_current,warning)
  }
 }
 if(length(warnings_current) > 0) {
  print(warnings_current)
 }
}

# Running the script ------------------------------------------------------

all_files <- list.files(path = '.',pattern = '*\\.pdb')

all_pdbs <- sapply(all_files,FUN = process_pdb) %>%
 t() %>%
 data.frame()
colnames(all_pdbs) <- c('file','mon1','mon2','chain1','chain2')
rownames(all_pdbs) <- NULL
write.csv(all_pdbs,file = 'all_pdb.csv')

unique_mon1 <- unique(all_pdbs$mon1)
unique_mon2 <- unique(all_pdbs$mon2)

plot_nma(all_pdbs$mon1,unique_mon1,all_pdbs,all_pdbs$chain2,all_pdbs$mon2)
plot_nma(all_pdbs$mon2,unique_mon2,all_pdbs,all_pdbs$chain1,all_pdbs$mon1)
