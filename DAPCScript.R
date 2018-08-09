#This code performs Discriminant Analysis of Principal Components (DAPC), using adegenet 2.0.0
#This code finds the most supported number of clusters in the sample group and generates an assignment plot
#similar to STRUCTURE plots.  You will likely have to change the file and folder names and pathways,
#But can use the files and folders provided as an example.

#Install Packages if necessary and Libraries, coded below
library(adegenet)

#Set your working directory, which likely has a different name and pathway
setwd("~/LINA/Lina Structure without less than 10")

define_and_plot_cluster <- function(file_name, n_max_clusters){
  #Define your file, which has alleles with 3 digits (ncode argument)
  obj1 <- import2genind(file_name, ncode=3)
  obj1
  group <- find.clusters(obj1, max.n.clust = n_max_clusters)
  #Retain 200 PCs (or over all) and choose the number of clusters on the elbow of the BIC graph
  dapc1 <- dapc(obj1, group$grp)
  #Retained informative PCs, all eigenvalues retained (so 15 to be safe)
  compoplot(dapc1,
            txt.leg=paste("Cluster", 1:9), lab="",
            ncol=1, xlab="individuals", col=rainbow(9))
  return(print("Done!"))
}

#Run the above function for each file desired, with the number of maximum clusters defined
#Typically you will want the number of maximum clusters the number of sites sampled + 1
#Look at the comments in the above function for suggestions on how to respond to the prompts in the console
define_and_plot_cluster(file_name = "Example_file.gen", n_max_clusters = 41)
