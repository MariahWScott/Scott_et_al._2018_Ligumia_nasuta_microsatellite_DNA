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
  #Retained informative PCs (we used 100), all eigenvalues retained (so 8 to be safe)
  compoplot(dapc1,
            txt.leg=paste("Cluster", 1:9), lab="",
            ncol=1, xlab="individuals", col=rainbow(9))
  #Multiple Scatter Plots of different formats
  scatter(dapc1, posi.da="bottomright", bg="white", pch=20, cstar=0, 
          col=rainbow(9), scree.da=TRUE, posi.pca="bottomleft", cell=0, 
          cex=1.5, clab=0, txt.leg=paste("Cluster",1:9), scree.pca=TRUE)
  scatter(dapc1, posi.da="bottomright", bg="white",
          pch=20, cstar=0, col=rainbow(9), scree.pca=TRUE,
          posi.pca="bottomleft", clab=0.50)
  scatter(dapc1, pch=20, cstar=0, col=rainbow(9), clab=0.5)
  #Need to make number of colors in scatter plot match the number of PC clusters used
  return(print("Done!"))
}

#Run the above function for each file desired, with the number of maximum clusters defined
#Typically you will want the number of maximum clusters the number of sites sampled + 1
#Look at the comments in the above function for suggestions on how to respond to the prompts in the console
define_and_plot_cluster(file_name = "GenePop_April_9_K_is_57.gen", n_max_clusters = 41)
