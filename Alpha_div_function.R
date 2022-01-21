##alpha div function 

alpha_div <- function(phylo, depth, trials){
  ##adopted from 'http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#alpha_diversity'
  nsamp = nsamples(phylo)
  richness <- matrix(nrow = nsamp, ncol = trials)
  row.names(richness) <- sample_names(phylo)
  ShanIn <- matrix(nrow = nsamp, ncol = trials)
  row.names(ShanIn) <- sample_names(phylo)
  set.seed(4200)
  
  for (i in 1:trials) {
    # Subsample
    r <- rarefy_even_depth(phylo, sample.size = depth, verbose = FALSE, replace = TRUE)
    
    # Calculate richness
    rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
    richness[ ,i] <- rich
    
    # Calculate ShanIn
    shan <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
    ShanIn[ ,i] <- shan
  }
  
  # Create a new dataframe to hold the means and standard deviations of richness estimates
  SampleID <- row.names(richness)
  mean <- apply(richness, 1, mean)
  sd <- apply(richness, 1, sd)
  measure <- rep("Observed", nsamp)
  rich_stats <- data.frame(SampleID, mean, sd, measure)
  
  # Create a new dataframe to hold the means and standard deviations of ShanIn estimates
  SampleID <- row.names(ShanIn)
  mean <- apply(ShanIn, 1, mean)
  sd <- apply(ShanIn, 1, sd)
  measure <- rep("Shannon", nsamp)
  shan_stats <- data.frame(SampleID, mean, sd, measure)
  
  alpha <- rbind(rich_stats, shan_stats)
  s <- data.frame(sample_data(phylo))
  alphadiv <- merge(alpha, s, by = "SampleID") 
  
  return(alphadiv)
  
}