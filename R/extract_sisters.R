extract_sisters <-
function(tree, sis_age=FALSE, mol_clock=NULL, crown_age=NULL) {
  # find which nodes lead to tips
  tip_node = tree$edge[tree$edge[,2] <= length(tree$tip.label), 1]
  # of these, which appear exactly twice
  pair_node = unique(tip_node[sapply(tip_node, function(y) sum(tip_node == y) == 2)])
  # find one of the two indices for each of these "cherry" nodes
  edge.index=sapply(pair_node, function(x) which(tree$edge[,1]==x))[1,]
  # find the taxa names for each pair
  sp_names = t(sapply(edge.index, function(x) {
  tips = tree$edge[tree$edge[,1]==tree$edge[x,1],2]
  c(tree$tip.label[tips[1]], tree$tip.label[tips[2]])
  }))
  # Estimate ages of pairs
  if(sis_age==TRUE) {
  	warning("Just a reminder to please be aware: extract_sisters assumes your tree is ultrametric")
  	if(is.null(mol_clock)==FALSE & is.null(crown_age)==FALSE) {
  	stop("Please provide either a molecular clock OR a crown age 
  	OR neither, but not both")
    }
    # find the branch length corresponding to each index and calculate the age of the pair
    if(is.null(mol_clock) & is.null(crown_age)) {
  	  pair_ages = cbind(tree_node=pair_node, pair_age=tree$edge.length[edge.index])
    }
    if(is.null(mol_clock)==FALSE) {
      pair_ages = cbind(tree_node=pair_node, pair_age=tree$edge.length[edge.index]/mol_clock)
    }
    if(is.null(crown_age)==FALSE) {
  	  # get crown node index
  	  crown = min(tree$edge[,1])
  	  # get the index of all the branches leading to the crown nodes (up to 50 branches)
  	  row1 = which(tree$edge.length==max(tree$edge.length))
      if(tree$edge[row1,1]==crown) {
        cr_length=tree$edge.length[row1]
      } else {
        rows=vector()
        rows[1]=row1
        for(i in 2:50) {
  	      row = which(tree$edge[,2]==tree$edge[rows[(i-1)],1])
  	      if(length(row)>0) rows[i] = row
  	      if(length(row)==0) rows[i] = NA
        }
        # sum the branch lengths from a tip to the node
        cr_length = sum(tree$edge.length[rows[is.na(rows)==FALSE]])
      }
      # calculate mol.clock rate based on crown age and crown branch lengths
  	  rate = cr_length/crown_age
  	  # calcualte age of each pair based on branch lengths and this rate
      pair_ages = cbind(tree_node=pair_node, pair_age=tree$edge.length[edge.index]/rate)
    }
  # wrap everything up in a dataframe
  res = data.frame(sp1 = sp_names[,1], sp2 = sp_names[,2], pair_ages, stringsAsFactors=FALSE)
  } else {
  # wrap everything up in a dataframe
  res = data.frame(sp1 = sp_names[,1], sp2 = sp_names[,2], stringsAsFactors=FALSE)
  }
  return(res)
}