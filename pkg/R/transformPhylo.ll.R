transformPhylo.ll <- function(y, phy, model=NULL, kappa=NULL, lambda=NULL, delta=NULL, nodeIDs=NULL, branchRates=NULL, cladeRates=NULL) {
	
	switch(model,		  
		   
		   "bm" = {
					vo <- likTraitPhylo(y, phy)
					},
		   
		   "kappa" = {
					transformPhy <- transformPhylo(phy=phy, model="kappa", kappa=kappa)
					},
		   
		   "lambda" = {
					transformPhy <- transformPhylo(phy=phy, model="lambda", lambda=lambda)
					},
		   
		   "delta" = {
					transformPhy <- transformPhylo(phy=phy, model="delta", delta=delta)
					},
					
			"free" = {
					transformPhy <- transformPhylo(phy=phy, model="free", branchRates=branchRates)
					},
		   
		   "clade" = {
					transformPhy <- transformPhylo(phy=phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates)
					}
		   
		   )
	
	return(likTraitPhylo(y=y, phy=transformPhy))
}
			
