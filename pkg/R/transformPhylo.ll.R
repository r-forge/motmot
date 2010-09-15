transformPhylo.ll <- function(y, phy, model=NULL, kappa=NULL, lambda=NULL, delta=NULL, alpha=NULL, psi=NULL, nodeIDs=NULL, rateType=NULL, branchRates=NULL, cladeRates=NULL) {
	
	switch(model,		  
		   
		   "bm" = {
					transformPhy <- phy
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
					transformPhy <- transformPhylo(phy=phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates, rateType=rateType)
					},
		   
		   "OU" = {
					transformPhy <- transformPhylo(phy=phy, model="OU", alpha=alpha)
					},
		   
		   "psi" = {
					transformPhy <- transformPhylo(phy=phy, model="psi", psi=psi)
					}
		   
		   )
	
	return(likTraitPhylo(y=y, phy=transformPhy))
}
			
