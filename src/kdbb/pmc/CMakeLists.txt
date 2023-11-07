add_library(PMC pmc.h
	pmc_headers.h
	pmc_input.h
	pmc_heu.h pmc_heu.cpp
	pmc_graph.h pmc_graph.cpp pmc_cores.cpp pmc_clique_utils.cpp
	pmc_utils.h pmc_utils.cpp
	pmc_maxclique.h pmc_maxclique.cpp
	pmcx_maxclique_basic.h pmcx_maxclique_basic.cpp
	pmcx_maxclique.h pmcx_maxclique.cpp
)


