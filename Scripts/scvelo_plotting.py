import session_info
import wget
import scvelo as scv


#>>> session_info.show()
#-----
#scvelo              0.2.3
#session_info        1.0.0
#wget                3.2
#-----
#Python 3.8.5


filename = wget.download("https://data.bioinfo.vbc.ac.at/knoblich.grp/OEichmueller/TSC2/OEichmueller.TSC2.loom")

adata = scv.read(filename, cache=False)
adata.obs['Cluster'] =  adata.obs['Cluster'].apply(str)

scv.settings.autosave = True
scv.settings.figdir = "OEichmueller_velocity/"
scv.settings.set_figure_params(format="pdf")


# Filtering, normalization and log transform
scv.pp.filter_and_normalize(adata)

# Computing first-/second-order moments per cell
scv.pp.moments(adata)

# Estimating gene-specific velocities using stochastic model of transcriptional dynamics
scv.tl.velocity(adata, mode='stochastic')

# Computing velocity graph based on cosine similarities
scv.tl.velocity_graph(adata)


# Plotting of velocities on the embedding
scv.pl.velocity_embedding_stream(adata, color="Cluster", basis='umap', save= ("static.scv.pl.velocity_embedding_stream.3d"),projection="3d",palette= { "1": '#F8766D', "2": '#EA8331', "3": '#D89000', "4": '#C09B00', "5": '#A3A500', "6": '#7CAE00', "7": '#39B600', "8": '#00BB4E', "9": '#00BF7D', "10": '#00C1A3', "11": '#00BFC4', "12": '#00BAE0', "13": '#00B0F6', "14": '#35A2FF', "15": '#9590FF', "16": '#C77CFF', "17": '#E76BF3', "18": '#FA62DB', "19": '#FF62BC', "20": '#FF6A98' } )

# Plotting spliced against unspliced expressions for EGFR
# Plotting velocity and expression on the embedding
scv.pl.velocity(adata,var_names="EGFR",color="velocity",colorbar=True,color_map="RdBu_r",save= ("scv.pl.velocity.EGFR"))
