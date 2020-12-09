# PLeXuS
PLeXuS - Protein Learning Exuberant Software. This is an app collecting simple games helping you to learn protein-protein interactions easier. It is heavily built upon Igraph package and has a visual interface made in Tkinter. 

Learn protein interactions in delineated signaling cascades from (KEGG, SIGNOR, Panther or reactome) or load your own gene sets and and highlight connections between them as stored in the comprehensive graph collected by <a href="https://github.com/saezlab/pypath"> Omipath</a>. 

<img src="https://github.com/culpritgene/PLeXuS/blob/master/Resources/PLEXUS_Sreenshot.png" width="620" height="700" />

Currently available modes

<b>Labirinth game</b>
Pave the path through a mesh of proteins selecting only edges supported by the biological data; help yourself using snippets of functional gene annotations from Uniprot or Ensemble masked on different levels. Store your corrent and incorrect guesses, building personal statstics. 

<b>Mines game</b>
Organized in the similar way as Labirinth game, but the objective is reversed - highlight only genes which do not have protein interaction with a given gene. 
