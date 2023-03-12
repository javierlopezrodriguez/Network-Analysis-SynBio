# Network-Analysis-SynBio
Network Analysis exercise from the Synthetic and Systems Biology course of UPM's Master in Computational Biology

The task: Considering E. coli's regulatory network (transcription factors versus genes) from RegulonDB, and comparing it with an appropiate random network, explain the following points:
- are autoregulation motifs overrepresented?
- are feed forward loops overrepresented?
- is the length of the longest cascade shorter or longer than expected at random?

The code is written in R, mainly using base R (instead of the igraph library). Only the script files are included (no data or results). The files are the following:
- Network.R: main script.
- NetworkFunctions.R: Implementation of the algorithms for finding and classifying autoregulation motifs, finding and classifying feed forward loops, and finding the length of the longest path using base R. There are also other helper functions.
    - Autoregulation motifs: distinguishing between positive and negative loops.
    - Feed forward loops: distinguishing between the four coherent and four incoherent loops.
- NetworkFunctionsTest.R: Tests for the main functions of NetworkFunctions.R
- ResultsPlotting.R: Generating result plots.
- StatisticalAnalysis.R: Analyze the statistical significance of the results.

Disclaimer: This is a very crude implementation of those algorithms, and there is a large room for improvement. Don't use this in any meaningful task without checking the igraph library first (which can do everything I do here in a more optimized fashion, and much more. Available here: https://igraph.org/r/)

Data used: RegulonDB's Escherichia coli Regulatory Network interactions (TF-Gene), Release: 10.9 Date: 06/29/2021, downloaded on October 12, 2021.


References used on this project:

Wong, E., Baur, B., Quader, S., & Huang, C. H. (2012). Biological network motif detection: principles and practice. Briefings in bioinformatics, 13(2), 202-215.

RegulonDB v 10.5: tackling challenges to unify classic and high throughput knowledge of gene regulation in E. coli K-12. Alberto Santos-Zavaleta, Heladia Salgado, Socorro Gama-Castro, Mishael Sánchez-Pérez, Laura Gómez-Romero, Daniela Ledezma-Tejeida, Jair Santiago García-Sotelo, Kevin Alquicira-Hernández, Luis José Muñiz-Rascado, Pablo Peña-Loredo, Cecilia Ishida-Gutiérrez, David A Velázquez-Ramírez, Víctor Del Moral-Chávez, César Bonavides-Martínez, Carlos-Francisco Méndez-Cruz, James Galagan, Julio Collado-Vides. Nucleic Acids Research, Volume 47, Issue D1, 08 January 2019, Pages D212?D220, doi: https://doi.org/10.1093/nar/gky1077

Madar, D., Dekel, E., Bren, A., & Alon, U. (2011). Negative auto-regulation increases the input dynamic-range of the arabinose system of Escherichia coli. BMC systems biology, 5(1), 1-9.

Jin G. (2013) Feed Forward Loop. In: Dubitzky W., Wolkenhauer O., Cho KH., Yokota H. (eds) Encyclopedia of Systems Biology. Springer, New York, NY. https://doi.org/10.1007/978-1-4419-9863-7_463

Ren X. (2013) Transcriptional Cascade. In: Dubitzky W., Wolkenhauer O., Cho KH., Yokota H. (eds) Encyclopedia of Systems Biology. Springer, New York, NY. https://doi.org/10.1007/978-1-4419-9863-7_471

Sedgewick, Robert; Wayne, Kevin Daniel (2011), Algorithms (4th ed.), Addison-Wesley Professional, pp. 661–666, ISBN 9780321573513.

Ioannidou, K., Mertzios, G. B., & Nikolopoulos, S. D. (2009, August). The longest path problem is polynomial on interval graphs. In International Symposium on Mathematical Foundations of Computer Science (pp. 403-414). Springer, Berlin, Heidelberg.

stackoverflow: How to find the longest simple path in a graph? https://stackoverflow.com/a/21880918

Mangan, S., Itzkovitz, S., Zaslaver, A., & Alon, U. (2006). The incoherent feed-forward loop accelerates the response-time of the gal system of Escherichia coli. Journal of molecular biology, 356(5), 1073-1081.