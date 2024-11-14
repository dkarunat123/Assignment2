University of Guelph - BINF 6210: Bioinformatics Software Tools

Phylogenetic Signal in Body Mass and Habitat Breadth Traits Among Mustelidae Species

Dilma Karunathilake (0969209)



Introduction: 

This project investigates the relationship between two traits, body size and habitat breadth, within the Mustelidae Family. Mustelidae is composed of a variety of species, including otters, weasels, ferrets, badger, wolverines, and ferrets (Britannica Mustelid, 2024). Understanding whether these traits are evolutionarily conserved within this Family is important as this can provide insight into the ecological and evolutionary pressures that shape species adaptation. This project will explore whether there is a significant correlation between body size and habitat specialization and whether these traits exhibit a phylogenetic signal, indicating evolutionary inheritance. The presence of a phylogenetic signal indicates that species that are closely related tend to exhibit similar trait values specifically because of their shared ancestry (Adams, 2014). To conduct this analysis, Mustelidae species sequence data from the cytochrome c oxidase subunit I (COI) gene, a commonly used DNA barcode marker for animals, combined with trait data from the PanTHERIA database (Jones et al., 2009). DNA sequence analysis methods such as multiple sequence alignment and phylogenetic tree construction are applied alongside statistical tests (Blomberg’s K and Pagel’s λ), to assess phylogenetic signals (Münkemüller et al., 2012).

The objective of this project is to answer the research question: "Are body size and habitat breadth correlated within the Mustelidae family, and do these traits exhibit significant phylogenetic signals?". This will allow us to determine if body size and habitat breadth are phylogenetically conserved and evaluate the relationship between the traits themselves. This project also aims to contribute to understanding evolutionary patterns within the Mustelidae Family and demonstrates how DNA sequence data and phylogenetic comparative methods can address ecological and evolutionary questions.



Results & Discussion: 

In this project, the phylogenetic structure and trait correlation within the Mustelidae Family, focusing on body mass and habitat breadth, was examined. Figure 1, a boxplot of adult body mass, shows the variation in size among Mustelidae species, with a median around 3750 grams and several extreme values, highlighting the diversity in body sizes within this Family. Figure 2, a histogram of habitat breadth, reveals that most species occupy a narrow range of habitats, indicating a trend toward habitat specialization, while a few species demonstrate greater habitat generalism. The scatter plot in Figure 3 shows a positive, albeit weak, correlation between body mass and habitat breadth, suggesting that larger species may have adapted to a broader range of habitats. However, the shallow slope of the regression line indicates that this relationship is not particularly strong, meaning that other ecological or phylogenetic factors may be influencing this association. The phylogenetic analyses used in this project also further support these observations. Figure 4, a phylogenetic tree of 16 species of Mustelidae, illustrates the evolutionary relationships among the species based on COI-5P sequence data. Statistical tests indicate significant phylogenetic signals for both traits. Blomberg’s K test resulted in a K value of 1.37 (p = 0.001) for body mass and 1.23 (p = 0.029) for habitat breadth, showing that both traits are phylogenetically conserved. Pagel’s λ test confirmed these findings as well, with λ values close to 1 for both traits (body mass: λ = 0.999, p = 0.011, habitat breadth: λ = 0.999, p = 0.022), suggesting that closely related species share similar body sizes and habitat breadths. 

While the dataset was sufficient to detect these patterns, a limitation of this study is the modest sample size, as only 16 species had both sequence and trait data available, which may be limiting broader generalizations within the Family. Future work could expand the dataset by including additional Mustelidae species or traits, and possibly reveal further ecological or evolutionary patterns. Additionally, analyzing different gene regions other than COI-5P might provide a more comprehensive understanding of trait evolution in relation to phylogenetic signals within this diverse Family.



References: 

Adams, D. C. (2014). A generalized K statistic for estimating phylogenetic signal from shape and 
other high-dimensional multivariate data. Systematic Biology, 63(5), 685–697. https://doi.org/10.1093/sysbio/syu030

Aly S. M. (2014). Reliability of long vs short COI markers in identification of forensically 
important flies. Croatian medical journal, 55(1), 19–26. https://doi.org/10.3325/cmj.2014.55.19

Britannica. (2024). Mink. In Encyclopedia Britannica. Retrieved October 23, 2024, from 
https://www.britannica.com/animal/mink#ref158576

Britannica. (2024). Mustelid. In Encyclopedia Britannica. Retrieved October 23, 2024, from 
https://www.britannica.com/animal/mustelid

Jones, K. E., Bielby, J., Cardillo, M., Fritz, S. A., O'Dell, J., Orme, C. D. L., Safi, K., Sechrest, 
W., Boakes, E. H., Carbone, C., Connolly, C., Cutts, M. J., Foster, J. K., Grenyer, R., Habib, M., Plaster, C. A., Price, S. A., Rigby, E. A., Rist, J., Teacher, A., ... Purvis, A. (2009). PanTHERIA: A species-level database of life history, ecology, and geography of extant and recently extinct mammals. Ecological Archives, E090-184. https://doi.org/10.1890/08-1494.1

Münkemüller, T., Lavergne, S., Bzeznik, B., Dray, S., Jombart, T., Schiffers, K., & Thuiller, W. 
(2012). How to measure and test phylogenetic signal. Methods in Ecology and Evolution, 3(6), 743–756. https://doi.org/10.1111/j.2041-210X.2012.00196.x

Nugent, C. M., Elliott, T. A., Ratnasingham, S., & Adamowicz, S. J. (2020). COI-PL: An R 
package for cytochrome c oxidase I (COI) DNA barcode data cleaning, translation, and error evaluation. Genome, 63(6), 291–305. https://doi.org/10.1139/gen-2019-0206

Revell, L. J., Harmon, L. J., & Collar, D. C. (2008). Phylogenetic signal, evolutionary process, 
and rate. Systematic biology, 57(4), 591–601. https://doi.org/10.1080/10635150802302427

Schwartz, S., Elnitski, L., Li, M., Weirauch, M., Riemer, C., Smit, A., NISC Comparative 
Sequencing Program, Green, E. D., Hardison, R. C., & Miller, W. (2003). MultiPipMaker and supporting tools: Alignments and analysis of multiple genomic DNA sequences. Nucleic acids research, 31(13), 3518–3524. https://doi.org/10.1093/nar/gkg579

Package References:

Pagès, H., Aboyoun, P., Gentleman, R., & DebRoy, S. (2024). Biostrings: Efficient manipulation 
of biological strings (R package version 2.72.1). Bioconductor. https://doi.org/10.18129/B9.bioc.Biostrings

Paradis, E., & Schliep, K. (2019). ape 5.0: An environment for modern phylogenetics and 
evolutionary analyses in R. Bioinformatics, 35, 526-528. https://doi.org/10.1093/bioinformatics/bty633

Revell, L. J. (2024). phytools 2.0: An updated R ecosystem for phylogenetic comparative 
methods (and other things). PeerJ, 12, e16505. https://doi.org/10.7717/peerj.16505

Schliep, K. P. (2011). phangorn: Phylogenetic analysis in R. Bioinformatics, 27(4), 592-593.

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L. D., François, R., Grolemund, G., 
Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T. L., Miller, E., Bache, S. M., Müller, K., Ooms, J., Robinson, D., Seidel, D. P., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K., & Yutani, H. (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686. https://doi.org/10.21105/joss.01686

Winter, D. J. (2017). rentrez: An R package for the NCBI eUtils API. The R Journal, 9(2), 520-
526.

Wright, E. S. (2016). Using DECIPHER v2.0 to analyze big biological sequence data in R. The R 
Journal, 8(1), 352-359.
