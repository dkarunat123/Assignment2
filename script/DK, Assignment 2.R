## 1. Installing/Loading Script Dependencies: ----

# Installing packages if not already installed (commented out)
# Uncomment if packages are not installed in your environment
# install.packages("tidyverse")
# install.packages("rentrez")
# BiocManager::install("Biostrings")
# BiocManager::install("DECIPHER")
# install.packages("ape")
# install.packages("phangorn")
# install.packages("phytools")
# install.packages("leaflet")

# Loading required libraries
library(tidyverse)
library(rentrez)
library(Biostrings)
library(DECIPHER)
library(ape)
library(phangorn)
library(phytools)
library(styler)
library(leaflet)
library(dplyr)
library(viridis)


## 2a. Searching For & Fetching NCBI Sequence Data: ----

# I will search for and fetch only sequences that fall within the typical COI-5P length remain.

# The main marker that is used in animal DNA barcoding is a fragment of the 5' end of COI (cytochrome c oxidase subunit I) gene called COI-5P (Nugent et al., 2020). The length of this fragment is approximately 657 base pairs long (Nugent et al., 2020). So the sequences I will keep will be in the 600-700 bp range.

# Defining the search query for Mustelidae COI sequences.
Search_Query_Test <- entrez_search(
  db = "nuccore",
  term = "Mustelidae[ORGN] AND COI[GENE] AND 600:700[SLEN]"
)

# Checking how many sequences were found.
Search_Query_Test$count # There are 187 hits in total (which matches the filtering done in the original script).
maxHits <- Search_Query_Test$count

# Now we know that there are 187 hits in total. However, this is too many sequences to fetch at once (Error -> HTTP failure 414, the request is too large), so we have to store the search results on NCBI's servers and retrieve them in using both the use_history argument in entrez_search() and the web_history argument in entrez_fetch().
# entrez_search returns an object with a web_history field that contains both WebEnv and QueryKey so we don't need to include them in the code separately.


# Searching with the argument 'use_history'.
Search_Query <- entrez_search(
  db = "nuccore",
  term = "Mustelidae[ORGN] AND COI[GENE] AND 600:700[SLEN]",
  use_history = TRUE,
  retmax = maxHits
)

# Fetching sequences in batches using the argument 'web_history'.
Batch_Size <- 100

for (seq_start in seq(1, maxHits, Batch_Size)) {
  Mustelidae_COI_fetch <- entrez_fetch(
    db = "nuccore",
    web_history = Search_Query$web_history,
    rettype = "fasta",
    retmax = Batch_Size,
    retstart = seq_start - 1
  )
  cat(Mustelidae_COI_fetch, file = "../data/MustelidaeCOI.fasta", append = TRUE)
  cat(seq_start + Batch_Size - 1, "sequences downloaded\r")
}

Mustelidae_Sequences <- readDNAStringSet("../data/MustelidaeCOI.fasta", format = "fasta")
Mustelidae_Sequences

# To ensure my data was filtered out correctly, I conducted a check with the filtered data from the original code (this has now been deleted in my updated script):
# all.equal(mean(df_Seq_Filtered$Sequence_Length), mean(width(Mustelidae_Sequences)))
# Output: TRUE


# Reading the FASTA data back from the file using readDNAStringSet from the Biostrings package.
Mustelidae_Sequences <- readDNAStringSet("../data/MustelidaeCOI.fasta", format = "fasta")

# Checking that the command worked.
Mustelidae_Sequences
class(Mustelidae_Sequences) # Yes, a DNAStringSet object was created.



## 2b. Getting Rid of Duplicate Species Sequence: ----

# First, we need to get the names (genus, species) associated with each sequence in the data set.

# Putting the sequences from Mustelidae_Sequences into a dataframe. The columns will consists of the sequence header (Title) and the sequence (Sequence).
df_Mustelidae_Seq <- data.frame(Title = names(Mustelidae_Sequences), Sequence = paste(Mustelidae_Sequences))
View(df_Mustelidae_Seq)

# Making a new column called 'Species_Name' and adding it to the dataframe. The code will extract the genus (2L) and species (3L) from the Title column and put it into this new column.
df_Mustelidae_Seq$Species_Name <- word(df_Mustelidae_Seq$Title, 2L, 3L)

# Checking to make sure that this new column has appeared in the data frame and holds the genus and species names for each sequence.
View(df_Mustelidae_Seq)
df_Mustelidae_Seq$Species_Name


# Since there are duplicates (more than one available sequence file belonging to a single species), we need to choose one sequence per species for our downstream analysis.

# I will first select the longest available sequence for each species. This is because Neighbor-Joining Phylogenetic Trees (like the one that I will create later in the script) cluster species in accordance with their taxonomic classification better when longer fragments are used, compared to when shorter fragments are used (Aly, 2014).

# Obtaining the lengths of each sequence and placing this information in a new column called Sequence_Length in the data frame df_Mustelidae_Seq.
df_Mustelidae_Seq$Sequence_Length <- nchar(df_Mustelidae_Seq$Sequence)
View(df_Mustelidae_Seq)

# Grouping by Species_Name and selecting the longest sequence for each species.
df_Species.Seq_Unique <- df_Mustelidae_Seq %>%
  group_by(Species_Name) %>%
  slice_max(order_by = Sequence_Length, n = 1, with_ties = FALSE) %>%
  ungroup()
# slice_max() is used to select the rows with the maximum value in Sequence_Length. Its argument n = 1 makes sure that only one row is selected for each species. with_ties = FALSE makes sure that if there are ties (i.e., multiple sequences with the same maximum length), only one row is selected (arbitrarily).

# Checking that there is only one row per species now.
View(df_Species.Seq_Unique)


# "Neovison vison" and "Neogale vison" are both used to refer to the American Mink (Britannica, 2024). Therefore, I will choose the one with the longer sequence, which is Neogale vison (687 > 657).

# Removing the row belonging to "Neovison vison".
df_Species.Seq_Unique <- df_Species.Seq_Unique %>%
  filter(Species_Name != "Neovison vison")

# Checking that "Neovison vison" is no longer a row in the data frame.
View(df_Species.Seq_Unique)


# Getting the number of unique species & accompanying sequences in the cleaned data set:
length(df_Species.Seq_Unique$Species_Name)
# There are now 20 Mustelidae species in this data set with their accompanying sequences for COI-5P.



## 3a. Reading-in & Looking at the PanTHERIA Traits data: ----

# This data was obtained from Courselink from the BINF6210 lecture on 26th September 2024 and is originally derived from PanTHERIA: a species-level database of life history, ecology, and geography of extant and recently extinct mammals by Jones et al. (2009).

df_Traits <- read_tsv("../data/PanTHERIA.tsv")

# To get an idea of how the column names look & what traits and types of data this file contains we use names() and summary().
names(df_Traits)
summary(df_Traits)

# Checking that the data frame holds data on the Mustelidae family in the second column (MSW05_Family).
"Mustelidae" %in% unique(df_Traits[[2]]) # It does.

# Originally, this block of code (handling the PanTHERIA data) was done first because I wanted to make sure that my Family of choice (Mustelidae) was present in this data set and also get an idea of what trait data I could possibly use.


## 3b. Cleaning the Traits Data Column Names: ----

# By cleaning the column names, they become easier to work with:

# Column names as they appear in the current data frame.
Colnames_Original <- as.vector(names(df_Traits))
Colnames_Original

# Cleaning the first 5 column names (that all start with MSW05_) followed by the rest of the column names (that all start with a 1-2 digit number, followed by a -, another number, and then a _)
Colnames_Clean <- Colnames_Original %>%
  str_replace("^MSW05_", "") %>%
  str_replace("^[0-9]+-[0-9]+_", "")

# Checking that the column names have in fact been cleaned:
Colnames_Clean

# Assigning the new, cleaned column names to the df_Traits data frame.
names(df_Traits) <- Colnames_Clean
names(df_Traits)

# We can remove a few objects since we won't need them any longer.
rm(Colnames_Original, Colnames_Clean)


# This is a glimpse into the current state of the data frame df_Traits (checking that the new, cleaned column names have been assigned):
head(df_Traits)


# Since our research question is: Are body size and habitat type correlated within Mustelidae while accounting for phylogeny?, we will specifically want to look at the following trait data for Mustelidae:
# 1. Adult Body Mass
# 2. Habitat Breadth (A proxy for habitat generalism/specialization)




## 3c. Creating a Traits Data Frame With Only Required Data: ----

# We only need data for Mustelidae, so we can filter for rows where the Family column read Mustelidae. This specific data will be put into the new data frame Mustelidae_Traits.
df_Mustelidae_Traits <- df_Traits %>%
  filter(Family == "Mustelidae")

df_Mustelidae_Traits

# Checking if the new data frame only holds Mustelidae data.
unique(df_Mustelidae_Traits$Family) # Yes, it does.


# Now we want to only keep the columns for our traits of interest: Adult Body Mass & Habitat Breadth. We will also keep the Binomial (combined Genus and Species) column, as it identifies the actual organisms within the family and will be needed for downstream analysis. Additionally, we want to keep the columns for Latitude and Longitude because they will be used later in the code for the formation of an interactive map.
df_Mustelidae_Selected <- df_Mustelidae_Traits %>%
  select("Binomial", "AdultBodyMass_g", "HabitatBreadth", "GR_MidRangeLat_dd", "GR_MidRangeLong_dd")

# Check that this new data frame has only 5 columns with the names above.
head(df_Mustelidae_Selected)

# Get an idea of the data.
summary(df_Mustelidae_Selected)

# From summary(df_Mustelidae_Selected), we saw that missing values in this data are represented by -999 or -999.0. We want to replace these values with NA so that functions can be applied to the data later to remove them.
df_Mustelidae_Selected[df_Mustelidae_Selected == -999] <- NA

# Check that cells with -999 and -999.0 values have been replaced by NAs.
View(df_Mustelidae_Selected)
# Check how many NA values each column has.
summary(df_Mustelidae_Selected)

# Now we can remove rows containing NA values using drop_na() (from the tidyr package), as they may interfere with downstream analysis.
df_Mustelidae_Selected_Traits <- df_Mustelidae_Selected %>%
  drop_na(AdultBodyMass_g, HabitatBreadth)

View(df_Mustelidae_Selected_Traits)

# Double-check that there are no NA values left for Body Mass and Habitat Breadth in the data frame.
summary(df_Mustelidae_Selected_Traits)




## 4. Merging Sequence & Trait Data: ----

# These are the two sets of data that we want to merge:
# View(df_Species.Seq_Unique)
# View(df_Mustelidae_Selected_Traits)

# We want to see what species are common between the two data sets because we only want to keep those entries for analysis. The names being compared are found in df_Species.Seq_Unique$Species_Name and df_Mustelidae_Selected_Traits$Binomial.

# Merging Trait and Sequence Data with dplyr function inner.join().
df_Merged_Data <- df_Species.Seq_Unique %>%
  inner_join(df_Mustelidae_Selected_Traits, by = c("Species_Name" = "Binomial"))

# Check how many species end up being common between the two data sets.
df_Merged_Data # There are 16 species in common.

# Rearranging df_Merged_Data for better organization and readability. and removing the column Sequence_Length as we no longer need this information.
df_Merged_Data <- df_Merged_Data[, c("Title", "Species_Name", "Sequence", "AdultBodyMass_g", "HabitatBreadth", "GR_MidRangeLat_dd", "GR_MidRangeLong_dd")]

# Checking that the columns have now been rearranged and that Sequence_Length is now gone.
df_Merged_Data


## 5. Visualizing the Selected Mustelidae Traits: ----

# Trait 1: Creating a Boxplot for Adult Body Mass (g).
boxplot(df_Merged_Data$AdultBodyMass_g,
  main = "Figure 1: Boxplot of Adult Body Mass in Mustelidae",
  ylab = "Adult Body Mass (g)",
  col = "lightblue"
)

# Trait 2: Creating a Histogram for Habitat Breadth.
hist(df_Merged_Data$HabitatBreadth,
  main = "Figure 2: Histogram of Habitat Breadth in Mustelidae",
  xlab = "Habitat Breadth",
  col = "lightcoral"
)



## 6a. Aligning (MSA) the Unique Mustelidae Sequences: ----

# Converting the sequences to a DNAStringSet object so it is able to be used in the alignment algorithms & functions in packages like muscle and DECIPHER (BrowseSeqs()) (as we will be using below.)
DNAStringSet_Seqs <- DNAStringSet(df_Merged_Data$Sequence)

# Assigning species names (Species_Name) to the DNAStringSet object
names(DNAStringSet_Seqs) <- df_Merged_Data$Species_Name

# Performing a full multiple sequence alignment.
df_Alignment_Full <- DNAStringSet(
  muscle::muscle(
    DNAStringSet_Seqs,
    log = "alignment_log.txt",
    verbose = TRUE
  ),
  use.names = TRUE
)

# Viewing the full alignment.
BrowseSeqs(df_Alignment_Full)

# Saving the aligned sequences to a file.
writeXStringSet(df_Alignment_Full, "../data/Mustelidae_Full_Aligned_Sequences.fasta")



## 6b. Inspecting the Alignment: ----

# Checking the total number of sequences in the alignment and making sure they match with the number we saw in the merged data frame (16).
length(df_Alignment_Full) # Yes, it matches.

# Counting the number of gaps ("-") in each sequence of the alignment.
Gap_Counts <- lapply(
  as.character(df_Alignment_Full),
  str_count,
  "-"
)

Gap_Counts

# Calculating summary statistics of the gap counts.
summary(unlist(Gap_Counts))

# Minimum (42 gaps): The sequence with the fewest gaps has 42.
# 1st Quartile (51 gaps): 25% of sequences have 51 or fewer gaps.
# Median (82 gaps): Half of the sequences have 82 or fewer gaps.
# Mean (71.75 gaps): On average, sequences contain about 72 gaps.
# 3rd Quartile (83.25 gaps): 75% of sequences have 83 or fewer gaps.
# Maximum (111 gaps): The sequence with the most gaps has 111.

# On average, each sequence in the alignment contains about 72 gaps (insertions/deletions). Since all the sequence fall within the range of 600-700 bp in length, then having around 72 gaps means roughly 10-12% of the alignment positions are gaps. This is a relatively high percentage, especially since the species are all within the same Family.

# However, when we view df_Alignment_Full, we can see that all the gaps are at the ends of the sequences. End gaps are commonly seen in alignment and are generally not indicative of major evolutionary events (like insertions/deletions) (Schwartz et al., 2003).

# The middle region of the alignment however, contains no gaps for the 16 Mustelidae species, indicating a good quality alignment.


## 7. Constructing a Phylogenetic Tree for Mustelidae: ----

# Convert aligned sequences to DNAbin format for tree construction.
Alignment_DNAbin <- as.DNAbin(df_Alignment_Full)

# Compute a distance matrix using a chosen substitution model.
Distance_Matrix <- dist.dna(
  Alignment_DNAbin,
  model = "TN93",
  pairwise.deletion = TRUE
)

head(Distance_Matrix)

# Constructing the tree using neighbour joining.
Tree <- nj(Distance_Matrix)

# Plotting the tree.
plot(Tree, main = "Figure 4. Phylogenetic Tree of Mustelidae", cex = 0.7)



## 8. Is There a Relationship Between the Two Traits?: ----

# We can observe our two traits: Habitat Breadth and Adult Body Mass (g) on an interactive map to see if they exhibit geographic trends and whether we can observe whether there is a pattern observed between the two traits. Given that habitat breadth is count data and only ranges between 1 to 2 in our filtered data, we cannot necessarily observe a correlation between our two traits, but we can lay out the data on a map to have a better understanding of the species in our data set and the traits they exhibit. Important note: Meles meles is not available on the map as it had missing geographic data.

# Prepare the data by making Habitat Breadth a categorical variable so it can be visualised on the map
mapData <- df_Merged_Data %>%
  mutate(HabitatBreadth = as.factor(HabitatBreadth))

# Define a color palette for Adult Body Mass and Habitat Breadth using two different palettes, with Adult Body Mass treated as a continuous variable but Habitat Breadth treated as a categorical variable
adult_body_mass_palette <- colorNumeric(
  palette = "RdBu",
  domain = mapData$AdultBodyMass_g
)

habitat_palette <- colorFactor(
  palette = "viridis",
  domain = mapData$HabitatBreadth
)

# Create the interactive map using the map from OpenStreetMap through the leaflet package.The radius and the border colours of the circles on our map are indicators of Adult Body Mass whereas the colour filled inside the circle is an indicator of Habitat Breadth. Ensure species names are visible when hovered over on the interactive map. Legends are provided to provide insight.
Interactive_Map <- leaflet(mapData) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~GR_MidRangeLong_dd, lat = ~GR_MidRangeLat_dd,
    radius = ~ sqrt(AdultBodyMass_g) * 0.4,
    color = ~ adult_body_mass_palette(AdultBodyMass_g),
    opacity = 1,
    fillColor = ~ habitat_palette(HabitatBreadth),
    fillOpacity = 0.5,
    popup = ~ paste(
      "Species:", Species_Name, "<br>",
      "Adult Body Mass:", AdultBodyMass_g, "g"
    ),
    label = ~Species_Name,
    labelOptions = labelOptions(noHide = F, direction = "top", style = list("font-weight" = "bold"))
  ) %>%
  addLegend(
    position = "bottomright",
    pal = habitat_palette,
    values = mapData$HabitatBreadth,
    title = "Habitat Breadth"
  ) %>%
  addLegend(
    position = "bottomleft",
    pal = adult_body_mass_palette,
    values = mapData$AdultBodyMass_g,
    title = "Adult Body Mass"
  )

# Display the interactive map
Interactive_Map


## 9. Do These Two Traits Have a Phylogenetic Signal? ----

# Blomberg’s K and Pagel’s λ are statistical tests that can be used to measure a phylogenetic signal, which is a key aspect of answering my research question: Are body size and habitat type correlated within Mustelidae while accounting for phylogeny?


# Making sure that the species names in the tree and the trait data match by setting row names in the trait data to the species names.
row.names(df_Merged_Data) <- df_Merged_Data$Species_Name

# Creating a function to run both Blomberg’s K and Pagel’s λ tests for a given trait
run_tests <- function(trait, tree) {
  # Running Blomberg’s K test
  K_result <- phylosig(tree, trait, method = "K", test = TRUE)
  print(K_result)

  # Running Pagel’s lambda test
  lambda_result <- phylosig(tree, trait, method = "lambda", test = TRUE)
  print(lambda_result)
}

# Trait 1: Adult Body Mass
adult_body_mass <- setNames(df_Merged_Data$AdultBodyMass_g, row.names(df_Merged_Data))
cat("Blomberg's K and Pagel's lambda result for Adult Body Mass respectively are:\n")
run_tests(adult_body_mass, Tree)

# Trait 2: Habitat Breadth
habitat_breadth <- setNames(df_Merged_Data$HabitatBreadth, row.names(df_Merged_Data))
cat("Blomberg's K and Pagel's lambda result for Habitat Breadth respectively are:\n")
run_tests(habitat_breadth, Tree)




## Brief Discussion of the Results of the Tests:

# Both traits demonstrate significant phylogenetic signal, with Adult Body Mass showing slightly stronger evidence of phylogenetic influence than Habitat Breadth.
# This supports the idea that body mass and habitat preferences are inherited traits within the Mustelidae family, yet body mass appears to be more strongly conserved.
# (Discussed in more detail in my report.)
