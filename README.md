# AbdulJamaliAssignment1PartC
Assignment 3

#I, Abdul Jamali, hereby declare that this is my original work and no one else has worked on it.

#The following is based off of a television show I watched over the summer, I could not find a credible source for this however there are many google search results indicating that croatia is home to many lethal jellyfish species.

#Assignment 1 part C:
#Determine the global distribution of Jellyfish and gain insight on Croatia's unique jellyfish species, which have said to be extremely deadly in only one region of the country. Reports have been made that Croatia has a vast array of jellyfish species, tourists and locals report deaths from contact to poisonous jellyfish located in one specific region of Croatia. How can we use R to find the most unique species make up of Croatia's jellyfish, and how can we understand them more? That is the goal of this code.

#install all necessary packages to work with
install.packages("tidyverse")
install.packages("vegan")
install.packages("ape")
install.packages("iNEXT")
install.packages("RSQLite")

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("muscle")
biocLite("DECIPHER")

#call all libraries that are necessary to work with
library("tidyverse")
library("vegan")
library("Biostrings")
library("ape")
library("muscle")
library("DECIPHER")
library("iNEXT")


#obtain data for jellyfish by using the bold api tool
Scyph <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Scyphozoa&format=tsv")

#write the data into the working directory
write_tsv(Scyph, "Scyphozoa_BOLD_data.tsv")

#review variables in Scyph
names(Scyph)


######
## JB Edit #1----------
# your initial data import using the bold api tool loaded the dataframe and saved it as an object (Scyphozoa).  You saved a local copy using write_tsv, but then read the object back into the workspacee and assigned it to a new object (Scyph).  If you use the identical() function, it shows that the two objects Scyphozoa and Scyph are actually the same dataframe.  You did it right the first time by saving the BOLD tsv to a variable, but then theres no need to re-write the data as a new object.  Also then you won't need to remove the duplicate object later in the script. 
######


#remove all data with NA in the BIN listings
Scyph1 <- Scyph[(grep(":", Scyph$bin_uri)),]

#create a new data object with only bins, bin frequencies and their respective countries
Scyph3 <- Scyph1 %>%
  group_by(country, bin_uri) %>%
  count(bin_uri)

#using the tidyverse feature spread, spread out the data amongst collums instead of rows
Scyph4 <- spread(Scyph3, bin_uri, n)

#replace all "NA" data with numeric "0"
Scyph4[is.na(Scyph4)] <- 0

#remove rownames as vegan can only manipulate numeric data
Scyph4 <- Scyph4 %>%
  remove_rownames %>%
  column_to_rownames(var="country")

#create a species accumulation curve for a global perspective of jellyfish discovery per sites
Scyph.accum <- specaccum(Scyph4)

#plot the species acuumulation curve
plot(Scyph.accum)

#Using vegdist and spantree, create a spantree plot of the dissimilarities in species variation by country (Croatia is close to multiple other countries such as Italy, in this spantree... this indicates that these lethal jellyfish may be found in other countries as well.)
vegScyph <- vegdist(Scyph4)
treeScyph <- spantree(Scyph4)
plot(treeScyph, cmdscale(vegScyph), type="t") 

#dive deeper and find only records from Croatia
ScyphCroatia <- Scyph[(grep("Croatia", Scyph$country)),]

ScyphCroatia1 <- ScyphCroatia %>%
  group_by(region) %>%
  summarize(count = length(processid)) %>%
  arrange(desc(count))

ScyphCroatia2 <- ScyphCroatia %>%
  group_by(region, bin_uri) %>%
  count(bin_uri)

ScyphCroatia3 <- spread(ScyphCroatia2, bin_uri, n)

ScyphCroatia3[is.na(ScyphCroatia3)] <- 0

#remove rownames for the regian and prepare for vegan to manipulate data
ScyphCroatia3 <- ScyphCroatia3 %>%
  remove_rownames %>%
  column_to_rownames(var="region")

#create another spantree plot to determine the most unique region of jellyfish distributio (Bay of Stron)
vegScyphCroat <- vegdist(ScyphCroatia3)
treeScyphCroat <- spantree(ScyphCroatia3)
plot(treeScyphCroat, cmdscale(vegScyphCroat), type="t")

#Only look for records from Bay of Ston
ScyphBayofSton <- Scyph[(grep("Bay of Ston", Scyph$region)),]

#determine the species of Jellyfish found in the Bay of Ston region (Aurelia sp. 7-MND-2005)
unique(ScyphBayofSton$species_name)

#convert the nucleotide sequences into data objects that can be manipulated by biostrings
ScyphBayofSton$nucleotides <- DNAStringSet(ScyphBayofSton$nucleotides)

#determine the nucleotide frequencies and store into a data object (it is apparent that there is almost double the amount of thymine than other nucleotides in the sequences... this may be indicative of something special)
ScyphBayofSton.NucFreq <- as.data.frame(letterFrequency(ScyphBayofSton$nucleotides, letters = c("A", "C", "G", "T")))

######
# JB Edit #2------
## I think this dataframe with the nucleotide counts and T-Frequency is the most interesting part of the analysis.  The problem is that when I view this dataframe I'm not entirely sure what I'm looking at.  I'll add more information to this dataframe and store it in a new object (ScyphBayofSton.NucFreq_jbEdit).  I'll start by using rownames() to add the Process ID identifiers:
ScyphBayofSton.NucFreq_jbEdit <- ScyphBayofSton.NucFreq
rownames(ScyphBayofSton.NucFreq_jbEdit) <- ScyphBayofSton$processid

## I'll also add the bin_uri, the region, and the sequence itself:
ScyphBayofSton.NucFreq_jbEdit$bin_uri <- ScyphBayofSton$bin_uri
ScyphBayofSton.NucFreq_jbEdit$Region <- ScyphBayofSton$region
ScyphBayofSton.NucFreq_jbEdit$COI_Sequence <- ScyphBayofSton$nucleotides
## Re-arrange columns
ScyphBayofSton.NucFreq_jbEdit <- ScyphBayofSton.NucFreq_jbEdit[,c(7,6,1,2,3,4,5,8)]
######


#calculate the thymine proportion and add to the end of the data object.
ScyphBayofSton.NucFreq <- ScyphBayofSton.NucFreq %>%
  mutate(Tproportion = ((T) / (A + T + G + C)))

#calculate the mean of the Tproportion (0.347, meaning 34.7% of the sequence is thymine)
mean(ScyphBayofSton.NucFreq$Tproportion)


######
# JB Edit #3-------------
## using the spantree to find a discrepency between species of the Bay of Ston location and the other 8 regions is interesting.  I'm going to try to expand on it.  You did a quick nucleotide analysis by computing the A,C,T, and G count and the T-Frequencies for the Bay of Ston region. I'm going to take this a step further by computing A,C,T,G counts and T-frequencies for two additional regions (South Adriatic and Veliko Jezero,Mljet), then see how they compare to the Bay of Ston data.

## South Adriatic Region
## Filter for Records in South Adriatic Region
Scyph_s.adriatic <- Scyph %>%
  filter(region == "South Adriatic") %>%
  print()
## Convert to DNAbin
Scyph_s.adriatic$nucleotides <- DNAStringSet(Scyph_s.adriatic$nucleotides)
## Calculate Nucleotide Counts and T Frequency
Scyph_s.adriatic.NucFreq <- as.data.frame(letterFrequency(Scyph_s.adriatic$nucleotides, letters = c("A", "C", "G", "T"))) 
Scyph_s.adriatic.NucFreq <- Scyph_s.adriatic.NucFreq %>% mutate(Tproportion = ((T) / (A + T + G + C)))
# Set processid as Rownames
rownames(Scyph_s.adriatic.NucFreq) <- Scyph_s.adriatic$processid

## Veliko Jezero,Miljet Region
## Filter for Records in VJM Region
Scyph_vjm <- Scyph %>%
  filter(region == "Veliko Jezero,Mljet") %>%
  print()
## Convert to DNAbin
Scyph_vjm$nucleotides <- DNAStringSet(Scyph_vjm$nucleotides)
## Calculate Nucleotide Counts and T Frequency
Scyph_vjm.NucFreq <- as.data.frame(letterFrequency(Scyph_vjm$nucleotides, letters = c("A", "C", "G", "T")))
Scyph_vjm.NucFreq <- Scyph_vjm.NucFreq %>% mutate(Tproportion = ((T) / (A + T + G + C)))
# Set processid as Rownames
rownames(Scyph_vjm.NucFreq) <- Scyph_vjm$processid

## we can then compare these 3 dataframes:
ScyphBayofSton.NucFreq
Scyph_s.adriatic.NucFreq
Scyph_vjm.NucFreq
######

#With careful manipulation of the BOLD data on jellyfish (schypho), it can be determined that the global distribution of jellyfish is broken down into three main areas. One of these distributions includes Croatia, where an abundance of lethal jellyfish have been reported. To dive into the unique regions of Croatia where jellyfish are found, we ran a spantree for the regions and found that the Bay of Ston had the highest dissimilarity. We then dived deeper into the species found in this region, Aurelia sp. 7-MND-2005, and determined certain genetic traits about their nucleotide sequences that may lead the way to understanding these species, and hopefully determining whether they are the lethal and poisonous species we are trying so hard to understand.
