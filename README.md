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
Scyphozoa <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Scyphozoa&format=tsv")

#write the data into the working directory
write_tsv(Scyphozoa, "Scyphozoa_BOLD_data.tsv")

#read the data from the file into R and store in data object Scyph
Scyph <- read_tsv("Scyphozoa_BOLD_data.tsv")

#review variables in Scyph
names(Scyph)

#remove Scyphozoa
rm(Scyphozoa)

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

#calculate the thymine proportion and add to the end of the data object.
ScyphBayofSton.NucFreq <- ScyphBayofSton.NucFreq %>%
  mutate(Tproportion = ((T) / (A + T + G + C)))

#calculate the mean of the Tproportion (0.347, meaning 34.7% of the sequence is thymine)
mean(ScyphBayofSton.NucFreq$Tproportion)


#With careful manipulation of the BOLD data on jellyfish (schypho), it can be determined that the global distribution of jellyfish is broken down into three main areas. One of these distributions includes Croatia, where an abundance of lethal jellyfish have been reported. To dive into the unique regions of Croatia where jellyfish are found, we ran a spantree for the regions and found that the Bay of Ston had the highest dissimilarity. We then dived deeper into the species found in this region, Aurelia sp. 7-MND-2005, and determined certain genetic traits about their nucleotide sequences that may lead the way to understanding these species, and hopefully determining whether they are the lethal and poisonous species we are trying so hard to understand.
