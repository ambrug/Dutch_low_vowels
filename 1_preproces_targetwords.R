rm(list=ls())
library(tidyverse)
'%nin%' = Negate('%in%')

# set directory to output files
filedir <- paste0("/Users/annabruggeman/sciebo/Dutch_aa/Ranalysis", "/outputs_F3_3points/")

## initiate empty df
allv <- data.frame()
## append all CGN acoustic extraction output files
for (i in list.files(filedir, "*.txt")) {
  output <- read_csv(paste0(filedir, i))
  allv <- rbind(allv, output)
}

#############################################################
########## GET EXTERNAL INFO                 ################
#############################################################

##############################v
######## speaker info
##############################
##  sex1 = male; sex2 = female; sexX = unknown
speak <- read_tsv("/Users/annabruggeman/Documents/CGN_data/CGN_2.0.3/data/meta/text/speakers.txt")
speak[speak=="sex1"] <- "m"
speak[speak=="sex2"] <- "f"
speak[speak=="sexX"] <- "x"

#select only relevant columns
speak2 <- speak %>%
  select(ID, sex, birthYear, eduPlace, language, birthRegion, resRegion, eduLevel)
## merge speaker info with aa
allv <- merge(x = allv, y = speak2, by.x = "speaker", by.y = "ID", all.x = TRUE)

##############################v
######## WORD FREQ INFO
##############################
subtlex <- read_tsv("/Users/annabruggeman/sciebo/Dutch_aa/Ranalysis/SUBTLEX-NL.txt")
subtlex$Word <- as.factor(subtlex$Word)
# select only useful columns to be merged
subtlexmerge <- subtlex %>%
  select(Word, SUBTLEXWF, Lg10WF, SUBTLEXCD, Lg10CD)

## remove punctuation at end of words; ? and .
# original words contain punct
sum(grepl(".", allv$word, fixed = TRUE))
# remove them
allv$cleanword <- gsub('[.,!?]', '', allv$word)
# no more
sum(grepl(".", allv$cleanword, fixed = TRUE))

# this merge takes a long time
allvfreq <- merge(x = allv, y = subtlexmerge, by.x = "cleanword", by.y = "Word", all.x = TRUE)


##############################v
######## polysyllabic target words
##############################
targets_poly <- read.csv("/Users/annabruggeman/sciebo/Dutch_aa/Ranalysis/wordlist_a_poly.csv", sep=";")

# take only target words with target stress patterns
tp <- targets_poly %>% 
  filter(cat_v1 %nin% "N/A")

levels(as.factor(tp$cat_v1))

## check number of polysyllabic target words
nrow(tp)
tp %>%
  group_by(cat_v1,stress.pattern) %>%
  summarise(total=n())

## MERGE allvfreq with the polysyllabic info
allvfreqtp <- merge(x = allvfreq, y = tp, by.x = "cleanword", by.y = "word", all.x = TRUE)

## crosscheck CGN label (seg) with phonological label (cat_v1) - not identical
allvfreqtp %>%
  group_by(cat_v1, stress.pattern, seg) %>%
  summarise(count=n())

##############################v
######## monosyllabic target words
##############################
# get list of monosyllabic target words
targets <- read.table("/Users/annabruggeman/sciebo/Dutch_aa/Ranalysis/wordlist_a_mono_lemalph.txt", skipNul=T, header=T)
targets$stress.pattern <- "S"
  
# get number of monosyllabic target words
length(unique(targets$word)) 
length(unique(tolower(targets$word)))

targets %>%
  mutate(lowerword = as.factor((tolower(word)))) %>%
  distinct(lowerword, .keep_all = T) %>%
  group_by(cat_v1) %>%
  summarise(countn = n())
## 517 a
## 324 aa

## if allfreqtp dataframe "cleanword" matches a word in targets
## - cat_v1 should be taken from target dataframe
allvfreqtptm <- merge(x = allvfreqtp, y = targets, by.x ="cleanword", by.y = "word", all.x = TRUE)
# merge cat_v1.x and cat_v1_y
# merge stress.pattern.x and stress.pattern.y
allvfreqtptm <- allvfreqtptm %>%
  unite(col="stress.pattern",c("stress.pattern.x", "stress.pattern.y"), na.rm=T) %>%
  unite(col="cat_v1",c("cat_v1.x", "cat_v1.y"), na.rm=T)
# set number of syllables to 1
allvfreqtptm$nsyll[allvfreqtptm$cleanword %in% targets$word] <- 1

# check
colnames(allvfreqtptm)
levels(as.factor(allvfreqtptm$stress.pattern))
levels(as.factor(allvfreqtptm$cat_v1))

#################################################################
########## CORRECT SEG INFO ########################################
#################################################################
#CGN label seg a should only be a1-stressed or a1-unstressed - but seems like CGN has annotation issues there a/A is used a bit inconsistently
#create new column with phonological info
allvfreqtptm$phonolseg <- substr(allvfreqtptm$cat_v1, 1,2)
allvfreqtptm[allvfreqtptm$phonolseg == "a1", ]$phonolseg <- "A"
allvfreqtptm[allvfreqtptm$phonolseg == "aa", ]$phonolseg <- "a"

levels(as.factor(allvfreqtptm$phonolseg))

### SAVE
save(allvfreqtptm, file = "allvfreqtptm.RData")

#################################################################
########## MONOSYLLABLES ########################################
#################################################################
## select only rows with A/aa
aafreq <- allvfreqtptm %>%
  filter(seg == "a" | seg== "A")
## 2.1 million

# create vector with all mono targetwords
targetvect_lower <- unique(tolower(targets$word)) #841 items
targetvect <- targets$word #934 items

## select only rows with targetwords
## 1.29million + based on all folders
monaafreq <- aafreq %>%
  filter(cleanword %in% targetvect)
length(unique(monaafreq$cleanword[monaafreq$cleanword %in% targetvect])) ##830 of the target items occur in the extracted data

### SAVE
save(monaafreq, file = "monaa.RData")

############## INFOS ABOUT ACTUAL CGN WORD OCCURRENCES
length(unique(monaafreq$cleanword)) #889 unique words (incl capitals)
length(unique(tolower(monaafreq$cleanword))) #816 (excl capitalisation diffs)

# get number of monosyllabic words
monowords_targetoverview <- monaafreq %>%
  group_by(cleanword) %>%
  summarise(count=n()) %>%
  arrange(count)

sum(is.na(monaafreq$SUBTLEXWF)) #6826 vowels have no matching word freq
sum(is.na(monaafreq$SUBTLEXWF))/nrow(monaafreq) # but this is about 0.5 per cent of total vowels

unique(monaafreq[which(is.na(monaafreq$SUBTLEXWF)),]$cleanword) ## lots of capitalised words dont have Freq info.

#################################################################
########## POLYSYLLABLES ########################################
#################################################################
## aafreq is dataframe with vowel a/aa

## select only rows with CGN code a/A and a matching targetword
tpaafreq <- aafreq %>%
  filter(cleanword %in% tp$word)

# get number of polysyllabic words
polywords_targetoverview <- tpaafreq %>%
  group_by(cleanword) %>%
  summarise(count=n()) %>%
  arrange(count)

### consider excluding words that have schwa in them (cat_v2 and cat_v3)? This reduces it by nearly half. Dont do it here yet.
tpaafreq_subs <- tpaafreq %>%
  filter(cat_v2 != "@") %>%
  filter(cat_v3 != "@")

### SAVE
save(tpaafreq, file = "tpaa.RData")
