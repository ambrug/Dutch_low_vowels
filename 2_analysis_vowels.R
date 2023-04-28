rm(list=ls())
library(tidyverse)
library(phonR)
library(lme4)
library(emmeans)
library(MuMIn)

'%nin%' = Negate('%in%')

outpath = "/Users/annabruggeman/sciebo/Talks-Pubs/2023-08 ICPhS Dutcha/"

load("monaa.Rdata") # 1.42 million target a monosylls
load("tpaa.Rdata") # 76k target a polysylls
load("allvfreqtptm.Rdata") #10.38 million items. ALL EXTRACTED CORPUS VOWELS.

# create a new column with the F1 and F2 measured at midpoint for all the monophtongal, and at 25% for /eː, øː, oː/
allvfreqtptm$f1mix <- 99999
allvfreqtptm$f2mix <- 99999
allvfreqtptm$seg <- as.factor(allvfreqtptm$seg)
levels(allvfreqtptm$seg)

allvfreqtptm <- allvfreqtptm %>%
mutate(f1mix = case_when(seg %in% c("e", "o", "2") ~ f1_25, 
                         seg %nin% c("e", "o", "2") ~ f1_50),
       f2mix = case_when(seg %in% c("e", "o", "2") ~ f2_25, 
                         seg %nin% c("e", "o", "2") ~ f2_50))


########################################################################
#### LOB NORMALISATION + PREPROCES FOR VOWEL SPACE
########################################################################
oldnofrow <- nrow(allvfreqtptm)

### NB this Lobanov-scaling might have some off-labels, based on the quality of CGN transcription. 
## words like Baas sometimes has A and sometimes a etc.
lob <- allvfreqtptm %>%
  # filter removing missing variables for VL/NL
  filter(!is.na(vlnl)) %>%
  group_by(speaker) %>%
  #remove cases with following /r/ or /l/
  filter(nextseg %nin% c("r", "l")) %>%
  mutate(f1lob = as.numeric(normLobanov(f1mix)), ## same as scale
         f2lob = as.numeric(normLobanov(f2mix))) 
## result: over 8 million tokens

newnofrow <- nrow(lob)
## (newnofrow - oldnofrow)/oldnofrow --> minus 20%

## get distribution per folder
lob %>%
  group_by(folder) %>%
  summarise(count=n())

## REMOVE OUTLIERS per vowel category
lob <- lob %>%
  ungroup() %>%
  group_by(seg) %>%
  mutate(f1_scaledvowel = scale(f1mix),
         f2_scaledvowel = scale(f2mix)) %>%
  filter(between(f1_scaledvowel, -3, 3)) %>%
  filter(between(f2_scaledvowel, -3, 3))

post3sd <- nrow(lob)
## (post3sd - newnofrow)/newnofrow --> minus 3 more per cent

## add stress info for those words that have it
lob$stress <- c(as.character(NA))
lob[grepl("-unstressed", lob$cat_v1),]$stress <- "unstressed"
lob[grepl("-stressed", lob$cat_v1),]$stress <- "stressed"
lob$stress <- as.factor(lob$stress)                        

# add poly/mono for words that have it
lob$monopoly <- c(as.character(NA))
lob[nchar(lob$stress.pattern) > 1,]$monopoly <- "poly"
lob[nchar(lob$stress.pattern) == 1,]$monopoly <- "mono"
lob$monopoly <- as.factor(lob$monopoly)                     

#### MEANS per vowel
lobmean <- lob %>%
  # first get speaker/vowel means
  group_by(seg, speaker, vlnl) %>%
  summarise(meanf1=mean(as.numeric(f1lob)),
            meanf2=mean(as.numeric(f1lob)),
            count=n()) %>%
  # checked with various filters: if taking means only for segs that have 10 or more tokens per speaker it still doesnt affect the overall vowel space means
  # then get means per VlNl group
  group_by(seg, vlnl) %>%
  summarise(allmeanf1=mean(meanf1, na.rm=T),
            allmeanf2=mean(meanf2, na.rm=T)) 

## get numbers 
sumslob <- lob %>%
  group_by(vlnl) %>%
  summarise(countspeaker = length(unique(speaker)),
            countitem = n())

# 5.05 mil NL, 2480 speakers
# 3.03 mil, 1254 speakers

## modify names of factor levels
lob$vlnl <- as.factor(lob$vlnl)
lob$seg <- as.factor(lob$seg)
levels(lob$seg)
levels(lob$vlnl) <- c("Netherlands", "Belgium")

########################################################################
#### PLOT VOWEL SPACE
########################################################################

ggplot(data=lob,
       aes(x = f2lob, y = f1lob, 
           color = seg, fill=seg)) +
  geom_text(data=lob
            %>% group_by(seg, vlnl) %>% 
              summarise_at(vars(f1lob:f2lob), mean, na.rm = TRUE), aes(label = seg), size = 9) +
  facet_wrap(~vlnl) +
  scale_colour_manual(values=c("darkgrey", "chartreuse2", 
                               "darkgoldenrod2", "darkgoldenrod2",
                               "red", "red",
                               "cyan3", "cyan3",
                               "blue3", "blue3",
                               "deeppink",
                               "darkmagenta", "darkmagenta")) +
  scale_x_reverse(name="Mean Lobanov F2") + 
  scale_y_reverse(name="Mean Lobanov F1") +
  theme_bw() +
  theme(legend.position="none",
        axis.title=element_text(size=15),
        axis.text=element_text(size=12),
        strip.background.x=element_rect(fill="white"),
        strip.text.x=element_text(size=20))
  
ggsave(paste0(outpath, "vowelspace_lob.eps"), width = 8, height = 5)


#######################################################
############ FULL DATA A/AA: 4X STRESS LEVELS ##########
#######################################################
lob$stress.pattern <- as.factor(lob$stress.pattern)
lob$phonolseg <- as.factor(lob$phonolseg)
levels(lob$stress.pattern)

lobaa <- droplevels(lob %>% 
                      # filter all corpus As
                      filter(seg %in% c("a", "A"))) %>%
                      # filter for only the ones that were manually annotated (phonolseg) for initial target position
                      filter(phonolseg %in% c("a", "A")) %>%
                      # exclude "wwwSw", "wwwS" so that unstressed is only primary-stress adjacent
                      filter(stress.pattern %nin% c("wwwSw", "wwwS")) %>%
                      mutate(stresstype = case_when(stress.pattern %in% c("wS", "wSw", "wSww", "wSwww") ~ "unstressed",
                                                    stress.pattern %in% c("Sw", "Sww") ~ "primary (poly)",
                                                    stress.pattern %in% c("S") ~ "primary (mono)", 
                                                    stress.pattern %in% c("wwS", "wwSw", "wwSww") ~ "secondary"))

## get numbers of target low vowels
## by NL/belgium
lobaa %>%
  group_by(vlnl) %>%
  summarise(count=n())
# by vowel
lobaa %>%
  group_by(phonolseg) %>%
  summarise(count=n())
# by stress pattern
lobaa %>%
  group_by(seg, stress.pattern) %>%
  summarise(count=n())

lobaa$stresstype <- as.factor(lobaa$stresstype)

## get schwa
lob_schwa <- droplevels(lob %>% 
                          filter(seg %in% c("@")))
lob_schwa$phonolseg <- "schwa"
lob_schwa$stresstype <- "unstressed"

################ GET MEANS PER SPEAKER FOR MODELS AND OVERALL MEANS FOR PLOT

#### SCHWA ####
################
## per speaker, check. filter for 4 or more tokens per speaker
lob_schwa_speakers <- lob_schwa %>%
  group_by(vlnl,phonolseg,stresstype,speaker,eduPlace) %>%
  summarise(meanf1=mean(f1lob),
            meanf2=mean(f2lob),
            count=n()) %>% 
  filter(count > 3)

# across the board
lob_schwa_mean <- lob_schwa_speakers %>%
  group_by(vlnl) %>%
  summarise(meanf1all= mean(meanf1),
            meanf2all= mean(meanf2))

#### A/AA ####
################
### per speaker, filter for 4 or more tokens per speaker
### ALL LM data below uses this filter!! Not mentioned in paper because no space.
lobaa_multstress_speakers <- lobaa %>%
  group_by(phonolseg,speaker,vlnl,stresstype,eduPlace) %>%
  summarise(meanf1=mean(f1lob),
            meanf2=mean(f2lob),
            count=n()) %>% #14558 rows (ie speaker means)
  filter(count > 3) #9630

# check total numbers of speakers per variety and a/aa
lobaa_multstress_speakers %>%
  group_by(phonolseg, vlnl) %>%
  summarise(speakers=length(unique(speaker)))


# across the board 
lobaa_multstress_mean <- lobaa_multstress_speakers %>%
  group_by(phonolseg,vlnl,stresstype) %>%
  ungroup(count,speaker) %>%
  summarise(meanf1all=mean(meanf1),
            meanf2all=mean(meanf2))

# reorder with unstressed left and stressed right
lobaa_multstress_mean$stresstype <- factor(lobaa_multstress_mean$stresstype, 
                                           levels=c('unstressed',
                                                    'primary (poly)',
                                                    'secondary',
                                                    'primary (mono)'))

#########################################################################
### GRAPH ALL TOGETHER means, 4x stress, a/aa/@ #########################
#########################################################################
ggplot() +
  geom_text(data=lob_schwa_mean,
           aes(x=meanf2all, y=meanf1all, label = "\U0259"), size = 12, colour="darkgrey") +
  geom_point(data=lobaa_multstress_mean,
             aes(x = meanf2all, y = meanf1all, 
                 shape = stresstype, fill=stresstype), size=7) +
  geom_line(data=lobaa_multstress_mean,
            aes(x = meanf2all, y = meanf1all, 
                group=phonolseg, linetype=phonolseg)) +
  scale_linetype_manual(values=c(1,2), guide="none") +
  scale_shape_manual(name="",
                     values=c(21,12,25,22)) +
  scale_fill_manual(name="",
                    values=c("white", "black","lightgrey","black")) +
  scale_colour_manual(name="",
                      values=c("white", "black","lightgrey", "black")) +
  facet_wrap(~vlnl, nrow=2) +
  annotate(geom="text", x=-0.83, y=0.5, label="\U0251", size=10) +
  annotate(geom="text", x=-0.43, y=0.8, label="a", size=10) +
  scale_x_reverse("Mean Lobanov F2") + 
  scale_y_reverse("Mean Lobanov F1") +
  theme_bw() +
  theme(legend.text=element_text(size=22),
        legend.position="top",
        strip.text.x = element_text(size = 26),
        strip.background.x=element_rect(fill="white")) +
  guides(colour=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         shape=guide_legend(nrow=2, byrow=TRUE))
ggsave(paste0(outpath, "a_aa_4xstress_4ormore.png"), width = 5, height = 10)


#######################################################
############ ANTWERPEN DATA A/AA: 4X STRESS LEVELS ##########
#######################################################

## get the Antwerpen A's (some NL speakers but stick with codes)
anta <- lobaa %>%
  group_by(phonolseg,speaker,vlnl,stresstype,eduPlace,birthYear) %>%
  summarise(meanf1=mean(f1lob),
            meanf2=mean(f2lob),
            count=n()) %>%
  filter(eduPlace %in% c(paste0("B-", 200:206)))

anta$eduPlace <- as.factor(anta$eduPlace)
anta$birthYear <- as.numeric(anta$birthYear)

## get the speaker means for schwa
lob_schwa_speakers_ant <- lob_schwa_speakers %>%
  filter(eduPlace %in% c(paste0("B-", 200:206)))

## get the overall means
lob_schwa_antmean <- lob_schwa_speakers_ant %>%
  ungroup() %>%
  summarise(meanf1all=mean(meanf1),
            meanf2all=mean(meanf2))

## merge A/AA with @ 
antall <- rbind(anta, lob_schwa_speakers_ant)

## by year, vowel, stress
antall_mean <- antall %>%
  group_by(phonolseg,stresstype) %>%
  summarise(meanf1all=mean(meanf1),
            meanf2all=mean(meanf2),
            count=n()) %>%
  filter(phonolseg %in% c("A", "a"))

antall_mean$stresstype <- as.factor(antall_mean$stresstype)
antall_mean$stresstype <- ordered(antall_mean$stresstype, 
                                  levels=c("unstressed", "primary (poly)",
                                           "secondary", "primary (mono)"))

## turn into factors
antall$speaker <- as.factor(antall$speaker)
antall$stresstype <- as.factor(antall$stresstype)
antall$eduPlace <- as.factor(antall$eduPlace)

# nof speakers
unique(antall$speaker) # 55
 #216 a/aa tokens
antall %>%
  filter(phonolseg != "schwa") %>%
  nrow()

######################################################################
### GRAPH ANTWERPEN means, 4x stress, a/aa/@ #########################
######################################################################
ggplot() +
  geom_text(data=lob_schwa_antmean,
            aes(x=meanf2all, y=meanf1all, label = "\U0259"), size = 12, colour="darkgrey") +
  geom_point(data=antall_mean,
             aes(x = meanf2all, y = meanf1all, 
                 shape = stresstype, fill=stresstype), size=7) +
  geom_line(data=antall_mean,
            aes(x = meanf2all, y = meanf1all, 
                group=phonolseg, linetype=phonolseg)) +
  scale_linetype_manual(values=c(1,2), guide="none") +
  scale_shape_manual(name="",
                     values=c(21,12,25,22)) +
  scale_fill_manual(name="",
                    values=c("white", "black","lightgrey","black")) +
  scale_colour_manual(name="",
                      values=c("white", "black","lightgrey", "black")) +
  annotate(geom="text", x=-0.79, y=0.5, label="\U0251", size=10) +
  annotate(geom="text", x=-0.38, y=0.4, label="a", size=10) +
  scale_x_reverse("Mean Lobanov F2") + 
  scale_y_reverse("Mean Lobanov F1") +
  theme_bw() +
  theme(legend.text=element_text(size=22),
        legend.position="top",
        strip.text.x = element_text(size = 26),
        strip.background.x=element_rect(fill="white")) +
  guides(colour=guide_legend(nrow=2, byrow=TRUE),
         fill=guide_legend(nrow=2, byrow=TRUE),
         shape=guide_legend(nrow=2, byrow=TRUE))
ggsave(paste0(outpath, "a_aa_4xstress_ant.png"), width = 5, height = 5)



######################################
## STATS ########################
######################################
# first get lobaa and schwa together
loba_schwa_speakers <- rbind(lob_schwa_speakers, lobaa_multstress_speakers)
loba_schwa_speakers$phonolseg <- as.factor(loba_schwa_speakers$phonolseg)
loba_schwa_speakers$stresstype <- as.factor(loba_schwa_speakers$stresstype)

################
## LINEAR MODELS
################

###########################################################################
######## MEAN MODEL WITH 1 MEAN PER SPEAKER/VOWEL (only cases of >3): R2 up to 0.75
###########################################################################
## subsets with means per speaker/vowel (a/aa + schwa)
loba2 <- droplevels(loba_schwa_speakers %>%
                      filter(phonolseg %in% c("A", "schwa")))
lobaa2 <- loba_schwa_speakers %>%
  filter(phonolseg %in% c("a", "schwa"))

######## SPLIT INTO VARIETIES
loba_schwa_speakers_nl <- loba_schwa_speakers %>%
  filter(vlnl %in% "Netherlands")
loba_schwa_speakers_vl <- loba_schwa_speakers %>%
  filter(vlnl %in% "Belgium")

#######################
## NETHERLANDS: a versus aa versus schwa * stress
#######################
f1_nl <- lmer(data=loba_schwa_speakers_nl,
              meanf1 ~ phonolseg * stresstype + 
                (1|speaker),
              #+ (0+speaker|phonolseg), #makes R terminate
              REML=F)
f1_nl_0 <- lmer(data=loba_schwa_speakers_nl,
                meanf1 ~ phonolseg + stresstype + 
                  (1|speaker), REML=F)
# interaction stress/segment, duh. x2(3)=1439.2, p<0.0001
anova(f1_nl, f1_nl_0)
summary(f1_nl)

f2_nl <- lmer(data=loba_schwa_speakers_nl,
              meanf2 ~ phonolseg * stresstype + 
                (1|speaker), REML=F)
f2_nl_0 <- lmer(data=loba_schwa_speakers_nl,
                meanf2 ~ phonolseg + stresstype + 
                  (1|speaker), REML=F)
# interaction here too x2(3)=74.9, p<0.0001
anova(f2_nl, f2_nl_0)

############# MULT COMP
## F1: most are pairwise different, save the ones that are NOT diff
nodiff_f1 <- as_tibble(emmeans(f1_nl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value > 0.05)
# a/aa very similar only in this case:
# a secondary=A unstressed est=0.05, SE=0.04

# info about schwa
diff_f1 <- as_tibble(emmeans(f1_nl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value < 0.05)
# a approaches schwa but is still different:
# a unstressed != schwa in terms of F1: est=0.12, SE=0.01, p=0.0000

## F2
nodiff_f2 <- as_tibble(emmeans(f2_nl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value > 0.05)
# main diffs:
# a secondary=A secondary (est=0.10, se=0.04)/A unstressed (est=0.07,se=0.03)

diff_f2 <- as_tibble(emmeans(f2_nl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value < 0.05)
# a unstressed != schwa in terms of F1: est=-0.25, SE=0.01, p=0.0000

### EXPLAINED VARIANCE NL
r.squaredGLMM(f1_nl) #62/70%
r.squaredGLMM(f2_nl) #55/71%

#######################
## BELGIUM: a versus aa versus schwa * stress
#######################
f1_vl <- lmer(data=loba_schwa_speakers_vl,
              meanf1 ~ phonolseg * stresstype + 
                (1|speaker), REML=F)
f2_vl <- lmer(data=loba_schwa_speakers_vl,
              meanf2 ~ phonolseg * stresstype + 
                (1|speaker), REML=F)
summary(f1_vl)
emmeans(f1_vl, ~ phonolseg+stresstype)

## F1: most are pairwise different
nodiff_f1_vl <- as_tibble(emmeans(f1_vl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value > 0.05)
diff_f1_vl <- as_tibble(emmeans(f1_vl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value < 0.05)
# a secondary = A primary (mono)/A secondary/A unstressed

## a secondary = A secondary/A unstressed/A primary
## F2
nodiff_f2_vl <- as_tibble(emmeans(f2_vl, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value > 0.05)
# a secondary= A primary poly/A secondary/A unstressed 

### EXPLAINED VARIANCE VL
r.squaredGLMM(f1_vl) #66/76%
r.squaredGLMM(f2_vl) #54/74%

######################################################################
############## JUST ANTWERPEN ##########################################

## MODELS ANTWERPEN: a versus aa (not schwa) * stress
f1_ant <- lmer(data=antall %>%
                 filter(phonolseg %in% c("A", "a")),
               meanf1 ~ phonolseg * stresstype +
                 (1|speaker), REML=F)
f2_ant <- lmer(data=antall %>%
                 filter(phonolseg %in% c("A", "a")),
               meanf2 ~ phonolseg * stresstype +
                 (1|speaker), REML=F)

## F1: a few comps across tense/lax are different
antf1_diff <- as_tibble(emmeans(f1_ant, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value < 0.05)
## 1 a primary (mono) - A primary (mono)    0.462 0.0872  171.    5.30 0.00000970   
## 3 a primary (mono) - A unstressed        0.549 0.158   181.    3.47 0.0148       
## 4 A primary (mono) - a primary (poly)   -0.530 0.112   180.   -4.73 0.000123     
## 6 a primary (poly) - A unstressed        0.617 0.167   175.    3.69 0.00720      
## 7 A primary (poly) - a unstressed        0.599 0.155   179.    3.85 0.00398 
## summary:
# a prim (mono) != A prim (mono), A unstressed
# a prim (poly) != A prim (poly), A unstressed
# a unstressed ! = A prim (poly)

# F2: 1 comp across tense/lax is different
antf2_diff <- as_tibble(emmeans(f2_ant, pairwise ~ phonolseg+stresstype)$contrasts) %>%
  filter(p.value < 0.05)
## 1 A primary (mono) - a unstressed   -0.318 0.0691  177.   -4.60 0.000213

############################################################
###################### RANDOM FORESTS ######################
############################################################

### create RF-specific dataset
lobaa_for_rf <- lobaa %>%
  select(phonolseg,seg,speaker,vlnl,folder,seg_dur,Lg10WF,resRegion,eduLevel,birthRegion,birthYear,stresstype,stress.pattern,eduPlace,cleanword,f1lob,f2lob)
lobaa_for_rf$birthYear <- as.numeric(as.character(lobaa_for_rf$birthYear))
lobaa_for_rf$seg_dur <- as.numeric(as.character(lobaa_for_rf$seg_dur))
lobaa_for_rf$resRegion <- as.factor(as.character(lobaa_for_rf$resRegion))
lobaa_for_rf$eduLevel <- as.factor(as.character(lobaa_for_rf$eduLevel))
lobaa_for_rf$birthRegion <- as.factor(as.character(lobaa_for_rf$birthRegion))
lobaa_for_rf$speaker <- as.factor(as.character(lobaa_for_rf$speaker))

## only tense a and means
lobaa_rf <- lobaa_for_rf %>%
  group_by(speaker, stresstype, phonolseg, vlnl, folder, Lg10WF, resRegion, eduLevel, birthRegion, birthYear, eduPlace, cleanword) %>%
  summarise(meanf1=mean(f1lob),
            meanf2=mean(f2lob),
            count=n()) %>% # 116k
  filter(count > 3,
         phonolseg == "a") #13k

#create dataset with no NAs in predictors
lobaa_rf_nona <- lobaa_rf %>%
  filter_at(vars(meanf1,meanf2,eduPlace,birthYear,birthRegion,resRegion,
                 eduLevel,Lg10WF), all_vars(!is.na(.))) #12k

sqrt(11)
RFpreds_all <- 'speaker + vlnl + 
                        folder +
                        stresstype +
                        birthYear + birthRegion +
                        resRegion + 
                        eduLevel + eduPlace +
                        cleanword + Lg10WF'
myFormula_f2 <- as.formula(paste0('meanf2 ~ ', RFpreds_all))

library(ranger)
##### on TENSE A ONLY
myforestf2 <- ranger(myFormula_f2,
                     data = lobaa_rf_nona,
                     importance = 'permutation', mtry = 3, num.trees = 2000)
r2 <- round(myforestf2$r.squared,2)

## tidy up variable importance
varimpf2 <- as.data.frame(sort(myforestf2$variable.importance))
varimpf2$predictor <- rownames(varimpf2)
colnames(varimpf2)[1] <- "importance"
varimpf2$predictor <- factor(varimpf2$predictor)
varimpf2$predictor <- fct_reorder(varimpf2$predictor, varimpf2$importance, .desc = F)
levels(varimpf2$predictor)[levels(varimpf2$predictor) %in% "cleanword"] <- "word"
levels(varimpf2$predictor)[levels(varimpf2$predictor) %in% "vlnl"] <- "Belg/Neth"
levels(varimpf2$predictor)[levels(varimpf2$predictor) %in% "stresstype"] <- "stress"
levels(varimpf2$predictor)[levels(varimpf2$predictor) %in% "folder"] <- "CGN folder"

### PLOT
ggplot(varimpf2) +
  geom_point(aes(y=predictor, x=importance), size=4) +
  scale_x_continuous("Relative variable importance") +
  annotate("text", label=paste0("R2=",r2), x = 0.053, y = 1.35, size=7) +
  theme_bw() +
  theme(axis.text.y =element_text(size=20),
        axis.title.x =element_text(size=15),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank())
ggsave(paste0(outpath, "varimp.eps"), width = 5, height = 4)


