## Anna Bruggeman, anna.bruggeman@uni-bielefeld.de
## last mod 05.12.2022
## Extracts acoustic infos from CGN wav files + TextGrids with corpus provided automatic word/segment annotation

# set base path for WAV and TGRID
cgnmainwav$ = "/Users/annabruggeman/Documents/CGN_data/CGN_2.0.3/data/audio/wav/"
cgnmaintgrid$ = "/Users/annabruggeman/Documents/CGN_data/CGN_2.0.3/data/annot/text/awd/"

## create lists for subdirectories
Create Strings as directory list: "dir_list", cgnmainwav$
nof_dirs = Get number of strings

Create Strings as tokens: "nl vl"
Rename: "nlvllist"

## list for the different speakers/tiers
Create Strings as tokens: "1 4 7"
Rename: "tierlist"

# loop through all folders
# 16=comp-o
# first folder a=1 is checksums, run from 2 onwards
for a from 16 to 16
	selectObject: "Strings dir_list"
	dir$ = Get string: 'a'
	appendInfoLine: "'dir$'"

	if left$(dir$,4) == "comp"

		# loop through both subfolders, nl/ and vl/
		for b from 1 to 2
			selectObject: "Strings nlvllist"
			vlnl$ = Get string: 'b'

			# delete existing output file
			deleteFile: "/Users/annabruggeman/sciebo/Dutch_aa/Ranalysis/output_'dir$'_'vlnl$'.txt"
			# create new output file
			outfile$ = "/Users/annabruggeman/sciebo/Dutch_aa/Ranalysis/output_'dir$'_'vlnl$'.txt"
			# write header of output file
			writeFileLine: "'outfile$'", "folder,file,vlnl,speaker,interval,word,word_start,word_dur,
					 ...seg,seg_dur,f1_25,f1_50,f1_75,f2_25,f2_50,f2_75,f3_25,f3_50,f3_75,nextseg,preseg"

			# set path to tgrids and create list
			tgridpath$ = "'cgnmaintgrid$''dir$'/'vlnl$'/"
			Create Strings as file list: "tgridlist", "'tgridpath$'*.TextGrid"
			nof_tgrids = Get number of strings
			appendInfoLine: "has 'nof_tgrids' tgrids"

			# set path to wavs and create list
			wavpath$ = "'cgnmainwav$''dir$'/'vlnl$'/"
			Create Strings as file list: "wavlist", "'wavpath$'*.wav"
			nof_wav = Get number of strings

			#if nof_tgrids <> nof_wav
			#	appendInfoLine: "not equal wav and tgrid for 'dir$'"
			#else
			#	appendInfoLine: "same no of wavs and tgrids for 'dir$'"
			#endif

			#########################
			# GO THROUGH TGRIDS
			#########################
			for i from 1 to nof_tgrids

				selectObject ("Strings tgridlist")
				current_file$ = Get string: 'i'
				name_prefix$ = current_file$ - ".TextGrid"

				# access TextGrid
				Read from file: "'tgridpath$''current_file$'"
				# access WAV
				Read from file: "'wavpath$''name_prefix$'.wav"

				# tgrid operations
				selectObject: "TextGrid 'name_prefix$'"

				#there is always a first speaker (tier 1-3)
				speaker1$ = Get tier name: 1
				nof_int3 = Get number of intervals: 3

				# check if there is a second/third speaker
				ntiers = Get number of tiers
				if ntiers == 3
					nofpasses = 1
				elsif ntiers == 6
					nofpasses = 2
					nof_int6 = Get number of intervals: 6
					speaker2$ = Get tier name: 4
				elsif ntiers == 9
					nofpasses = 3
					nof_int6 = Get number of intervals: 6
					speaker2$ = Get tier name: 4
					nof_int9 = Get number of intervals: 9
					speaker3$ = Get tier name: 7
				# only do the first three anyway if there are more than 3 diff speakers
				else
					nofpasses = 3
					nof_int6 = Get number of intervals: 6
					speaker2$ = Get tier name: 4
					nof_int9 = Get number of intervals: 9
					speaker3$ = Get tier name: 7
				endif

				wordcount = 0

				## extract vowels from both speakers, loop through tier 3 then optionally tier 6 and 9
				for c from 1 to nofpasses
				    selectObject: "Strings tierlist"
    				tierno$ = Get string: 'c'
					selectObject: "TextGrid 'name_prefix$'"
    				if tierno$ == "1"
						targetwordtier = 1
						targetsegtier = 3
    					nof_int = 'nof_int3'
						speaker$ = speaker1$
    				elsif tierno$ == "4"
						targetwordtier = 4
						targetsegtier = 6
    					nof_int = 'nof_int6'
						speaker$ = speaker2$
					elsif tierno$ == "7"
						targetwordtier = 7
						targetsegtier = 9
    					nof_int = 'nof_int9'
						speaker$ = speaker3$
    				endif
					appendInfoLine: "'name_prefix$' speaker: 'speaker$' on tier 'tierno$' has 'nof_int'"

                	#nof_filled_int = 0
					# there are some issues with final end of segment boundaries missing, throws an error. Ignore file-final 'a's with no end boundary.
                	for j from 1 to nof_int
						selectObject: "TextGrid 'name_prefix$'"
                		seg_label$ = Get label of interval: targetsegtier, j

                    				if seg_label$ == "a" or seg_label$ == "A" or seg_label$ == "i" or seg_label$ == "I" or seg_label$ == "E" or seg_label$ == "e" or seg_label$ == "o" or seg_label$ == "O" or seg_label$ == "y" or seg_label$ == "Y" or seg_label$ == "2" or seg_label$ == "u" or seg_label$ == "@"
                    					seg_start = Get start time of interval: targetsegtier, j
                    					seg_end = Get end time of interval: targetsegtier, j
                    					seg_dur = seg_end - seg_start
                    					seg_mid = seg_start + (seg_dur / 2)
                    					seg_1quart = seg_start + (seg_dur / 4)
                    					seg_3quart = seg_mid + (seg_dur / 4)

                    					select Sound 'name_prefix$'
                    					# Get formantx3 info at discrete timepoints
                    					Extract part: seg_start-0.05, seg_end+0.05, "rectangular", 1, "yes"
                    					To Formant (burg)... 0 5 5500 0.025 50

                    					f1_25 = Get value at time: 1, 'seg_1quart', "hertz", "linear"
                    					f1_50 = Get value at time: 1, 'seg_mid', "hertz", "linear"
                    					f1_75 = Get value at time: 1, 'seg_3quart', "hertz", "linear"

                    					f2_25 = Get value at time: 2, 'seg_1quart', "hertz", "linear"
                    					f2_50 = Get value at time: 2, 'seg_mid', "hertz", "linear"
                    					f2_75 = Get value at time: 2, 'seg_3quart', "hertz", "linear"

										f3_25 = Get value at time: 3, 'seg_1quart', "hertz", "linear"
                    					f3_50 = Get value at time: 3, 'seg_mid', "hertz", "linear"
                    					f3_75 = Get value at time: 3, 'seg_3quart', "hertz", "linear"

                    					selectObject: "Formant 'name_prefix$'_part"
                    					Remove
        								selectObject: "Sound 'name_prefix$'_part"
                    					Remove

                    					selectObject: "TextGrid 'name_prefix$'"
										if j < nof_int
                    						nextseg$ = Get label of interval: targetsegtier, j+1
										else
											nextseg$ = ""
										endif
										if j > 1
                    						preseg$ = Get label of interval: targetsegtier, j-1
										else
											preseg$ = ""
										endif

										# get matching word
										wordint = Get low interval at time: targetwordtier, seg_start+0.005
										word$ = Get label of interval: targetwordtier, wordint
                    					word_start = Get start time of interval: targetwordtier, wordint
                    					word_end = Get end time of interval: targetwordtier, wordint
                    					word_dur = word_end - word_start

                    					appendFile: "'outfile$'", "'dir$','name_prefix$','vlnl$','speaker$','j','word$','word_start:2','word_dur:3',
                    					...'seg_label$','seg_dur:3','f1_25:0','f1_50:0','f1_75:0','f2_25:0','f2_50:0','f2_75:0','f3_25:0','f3_50:0','f3_75:0','nextseg$','preseg$''newline$'"
   									endif

    				#for j from 1 to nof_int1
    				endfor

				# for c from 1 to nof_passes (speakers)
				endfor
				selectObject: "Sound 'name_prefix$'"
				plusObject: "TextGrid 'name_prefix$'"
				Remove
              		
              # for i to nof_tgrids
              endfor

		# for b from 1 to 2 (folder vl and nl)
		endfor

	# if left$(dir$)		
	endif
endfor

select all
Remove

appendInfoLine: "all acoustics extracted!"