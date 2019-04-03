#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.0.1
# Description: generate table with available GPSeq tracks.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

l = lapply(c("here", "argparser"), FUN = function(x) {
	if ( !suppressMessages(require(package = x, character.only = T)) )
		install.packages(pkgs = x)
	suppressMessages(library(package = x, character.only = T))
})

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('
', name = 'mkRecapTable.R')

# Define elective arguments
parser = add_argument(parser, arg = '--input-dir', short = '-i', type = class(""),
	help = 'Root folder (centrality_by_seq). Default: script folder',
	default = here(), nargs = 1)
parser = add_argument(parser, arg = '--cell-meta', short = '-m', type = class(""),
	help = 'Path to tsv table with seqrun|cell_type|synch_status', nargs = 1)

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

all_scores = c("prob_2p", "prob_f", "prob_g", "cor_2p", "cor_f", "cor_g", "roc_2p",
	"roc_f", "roc_g", "var_2p", "var_f", "ff_2p", "ff_f", "cv_2p", "cv_f",
	"svar_2p", "svar_f", "svar_g", "scv_2p", "scv_f", "scv_g", "sff_2p",
	"sff_f", "sff_g")

# FUNCTIONS ====================================================================

parse_root_folder = function(folder_label, meta = NA) {
	outl = list()

	if ( folder_label %in% list.files(input_dir) ) {
		curdir = file.path(input_dir, folder_label)

		# Parse sequencer run folders
		for ( seqrun in list.files(curdir) ) {
			seqrun_dir = file.path(curdir, seqrun)

			# Parser condition folders
			for ( conds in list.files(seqrun_dir) ) {
				conds_dir = file.path(seqrun_dir, conds)

				# Parse rank files
				rankfiles = list.files(conds_dir)
				rankfiles = rankfiles[grepl("settings", rankfiles)]

				if ( 0 != length(rankfiles) ) {
					outl[[length(outl)+1]] = do.call(rbind, lapply(rankfiles,
						FUN = function(x) parse_settings_file(file.path(conds_dir, x), meta)))
				} else {
					print(paste0("Cannot find any ranks in: ", conds_dir))
				}
			}
		}
	}

	return(do.call(rbind, outl))
}

extract_value = function(seed, settings) {
	output = settings[grepl(seed, settings)]
	if ( 0 != length(output) ) {
		output = unlist(strsplit(output, " "))
		output = output[length(output)]
	} else {
		output = NA
	}
	return(output)
}

extract_bool = function(seed, match, settings) {
	out = settings[grepl(seed, settings)]
	out = grepl(match, out)
	return(out)
}

parse_settings_file = function(file_path, meta = NA) {
	c = scan(file = file_path, what = "character", sep = "~", quiet = T)

	version = c[grepl("^ # GPSeq-centrality-estimate ", c)]
	version = unlist(strsplit(version, " "))
	version = version[length(version)]

	version_steps = as.numeric(unlist(strsplit(
		substr(version, 2, nchar(version)), ".", fixed = T)))
	if ( version_steps[1] <= 3 & version_steps[2] < 5 ) {
		out = parse_settings_file_before_3_5_0(file_path, meta)
	} else {
		out = parse_settings_file_after_3_5_0(file_path, meta)
	}

	return(out)
}

parse_settings_file_before_3_5_0 = function(file_path, meta = NA) {
	c = scan(file = file_path, what = "character", sep = "~", quiet = T)
	
	r_path = unlist(strsplit(basename(file_path), "settings"))
	r_path = paste0(r_path[1], "ranked", r_path[2])
	r_path = unlist(strsplit(r_path, '.', fixed = T))
	r_path = r_path[-length(r_path)]
	r_path = paste0(paste(r_path, collapse = "."), ".tsv")
	if ( !file.exists(file.path(dirname(file_path), r_path)) ) {
		cat(paste0("ERROR: missing estimates file for '", file_path, "'. Skipped!\n"))
		return(NULL)
	}

	version = c[grepl("^ # GPSeq-centrality-estimate ", c)]
	version = unlist(strsplit(version, " "))
	version = version[length(version)]

	date = c[grepl("^@", c)]
	date = unlist(strsplit(substr(date, 2, nchar(date)), " "))[1]

	group_size = extract_value("Group size", c)
	bin_size = extract_value("Bin size", c)
	bin_step = extract_value("Bin step", c)

	custom_bins = NA
	treatment = NA
	flag_fields = unlist(strsplit(unlist(strsplit(basename(file_path), "settings"))[1], ".", fixed = T))
	if ( 2 == sum(is.na(c(bin_size, bin_step))) ) {
		if ( 0 == sum(grepl("Bin bed", c)) ) {
			custom_bins = "chr-wide"
			
			if ( 2 > length(flag_fields) ) {
				treatment = NA
			} else {
				treatment = flag_fields[2]
			}
		} else {
			if ( 2 > length(flag_fields) ) {
				cat(paste0("ERROR: missing custom_bin prefix in '", file_path, "'. Skipped!\n"))
				return(NULL)
			}
			custom_bins = flag_fields[2]

			if ( 3 > length(flag_fields) ) {
				treatment = NA
			} else {
				treatment = flag_fields[3]
			}
		}
	} else {
		if ( 2 > length(flag_fields) ) {
			treatment = NA
		} else {
			treatment = flag_fields[2]
		}
	}

	outliers = c[grepl("Outliers", c)]
	outliers = unlist(strsplit(outliers, " "))
	outliers = outliers[length(outliers)-2]

	alphalim = extract_value("Alpha", c)
	outliers_rmAll = extract_bool("Outliers", "remove all", c)

	normalized = 1 == sum(grepl("Normalizing", c))

	masked_input = extract_value("Input mask", c)
	masked_output = extract_value("Output mask", c)

	scores = ""
	inscore = grepl("Included metrics", c)
	if ( 0 != sum(inscore) ) scores = gsub(",", ";", gsub(" ", "", c[which(inscore) + 1]))
	exscore = grepl("Excluded metrics", c)
	if ( 0 != sum(exscore) ) {
		scores = gsub(" ", "", c[which(exscore) + 1])
		scores = unlist(strsplit(scores, ","))
		scores = paste(all_scores[all_scores %in% scores], collapse = ";")
	}
	if ( 0 == nchar(scores) ) scores = paste(all_scores, collapse = ";")

	seq_run = basename(dirname(dirname(file_path)))
	cell_line = NA
	synch_status = NA
	if ( class(NA) != class(meta) ) {
		if ( seq_run %in% meta$seqrun ) {
			cell_line = meta$cell_line[meta$seqrun == seq_run]
			synch_status = meta$synch_status[meta$seqrun == seq_run]
		} else {
			if ( grepl("+", seq_run) ) {
				tmp = unlist(strsplit(seq_run, "+", fixed = T))[1]
				if ( tmp %in% meta$seqrun ) {
					cell_line = meta$cell_line[meta$seqrun == tmp]
					synch_status = meta$synch_status[meta$seqrun == tmp]
				} else {
					cat(paste0("WARNING: no cell metadata found for '", dirname(file_path), "'\n"))
				}
			}
		}
	}

	beds = gsub(" ", "", sub("\\([0-9]+\\)", "", c[(which(grepl("Bed files", c))+1):length(c)]))
	beds = unlist(lapply(beds, FUN = basename))

	libraries = paste(unlist(lapply(beds, FUN = function(x) {
		gsub("-", "", unlist(strsplit(unlist(strsplit(x, "(TK|JC)"))[2], "_"))[1])
	})), collapse = ";")
	conditions = paste(unlist(lapply(beds, FUN = function(x) {
		gsub("-", "", unlist(strsplit(unlist(strsplit(x, "(TK|JC)"))[2], "_"))[2])
	})), collapse = ";")

	out = data.frame(
		seq_run = seq_run,
		#file = basename(file_path),
		libraries = libraries,
		conditions = conditions,
		cell_line = cell_line,
		synch_status = synch_status,
		flag = treatment,
		date = date,
		#user = NA,
		version = version,
		custom_bins = custom_bins,
		bin_size = bin_size,
		bin_step = bin_step,
		group_size = group_size,
		bed_outliers = outliers,
		bed_alphalim = alphalim,
		bed_outliers_rmAll = outliers_rmAll,
		score_outliers = NA,
		score_alphalim = NA,
		normalized = normalized,
		masked_input = masked_input,
		masked_output = masked_output,
		estimates = scores,
		notes = NA,
		stringsAsFactors = F
	)

	return(out)
}

parse_settings_file_after_3_5_0 = function(file_path, meta = NA) {
	c = scan(file = file_path, what = "character", sep = "~", quiet = T)
	
	r_path = unlist(strsplit(basename(file_path), "settings"))
	r_path = paste0(r_path[1], "rescaled", r_path[2])
	r_path = unlist(strsplit(r_path, '.', fixed = T))
	r_path = r_path[-length(r_path)]
	r_path = paste0(paste(r_path, collapse = "."), ".tsv")
	if ( !file.exists(file.path(dirname(file_path), r_path)) ) {
		cat(paste0("ERROR: missing estimates file for '", file_path, "'. Skipped!\n"))
		return(NULL)
	}

	version = c[grepl("^ # GPSeq-centrality-estimate ", c)]
	version = unlist(strsplit(version, " "))
	version = version[length(version)]

	date = c[grepl("^@", c)]
	date = unlist(strsplit(substr(date, 2, nchar(date)), " "))[1]

	group_size = extract_value("Group size", c)
	bin_size = extract_value("Bin size", c)
	bin_step = extract_value("Bin step", c)

	custom_bins = NA
	treatment = NA
	flag_fields = unlist(strsplit(unlist(strsplit(basename(file_path), "settings"))[1], ".", fixed = T))
	if ( 2 == sum(is.na(c(bin_size, bin_step))) ) {
		if ( 0 == sum(grepl("Bin bed", c)) ) {
			custom_bins = "chr-wide"
			
			if ( 2 > length(flag_fields) ) {
				treatment = NA
			} else {
				treatment = flag_fields[2]
			}
		} else {
			if ( 2 > length(flag_fields) ) {
				cat(paste0("ERROR: missing custom_bin prefix in '", file_path, "'. Skipped!\n"))
				return(NULL)
			}
			custom_bins = flag_fields[2]

			if ( 3 > length(flag_fields) ) {
				treatment = NA
			} else {
				treatment = flag_fields[3]
			}
		}
	} else {
		if ( 2 > length(flag_fields) ) {
			treatment = NA
		} else {
			treatment = flag_fields[2]
		}
	}

	bed_outliers_id = which(grepl("Bed outliers", c))
	bed_outliers = unlist(strsplit(c[bed_outliers_id + 1], " "))
	bed_outliers = bed_outliers[length(bed_outliers)-2]
	bed_alphalim = unlist(strsplit(c[bed_outliers_id + 2], " "))
	bed_alphalim = bed_alphalim[length(bed_alphalim)]
	bed_outliers_rmAll = grepl("remove all", c[bed_outliers_id + 1])

	score_outliers_id = which(grepl("Score outliers", c))
	score_outliers = unlist(strsplit(c[score_outliers_id + 1], " "))
	if ( grepl("remove all", c[score_outliers_id + 1]) ) {
		score_outliers = score_outliers[length(score_outliers) - 2]
	} else {
		score_outliers = score_outliers[length(score_outliers)]
	}
	score_alphalim = unlist(strsplit(c[score_outliers_id + 2], " "))
	score_alphalim = score_alphalim[length(score_alphalim)]

	normalized = 1 == sum(grepl("Normalizing", c))

	masked_input = extract_value("Input mask", c)
	masked_output = extract_value("Output mask", c)

	scores = ""
	inscore = grepl("Included metrics", c)
	if ( 0 != sum(inscore) ) scores = gsub(",", ";", gsub(" ", "", c[which(inscore) + 1]))
	exscore = grepl("Excluded metrics", c)
	if ( 0 != sum(exscore) ) {
		scores = gsub(" ", "", c[which(exscore) + 1])
		scores = unlist(strsplit(scores, ","))
		scores = paste(all_scores[all_scores %in% scores], collapse = ";")
	}
	if ( 0 == nchar(scores) ) scores = paste(all_scores, collapse = ";")

	seq_run = basename(dirname(dirname(file_path)))
	cell_line = NA
	synch_status = NA
	if ( class(NA) != class(meta) ) {
		if ( seq_run %in% meta$seqrun ) {
			cell_line = meta$cell_line[meta$seqrun == seq_run]
			synch_status = meta$synch_status[meta$seqrun == seq_run]
		} else {
			if ( grepl("+", seq_run) ) {
				tmp = unlist(strsplit(seq_run, "+", fixed = T))[1]
				if ( tmp %in% meta$seqrun ) {
					cell_line = meta$cell_line[meta$seqrun == tmp]
					synch_status = meta$synch_status[meta$seqrun == tmp]
				} else {
					cat(paste0("WARNING: no cell metadata found for '", dirname(file_path), "'\n"))
				}
			}
		}
	}

	beds = gsub(" ", "", sub("\\([0-9]+\\)", "", c[(which(grepl("Bed files", c))+1):length(c)]))
	beds = unlist(lapply(beds, FUN = basename))

	libraries = paste(unlist(lapply(beds, FUN = function(x) {
		gsub("-", "", unlist(strsplit(unlist(strsplit(x, "(TK|JC)"))[2], "_"))[1])
	})), collapse = ";")
	conditions = paste(unlist(lapply(beds, FUN = function(x) {
		gsub("-", "", unlist(strsplit(unlist(strsplit(x, "(TK|JC)"))[2], "_"))[2])
	})), collapse = ";")

	out = data.frame(
		seq_run = seq_run,
		#file = basename(file_path),
		libraries = libraries,
		conditions = conditions,
		cell_line = cell_line,
		synch_status = synch_status,
		flag = treatment,
		date = date,
		#user = NA,
		version = version,
		custom_bins = custom_bins,
		bin_size = bin_size,
		bin_step = bin_step,
		group_size = group_size,
		bed_outliers = bed_outliers,
		bed_alphalim = bed_alphalim,
		bed_outliers_rmAll = bed_outliers_rmAll,
		score_outliers = score_outliers,
		score_alphalim = score_alphalim,
		normalized = normalized,
		masked_input = masked_input,
		masked_output = masked_output,
		estimates = scores,
		notes = NA,
		stringsAsFactors = F
	)

	return(out)
}

# RUN ==========================================================================

if ( !is.na(cell_meta) )
	cell_meta = read.delim(cell_meta, as.is = T, header = T)

t1 = parse_root_folder("FISH", cell_meta)
t2 = parse_root_folder("Seq", cell_meta)

opath = paste0(input_dir, "/GPSeq_centrality_track.tsv")
write.table(rbind(t1, t2), opath,
	quote = F, row.names = F, col.names = T, sep = "\t")

cat(paste0("Output written to: ", opath, "\n"))

# END --------------------------------------------------------------------------

################################################################################
