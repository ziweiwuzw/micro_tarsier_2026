#!/usr/bin/env Rscript
# Full-depth KMA profile and count-level rarefaction; base R only.
args <- commandArgs(TRUE)
if (length(args) != 6L) stop("Usage: Rscript profile_rarefaction.R COUNTS CONTIG_MAP GE PAIRED_READ_BASES SAMPLE_ID OUT_DIR")
counts_file <- args[1]; map_file <- args[2]
ge <- if (args[3] == "NA") NA_real_ else as.numeric(args[3])
paired_bases <- as.numeric(args[4]); sid <- args[5]; out <- args[6]
if (!is.na(ge) && (!is.finite(ge) || ge <= 0)) stop("Genome equivalents must be positive")
if (!is.finite(paired_bases) || paired_bases <= 0) stop("Paired read bases must be positive")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

mapping <- read.delim(map_file, check.names=FALSE, stringsAsFactors=FALSE)
observed <- read.delim(counts_file, check.names=FALSE, stringsAsFactors=FALSE)
if (!all(c("contig_id", "sgb_id", "length") %in% names(mapping))) stop("Invalid contig map")
if (!all(c("contig_id", "fragment_count") %in% names(observed))) stop("Invalid count table")
if (anyDuplicated(mapping$contig_id) || anyDuplicated(observed$contig_id)) stop("Duplicate contig ID")
if (any(!is.finite(mapping$length)) || any(mapping$length <= 0)) stop("Invalid contig length")
if (any(!is.finite(observed$fragment_count)) ||
    any(observed$fragment_count < 0) ||
    any(observed$fragment_count != floor(observed$fragment_count))) stop("Invalid integer fragment count")

position <- match(observed$contig_id, mapping$contig_id)
if (anyNA(position)) stop("KMA contig missing from reference map")
count <- numeric(nrow(mapping)); count[position] <- observed$fragment_count
sgb <- unique(mapping$sgb_id)
group <- match(mapping$sgb_id, sgb)
genome_length <- numeric(length(sgb))
genome_length[as.integer(rownames(rowsum(mapping$length, group, reorder=FALSE)))] <-
    as.numeric(rowsum(mapping$length, group, reorder=FALSE))
if (any(genome_length <= 0)) stop("Invalid SGB genome length")

aggregate_count <- function(x, groups) {
    answer <- numeric(length(sgb))
    if (!length(x)) return(answer)
    summed <- rowsum(x, groups, reorder=FALSE)
    answer[as.integer(rownames(summed))] <- as.numeric(summed)
    answer
}
profile <- function(n) {
    per_kb <- n / (genome_length / 1000)
    total <- sum(per_kb)
    tpm <- if (total > 0) 1e6 * per_kb / total else numeric(length(sgb))
    data.frame(
        sample_id=sid, sgb_id=sgb, length_bp=genome_length,
        fragment_count=n,
        fpkb_per_ge=if (is.na(ge)) rep(NA_real_, length(sgb)) else per_kb / ge,
        coverage_per_ge=if (is.na(ge)) rep(NA_real_, length(sgb)) else
            n * paired_bases / genome_length / ge,
        tpm=tpm, relative_abundance=tpm / 1e6
    )
}

full <- profile(aggregate_count(count, group))
write.table(full, file.path(out, "sgb_full.tsv"), sep="\t", quote=FALSE,
            row.names=FALSE)

depths <- as.numeric(strsplit(Sys.getenv("RARE_DEPTHS", "100000,250000,500000,1000000,2500000,5000000,10000000"), ",")[[1]])
reps <- as.integer(Sys.getenv("RARE_REPS", "5"))
seed <- as.integer(Sys.getenv("RARE_SEED", "20261008"))
if (any(!is.finite(depths)) || any(depths <= 0) ||
    any(depths != floor(depths)) || !is.finite(reps) || reps < 1L ||
    !is.finite(seed)) stop("Invalid rarefaction settings")

# Sequential hypergeometric draws reproduce sampling without replacement
# while avoiding a potentially enormous per-read vector.
draw_counts <- function(pop, size) {
    answer <- numeric(length(pop)); remaining <- sum(pop); draws <- size
    for (j in seq_along(pop)) {
        if (draws == 0) break
        answer[j] <- if (draws == remaining) pop[j] else
            rhyper(1, pop[j], remaining - pop[j], draws)
        draws <- draws - answer[j]; remaining <- remaining - pop[j]
    }
    if (sum(answer) != size) stop("Rarefaction count mismatch")
    answer
}

positive <- which(count > 0)
pool <- count[positive]
groups <- group[positive]
total <- sum(pool)
seed_offset <- sum(utf8ToInt(sid) * seq_along(utf8ToInt(sid)))
rarefied <- file(file.path(out, "rarefied_sgb.tsv"), "w")
writeLines("sample_id\tdepth\treplicate\tsgb_id\tfragment_count\ttpm", rarefied)
status <- list()
for (depth in depths) {
    for (replicate in seq_len(reps)) {
        feasible <- depth <= total
        detected <- NA_integer_
        if (feasible) {
            set.seed(as.integer((seed + seed_offset + depth + replicate * 104729) %%
                                .Machine$integer.max))
            sampled <- draw_counts(pool, depth)
            by_sgb <- aggregate_count(sampled, groups)
            p <- profile(by_sgb)
            hits <- which(p$fragment_count > 0)
            detected <- length(hits)
            if (length(hits)) {
                write.table(
                    data.frame(sample_id=sid, depth=depth, replicate=replicate,
                               sgb_id=p$sgb_id[hits], fragment_count=p$fragment_count[hits],
                               tpm=p$tpm[hits]),
                    rarefied, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
                )
            }
        }
        status[[length(status) + 1L]] <- data.frame(
            sample_id=sid, depth=depth, replicate=replicate,
            status=if (feasible) "OK" else "INSUFFICIENT_TEMPLATE_ASSIGNMENTS",
            detected_sgbs=detected
        )
    }
}
close(rarefied)
write.table(do.call(rbind, status), file.path(out, "rarefaction_status.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
