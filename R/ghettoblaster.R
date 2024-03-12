#' @title NCBI BLAST web client
#'
#' @description This package provides a simple interface to the NCBI BLAST web service.
#'
#' @docType package
#' @name ghettoblaster
#' @import dplyr httr2 stringr glue brio worrms future furrr purrr tidyr
#' @importFrom xml2 read_xml as_list
NULL

#' Submit a sequence to NCBI BLAST
#'
#' @param sequence A DNA sequence
#' @param program The BLAST program to use
#' @param max_num_seq The maximum number of sequences to return
#' @return A request ID
#' @export
submit_sequence <- function(sequence, program = "blastn", max_num_seq = 100) {
  body <- list(
    QUERY = sequence,
    db = "nucleotide",
    GENETIC_CODE = "1",
    JOB_TITLE = "Nucleotide Sequence",
    ADV_VIEW = "on",
    stype = "nucleotide",
    DBTYPE = "gc",
    DATABASE = "nt",
    NUM_ORG = "1",
    BLAST_PROGRAMS = "megaBlast",
    MAX_NUM_SEQ = max_num_seq,
    SHORT_QUERY_ADJUST = "on",
    EXPECT = "0.05",
    WORD_SIZE = "28",
    HSP_RANGE_MAX = "0",
    MATRIX_NAME = "PAM30",
    MATCH_SCORES = "1,-2",
    GAPCOSTS = "0 0",
    COMPOSITION_BASED_STATISTICS = "0",
    FILTER = "L",
    REPEATS = "repeat_9606",
    FILTER = "m",
    TEMPLATE_LENGTH = "0",
    TEMPLATE_TYPE = "0",
    SHOW_OVERVIEW = "on",
    SHOW_LINKOUT = "on",
    GET_SEQUENCE = "on",
    FORMAT_OBJECT = "Alignment",
    FORMAT_TYPE = "HTML",
    ALIGNMENT_VIEW = "Pairwise",
    MASK_CHAR = "2",
    MASK_COLOR = "1",
    DESCRIPTIONS = "100",
    ALIGNMENTS = "100",
    LINE_LENGTH = "60",
    NEW_VIEW = "false",
    NCBI_GI = "false",
    SHOW_CDS_FEATURE = "false",
    NUM_OVERVIEW = "100",
    QUERY_INDEX = "0",
    FORMAT_NUM_ORG = "1",
    CONFIG_DESCR = "ClustMemNbr,ClustComn,Ds,Sc,Ms,Ts,Cov,Eval,Idnt,AccLen,Acc",
    CLIENT = "web",
    SERVICE = "plain",
    CMD = "request",
    PAGE = "MegaBlast",
    PROGRAM = program,
    MEGABLAST = "on",
    DEFAULT_PROG = "megaBlast",
    DB_DISPLAY_NAME = "nt",
    SELECTED_PROG_TYPE = "megaBlast",
    SAVED_SEARCH = "true",
    NUM_DIFFS = "0",
    NUM_OPTS_DIFFS = "0",
    PAGE_TYPE = "BlastSearch",
    USER_DEFAULT_PROG_TYPE = "megaBlast",
    USER_DEFAULT_MATCH_SCORES = "0"
  )
  req <- request("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi")
  res <- do.call("req_body_form", c(list(.req = req), body)) %>%
    req_perform()
  response_body <- res %>%
    resp_body_string()
  rid <- str_match(response_body, "<input name=\"RID\" type=\"hidden\" value=\"([A-Z0-9]+)\"")[,2]
}

#' Wait until a BLAST request is ready
#'
#' @param rid A request ID
#' @param verbose Whether to print status messages
#' @return TRUE if the request is ready, FALSE in case of an unexpected status
#' @export
wait_until_ready <- function(rid, verbose = FALSE) {
  url <- glue("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}")
  while (TRUE) {
    res <- request(url) %>%
      req_perform() %>%
      resp_body_string()
    status <- str_match(res, "Status=([A-Z]+)")[,2]
    if (verbose) {
      message(status)
    }
    if (status == "READY") {
      return(TRUE)
    } else if (status == "WAITING") {
      Sys.sleep(10)
    } else {
      return(FALSE)
    }
  }
}

#' Get the results of a BLAST request
#'
#' @param rid A request ID
#' @return A data frame with the results
#' @export
get_results <- function(rid) {
  url <- glue("https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=XML2&RID={rid}")
  tf <- tempfile()
  request(url) %>%
    req_perform() %>%
    resp_body_raw() %>%
    brio::write_file_raw(path = tf)
  x <- as_list(read_xml(unzip(tf, glue("{rid}_1.xml"))))
  x$BlastOutput2$report$Report$results$Results$search$Search$hits %>%
    map(function(hit) {
      list(
        len = hit$len,
        num = hit$num,
        id = hit$description$HitDescr$id,
        accession = hit$description$HitDescr$accession,
        title = hit$description$HitDescr$title,
        taxid = hit$description$HitDescr$taxid,
        sciname = hit$description$HitDescr$sciname,
        bit_score = hit$hsps$Hsp$`bit_score`,
        score = hit$hsps$Hsp$score,
        evalue = hit$hsps$Hsp$evalue,
        identity = hit$hsps$Hsp$identity,
        query_from = hit$hsps$Hsp$`query-from`,
        query_to = hit$hsps$Hsp$`query-to`,
        hit_from = hit$hsps$Hsp$`hit-from`,
        hit_to = hit$hsps$Hsp$`hit-to`,
        align_len = hit$hsps$Hsp$`align-len`,
        qseq = hit$hsps$Hsp$qseq,
        hseq = hit$hsps$Hsp$hseq
      )
    }) %>%
    bind_rows() %>%
    unnest(cols = c(len, num, id, accession, title, taxid, sciname, score, evalue, identity, query_from, query_to, hit_from, hit_to, align_len, qseq, hseq))
}

#' Get WoRMS taxonomy for a list of species names
#'
#' @param species_names A character vector of species names
#' @return A data frame with the taxonomy, input names in column input
worms_for_names <- purrr::possibly(function(species_names) {
  wm_records_names(species_names, marine_only = FALSE) %>%
    setNames(species_names) %>%
    bind_rows(.id = "input")
}, otherwise = NULL)

#' Get WoRMS taxonomy for a list of species names
#'
#' @param species_names A character vector of species names
#' @return A data frame with the taxonomy, input names in column input
#' @export
get_taxonomy <- function(species_names) {
  name_batches <- split(species_names, as.integer((seq_along(species_names) - 1) / 50))
  plan(multisession, workers = 10)
  future_map(name_batches, worms_for_names, .options = furrr_options(seed = NULL)) %>%
    bind_rows()
}

#' Run a BLAST search
#'
#' @param sequence A DNA sequence
#' @param max_num_seq The maximum number of sequences to return
#' @param verbose Whether to print status messages
#' @param add_worms Whether to add WoRMS taxonomy information
#' @return A data frame with the BLAST results
#' @export
blast <- function(sequence, max_num_seq = 100, verbose = FALSE, add_worms = FALSE) {
  rid <- submit_sequence(sequence, max_num_seq = max_num_seq)
  status_ready <- wait_until_ready(rid, verbose)
  if (status_ready) {
    results <- get_results(rid)
    results <- results %>%
      mutate(across(c(len, score, evalue, identity, query_from, query_to, hit_from, hit_to, align_len), as.numeric)) %>%
      mutate(
        perc_identity = identity / align_len,
        query_cover = (query_to - query_from) / nchar(sequence)
      )
    if (add_worms) {
      species_names <- unique(results$sciname)
      matches <- get_taxonomy(species_names)
      results <- results %>%
        left_join(matches %>% select(input, phylum, class, order, family, genus, scientificname, rank), by = c("sciname" = "input")) %>%
        mutate(species = if_else(rank == "Species", scientificname, NA))
    }
    return(results)
  } else {
    stop("BLAST failed")
  }
}
