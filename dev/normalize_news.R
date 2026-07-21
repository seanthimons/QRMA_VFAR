normalize_news_text <- function(x) {
  x <- gsub("\u2018", "'", x, fixed = TRUE)
  x <- gsub("\u2019", "'", x, fixed = TRUE)
  replacements <- c(
    "depreciated" = "deprecated",
    "docuementation" = "documentation",
    "developement" = "development",
    "calculatiosn" = "calculations",
    "debuggin" = "debugging",
    "Beta Poisson" = "Beta-Poisson",
    "beta Poisson" = "beta-Poisson"
  )
  for (from in names(replacements)) {
    x <- gsub(from, replacements[[from]], x, fixed = TRUE)
  }
  x
}

normalize_autonewsmd_repo_list <- function(an) {
  text_columns <- c("summary", "message", "clean_summary")

  for (tag in names(an$repo_list)) {
    commits <- an$repo_list[[tag]]$commits
    for (column in intersect(text_columns, names(commits))) {
      commits[[column]] <- normalize_news_text(commits[[column]])
    }
    an$repo_list[[tag]]$commits <- commits
  }

  invisible(an)
}

normalize_news_file <- function(path = "NEWS.md") {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  writeLines(normalize_news_text(lines), path, useBytes = TRUE)
  invisible(path)
}

build_normalized_news <- function(repo_name, repo_path = getwd(), path = "NEWS.md") {
  if (!requireNamespace("autonewsmd", quietly = TRUE)) {
    stop("Package 'autonewsmd' is required to build NEWS.md.", call. = FALSE)
  }

  an <- autonewsmd::autonewsmd$new(
    repo_name = repo_name,
    repo_path = repo_path
  )
  an$generate()
  normalize_autonewsmd_repo_list(an)
  an$write(force = TRUE)
  news_path <- if (grepl("^([A-Za-z]:)?[/\\\\]", path)) path else file.path(repo_path, path)
  normalize_news_file(news_path)
}
