# Shared Word-table styling for the supplementary tables
# (04_build_summary_tables.R: Tables S-A/S-B; 08_watermass_regression.R:
# Table S-C), so every table has the same look: landscape page, Times New
# Roman, booktabs-style rules (heavy top/bottom, light under the header and
# between response blocks), a shaded bold-italic section row, and a PDF copy
# rendered by Word.

library(officer)
library(flextable)

TABLE_FONT <- "Times New Roman"
set_flextable_defaults(font.family = TABLE_FONT, font.size = 8, padding.top = 1,
                       padding.bottom = 1, padding.left = 3, padding.right = 3)
rule_heavy <- fp_border(color = "black", width = 1.5)
rule_light <- fp_border(color = "black", width = 0.5)

#' one table block. df is pre-formatted text; `block_col` is blank on the
#' continuation rows of a block, and a light rule separates the blocks.
#' header_labels and col_width are named by column key, in display order.
booktabs_ft <- function(df, header_labels, col_width, left_cols, section_label,
                        block_col = "group") {
  keys <- names(header_labels)
  ft <- flextable(df, col_keys = keys)
  ft <- set_header_labels(ft, values = as.list(header_labels))
  ft <- width(ft, j = keys, width = col_width[keys])
  ft <- align(ft, j = setdiff(keys, left_cols), align = "right", part = "all")
  ft <- align(ft, j = intersect(left_cols, keys), align = "left", part = "all")
  ft <- bold(ft, part = "header")
  ft <- valign(ft, valign = "bottom", part = "header")
  ft <- valign(ft, valign = "top", part = "body")
  ft <- border_remove(ft)

  # section label as a shaded, full-width first body row
  ft <- add_body_row(ft, values = section_label, colwidths = length(keys), top = TRUE)

  # light rule above each new block; +1 for the section row added above
  starts <- which(df[[block_col]] != "") + 1
  if (length(starts) > 1) ft <- hline(ft, i = starts[-1] - 1, border = rule_light, part = "body")
  ft <- bg(ft, i = 1, bg = "#EDEDED", part = "body")
  ft <- bold(ft, i = 1, part = "body")
  ft <- italic(ft, i = 1, part = "body")
  ft <- align(ft, i = 1, align = "left", part = "body")

  ft <- hline_top(ft, border = rule_heavy, part = "header")
  ft <- hline_bottom(ft, border = rule_light, part = "header")
  ft <- hline_bottom(ft, border = rule_heavy, part = "body")
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- font(ft, fontname = TABLE_FONT, part = "all")
  set_table_properties(ft, layout = "fixed")
}

table_txt <- function(size, bold = FALSE) fp_text(font.family = TABLE_FONT, font.size = size, bold = bold)
add_caption <- function(doc, label, text) {
  body_add_fpar(doc, fpar(ftext(paste0(label, " "), table_txt(10, bold = TRUE)),
                          ftext(text, table_txt(10)), fp_p = fp_par(padding.bottom = 4)))
}
add_note <- function(doc, text) {
  body_add_fpar(doc, fpar(ftext(text, table_txt(8)), fp_p = fp_par(padding.top = 3)))
}
add_gap <- function(doc) body_add_fpar(doc, fpar(ftext("", table_txt(3))))

landscape_doc <- function() {
  body_set_default_section(read_docx(), prop_section(
    page_size = page_size(orient = "landscape", width = 8.5, height = 11),
    page_margins = page_mar(top = 0.4, bottom = 0.4, left = 0.5, right = 0.5,
                            header = 0.3, footer = 0.3, gutter = 0)))
}

#' write the docx, plus a PDF copy rendered by Microsoft Word through the
#' doconv package (Windows or macOS) so it always matches the docx; the PDF
#' is skipped with a message if doconv or Word is not available
save_docx_and_pdf <- function(doc, docx_path) {
  print(doc, target = docx_path)
  cat("wrote", docx_path, "\n")
  pdf_path <- sub("\\.docx$", ".pdf", docx_path)
  if (requireNamespace("doconv", quietly = TRUE) && doconv::msoffice_available()) {
    doconv::docx2pdf(input = docx_path, output = pdf_path)
    cat("wrote", pdf_path, "\n")
  } else {
    message("PDF not written: needs the doconv package and Microsoft Word ",
            "(install.packages(\"doconv\")); ", basename(pdf_path), " may be out of date.")
  }
}
