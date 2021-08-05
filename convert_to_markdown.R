# This script converts a word document to Rmarkdown
require(rmarkdown)

# Enter file name and location 
file=here::here("documents","2021_07_31 Statistics Tech Memo.docx")

#enter markdown output name and location 
md_file = here::here("documents","2021_07_31_Statistics_Tech_Memo.md")

# render as github flavored markdown (gfm)
pandoc_convert(file,to="gfm",output = md_file, citeproc = TRUE, options=c("--extract-media=."))
