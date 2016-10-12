seqmonk_names := $(notdir $(wildcard ../reports/*))
seqmonk_names := $(filter-out README.txt,$(seqmonk_names))
seqmonk_names := $(filter-out w%,$(seqmonk_names))


%.html: %.Rmd
	Rscript -e "library(rmarkdown); render('$<', output_format='html_document')"

dirs:
	for n in $(seqmonk_names); do \
		if [ ! -e $$n ]; then \
			mkdir $$n && cp src/index.Rmd $$n; \
		fi \
	done

clean:
	rm -rf $(seqmonk_names)

rmd_files = $(filter-out src/%,$(wildcard */index.Rmd))
html_files = $(patsubst %.Rmd,%.html, $(rmd_files))

run: $(html_files)
