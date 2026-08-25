# Parse-annotations
Kristina Gagalova

``` r
library(tidyr)
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(stringr)
```

``` r
corr <- read.csv(gzfile("iwgsc_refseq_all_correspondances.csv.gz"), 
                 sep = " ")
anno_v2 <- read.csv(gzfile("iwgsc_refseqv2.1_functional_annotation.csv.gz"))
anno_v1 <- read.csv(gzfile("Triticum_aestivum.IWGSC.pep.all.anno2.tsv.gz"),
                    sep = "\t",
                    header = F)

anno_master <- anno_v1 %>%
  mutate(
    # TAIR projection
    tair_accession = if_else(
    str_detect(V3, "Projected from Arabidopsis thaliana"),
    str_extract(V3, "(?<=\\().*?(?=\\))"),
    NA_character_
    ),
    oryza_accession = if_else(
      str_detect(V3, "Projected from Oryza sativa"),
      str_extract(V3, "(?<=\\().*?(?=\\))"),
      NA_character_
    ),
    # Extract UniProt part only if UniProt exists
    uniprot = if_else(
      str_detect(V3, "UniProtKB"),
      str_trim(str_extract(V3, "UniProtKB.*$")),
      NA_character_
    ),

    # ✅ UniProt parsing (STRICT)
    db = if_else(
      str_detect(uniprot, "^UniProtKB"),
      str_extract(uniprot, "^[^;]+"),
      NA_character_
    ),
    uniprot_accession = if_else(
      str_detect(uniprot, "^UniProtKB") & str_detect(uniprot, ";Acc:"),
      str_extract(uniprot, "(?<=;Acc:).*"),
      NA_character_
    ),

    # ✅ NCBI parsing (separate)
    ncbi_accession = if_else(
      str_detect(V3, "^NCBI gene") & str_detect(V3, ";Acc:"),
      str_extract(V3, "(?<=;Acc:).*"),
      NA_character_
    )
  ) %>%
  select(-db, -uniprot) %>%
  rename_with(~ c("v1.1_trans", "function", "source"), .cols = 1:3) %>%
  mutate(
    `v1.1` = str_remove(`v1.1_trans`, "\\..*")
  )

# clean and join
corr_clean <- corr %>%
  distinct(`v1.1`, .keep_all = TRUE)
anno_master <- anno_master %>%
  left_join(
    corr_clean %>% select("v1.1", "v2.1"),
    by = "v1.1" #%>%
    #select(-`v1.1`)
  )


anno_master <- anno_master %>%
  left_join(
    anno_v2,
    by = c("v2.1" = "g2.identifier")
  ) %>%
  select(-`v1.1`)
```

    Warning in left_join(., anno_v2, by = c(v2.1 = "g2.identifier")): Detected an unexpected many-to-many relationship between `x` and `y`.
    ℹ Row 3 of `x` matches multiple rows in `y`.
    ℹ Row 90567 of `y` matches multiple rows in `x`.
    ℹ If a many-to-many relationship is expected, set `relationship =
      "many-to-many"` to silence this warning.

``` r
# ---------------------------------- merge with Norin Cadenza mappings
norin_rbh <- read.csv(gzfile("Norin-CSv2.1-rbh.tsv.gz"),
                    sep = "\t",
                    header = F)
colnames(norin_rbh)[1:2] <- c("norin", "v2.1")

cadenza_rbh <- read.csv(gzfile("Cadenza-CSv2.1-rbh.tsv.gz"),
                    sep = "\t",
                    header = F)
colnames(cadenza_rbh)[1:2] <- c("cadenza", "v2.1")

norin_rbh2 <- norin_rbh %>%
  select("norin", "v2.1")  %>%
  rename_with(~ c("norin", "v2.1"), .cols = 1:2) %>%
  mutate(
    `v2.1` = str_remove(`v2.1`, "\\..*")
  )  %>% 
  left_join(
    anno_master,
    by = c("v2.1")
  )
```

    Warning in left_join(., anno_master, by = c("v2.1")): Detected an unexpected many-to-many relationship between `x` and `y`.
    ℹ Row 1 of `x` matches multiple rows in `y`.
    ℹ Row 608346 of `y` matches multiple rows in `x`.
    ℹ If a many-to-many relationship is expected, set `relationship =
      "many-to-many"` to silence this warning.

``` r
cadenza_rbh2 <- cadenza_rbh %>%
  select("cadenza", "v2.1")  %>%
  mutate(
    `v2.1` = str_remove(`v2.1`, "\\..*")
  )  %>% 
  left_join(
    anno_master,
    by = c("v2.1")
  )
```

    Warning in left_join(., anno_master, by = c("v2.1")): Detected an unexpected many-to-many relationship between `x` and `y`.
    ℹ Row 1 of `x` matches multiple rows in `y`.
    ℹ Row 273333 of `y` matches multiple rows in `x`.
    ℹ If a many-to-many relationship is expected, set `relationship =
      "many-to-many"` to silence this warning.

``` r
# summary of annotations - CS v2.1
table(norin_rbh2$f.type)
```


           Blast-Hit-Accession              Gene Ontology 
                        132466                     344405 
    Human readable description                   InterPro 
                        123067                     387525 
                          Pfam 
                        170747 

``` r
norin_rbh2 %>%
  group_by(f.type) %>%
  summarise(n_unique = n_distinct(norin))
```

    # A tibble: 6 × 2
      f.type                     n_unique
      <chr>                         <int>
    1 Blast-Hit-Accession           63396
    2 Gene Ontology                 63396
    3 Human readable description    63396
    4 InterPro                      63396
    5 Pfam                          63396
    6 <NA>                          38659

``` r
table(cadenza_rbh2$f.type)
```


           Blast-Hit-Accession              Gene Ontology 
                         98899                     267866 
    Human readable description                   InterPro 
                         93854                     289434 
                          Pfam 
                        126192 

``` r
cadenza_rbh2 %>%
  group_by(f.type) %>%
  summarise(n_unique = n_distinct(cadenza))
```

    # A tibble: 6 × 2
      f.type                     n_unique
      <chr>                         <int>
    1 Blast-Hit-Accession           52840
    2 Gene Ontology                 52840
    3 Human readable description    52840
    4 InterPro                      52840
    5 Pfam                          52840
    6 <NA>                          36738

``` r
write.table(norin_rbh2,
            gzfile("Norin-CS-rbh-final.tsv.gz"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

write.table(cadenza_rbh2,
            gzfile("Cadenza-CS-rbh-final.tsv.gz"),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
```
