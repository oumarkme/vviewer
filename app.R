### Load required packages
packages <- c('data.table', 'tidyverse', 'Rsamtools', 'bslib', 'shiny', 'GenomicRanges', 'IRanges', 'httr', 'jsonlite', 'xml2')
for(pkg in packages) {
    library(pkg, character.only = TRUE)
}

### vcfFile
vcf_paths <- list.files('./data', pattern = '\\.vcf(\\.gz)?$', full.names = TRUE)
vcf_files <- basename(vcf_paths)
vcf_file_list <- setNames(as.list(vcf_paths), vcf_files)

### User Interface
ui = page_fillable(
    titlePanel("Variants viewer"),

    ## Option selection UI
    layout_columns(
            card(
                card_header("Choose file"),
                radioButtons(
                    'vcfFile',
                    'Select VCF file from the list',
                    vcf_file_list
                )
            ),
            card(
                card_header("Variant coordinate input"),
                navset_card_tab(
                    nav_panel("Physical position", 
                        textInput( 
                        "coordinate", 
                        "Coordinate (hg38)",
                        placeholder = "(ex: chr22:42130692 or chr22:42130692-42130692)",
                        value = NULL
                        )
                    ),
                    nav_panel("HGVS notation", 
                        textInput(
                            "refseq_coordinate",
                            "Identifier based on a HGVS notation",
                            placeholder = "(ex: NM_000106.6:c.100C>T)",
                            value = NULL
                        )
                    )
                ),
                textOutput("coord")
            ),

            card(
                card_header("Sample"),
                selectizeInput(
                    'full_sample_list',
                    'Select samples',
                    NULL,
                    multiple = TRUE
                ),
                p('Remain empty to select all samples.')
            ),
    fill=F, min_height='200px'),

    ## Action button
    actionButton("applyChanges", "Apply Changes", icon("refresh")),

    ## Variant Table UI
    layout_columns(
        card(
            card_header("Variants Table"),
            tableOutput("variant_table") 
        )
    )
)

### Server 
server = function(input, output, session){
    
    # Get vcf file path based on user selection
    vcf_file_path <- reactive({
        vcfFile = input$vcfFile
        req(vcfFile)
        return(vcfFile)
    })
    output$vcfFilePath = renderText({ vcf_file_path() })


    # Dynamically update the sample list when a VCF file is selected
    observeEvent(vcf_file_path(), {
        current_vcf_path <- vcf_file_path()
        req(current_vcf_path)

        # Read VCF header to get sample names
        samples <- as.character(Rsamtools::scanBcfHeader(current_vcf_path)[[current_vcf_path]]$Sample)
        
        # Update the selectInput choices
        updateSelectInput(session, "full_sample_list", choices = samples, selected = NULL)
    })


    # Get coordinate input (from either coordinate or refseq_coordinate)
    coordinate <- reactive({
        coord = input$coordinate
        refseq = input$refseq_coordinate

        # Prefer coordinate if provided, otherwise use refseq_coordinate
        if (!is.null(coord) && nzchar(coord)) {
            coord = gsub(",", "", coord)

        } else if (!is.null(refseq) && nzchar(refseq)) {
            ensembl_get = GET(paste0('https://rest.ensembl.org/vep/human/hgvs/', refseq), 
                content_type("application/json"))
            ensembl_result = fromJSON(toJSON(content(ensembl_get)))

            if (!is.null(ensembl_result)){
                coord = paste0('chr', ensembl_result$seq_region_name[[1]][1], ':', 
                        ensembl_result$start[[1]][1], '-',
                        ensembl_result$end[[1]][1])
            }
        }

        return(coord)
    })
    output$coord = renderText({ coordinate() })


    # Get the VCF file and read it
    variants_df <- eventReactive(input$`applyChanges`, {
        vcf_path <- vcf_file_path()
        range = coordinate()
        selected_samples = input$full_sample_list

        # Prepare range information
        splitRange = stringr::str_split(range, pattern = ":|-", simplify = T, n=3)[1, ]
        chr = splitRange[1]
        start = as.numeric(splitRange[2])
        end = ifelse(splitRange[3]=='', start, as.numeric(splitRange[3]))
        width = end - start + 1
        gRange = GenomicRanges::GRanges(seqnames=chr, IRanges::IRanges(start, end, width))

        # Load variants from the vcf file
        varTmp = unlist(Rsamtools::scanTabix(vcf_path, param=gRange))

        # Make the variant result a data frame
        if (length(varTmp) == 1) {
            varTmp = as.data.table(t(unlist(str_split(varTmp, '\t'))))
        } else if (length(varTmp) > 1){
            as.data.frame(varTmp)
            varTmp = fread(text=varTmp)
        } else {
            varTmp = NULL
        }

        # Variant info
        INFOS = varTmp %>% select(V1, V2, V4, V5)
        colnames(INFOS) = c("CHROM", "POS", "REF", "ALT")

        # Genotype
        GENOS = varTmp %>% 
            select(-V1, -V2, -V3, -V4, -V5, -V6, -V7, -V8, -V9) %>%
            mutate(across(everything(), ~str_extract(., "^[01\\./\\|]+")))
        colnames(GENOS) = as.character(Rsamtools::scanBcfHeader(vcf_path)[[vcf_path]]$Sample)

        # Select samples if provided
        if (!is.null(selected_samples) && length(selected_samples) > 0) {
            GENOS = GENOS %>% select(all_of(selected_samples))
        }

        # Return data frame
        varTmp = bind_cols(INFOS, GENOS)
        return(varTmp)
    })
    output$variant_table = renderTable({variants_df()}, striped = TRUE)

}

### Run the application
shinyApp(ui = ui, server = server)
