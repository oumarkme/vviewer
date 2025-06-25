### Load required packages
packages <- c('data.table', 'tidyverse', 'Rsamtools', 'bslib', 'shiny', 'GenomicRanges', 'IRanges')
for(pkg in packages) {
    if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos='https://cloud.r-project.org/')
        library(pkg, character.only = TRUE)
    }
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
                card_header("Physical position"),
                textInput( 
                    "coordinate", 
                    "Coordinate",
                    placeholder = "(ex: chr1:1000-2000 or chr1:1000)",
                    value = NULL
                )
            ),
            card(
                card_header("Sample"),
                selectizeInput(
                    'full_sample_list',
                    'Select samples',
                    NULL,
                    multiple = TRUE
                )
                
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


    # Get coordinate input and remove all commas
    coordinate <- reactive({
        coord = input$coordinate
        req(coord)
        coord = gsub(",", "", coord)
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
