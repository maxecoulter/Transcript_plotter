#library(shiny)
library(dplyr)
library(stringr)
library(Cairo)
options(shiny.usecairo=T)










source('genemodel_adjusted2.r')
#input <- 'genemodel_input.csv'
#transcripts <- read.csv(input)

#gene_info <- read.csv('Shiny_out.csv')
#rownames(gene_info) <- gene_info$gene

#domain_input <- 'genemodel_input_domain.csv'
#domains <- read.csv(domain_input)

#snpeff_input <- read.csv("snpeff_out.csv")

gene_choice <- "NA"

plot_snpeff <- function(snpeff=snpeff){
  if (nrow(snpeff) == 0){
    return()
  }
  for (i in 1:nrow(snpeff)){
    start <- as.numeric(snpeff[i,"position_start"])
    end <- as.numeric(snpeff[i,"position_end"])
    label <- snpeff[i,"label"]
    depth <- snpeff[i,"depth"]
    colour <- snpeff[i,"colour"]
    mutation.plot(start, end, text=label, col='black', drop=depth, haplotypes=c(colour))
    
  }
}

fig_res <- 200 #Figure resolution


ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    fileInput("input", "Choose CSV file for genemodel input (genemodel_input.csv)", accept = ".csv"),
    fileInput("shiny_input", "Choose CSV file for transcript information (shiny_out.csv)", accept = ".csv"),
    fileInput("domain_input", "Choose CSV file for interproscan input (optional, genemodel_input_domain.csv)", accept = ".csv"),
    fileInput("snpeff_input", "Choose CSV file for snpeff input (optional, snpeff_out.csv)", accept = ".csv"),
    sliderInput(inputId = "axis_text_size", "Axis text size:",
                min = 0, max = 2,
                value = 1, step = 0.1, width = "300px"),
    sliderInput(inputId = "legend_gap", "Legend gap:",
                min = 0, max = 10,
                value = 1, step = 0.5, width = "300px"),
    sliderInput(inputId = "legend_size", "Legend size:",
                min = 0.1, max = 4,
                value = 1.4, step = 0.1, width = "200px"),
    sliderInput(inputId = "width", "Plot width:",
                min = 50, max = 1000,
                value = 500, step = 10, width = "200px"),
    sliderInput(inputId = "height", "Plot height:",
                min = 50, max = 1000,
                value = 400, step = 10, width = "200px"),
   
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   br(),
   
   sliderInput(inputId = "t_width", "Transcript plot width:",
               min = 50, max = 1000,
               value = 500, step = 10, width = "200px"),
   sliderInput(inputId = "t_height", "Transcript plot height:",
               min = 50, max = 1000,
               value = 200, step = 10, width = "200px")
  ),
  mainPanel(
    numericInput("browser_start", "Start", 0, min = 0, max = NA, step = 10),
    numericInput("browser_end", "End", 0, min = 0, max = NA, step = 10),
    plotOutput("browser_plot", width=700, height=700),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    br(),
    selectInput(inputId = "gene", 
                  label = "Choose a gene", choices = gene_choice),
                
                selectInput(
                  inputId = "colours",
                  label = "Select a list of domain colours",
                  choices = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
                  multiple = TRUE,
                  selectize = TRUE,
                ),
                actionButton(inputId="colour_button", label = "Add colours manually?"),
                actionButton(inputId="button", label = "Add domains"),
                actionButton(inputId="off_button", label = "Remove domains"),
                
              
                plotOutput("gene_plot", width=500, height=400),
                downloadButton("download"),
                br(),
                br(),
                br(),
                br(),
                #Now view transcripts...
                selectInput(
                  inputId = "transcript",
                  label = "Select a transcript to view",
                  choices = c(),
                  multiple = FALSE,
                  selectize = TRUE
                ),
                actionButton(inputId="transcript_button", label = "View select transcript"),
                sliderInput(inputId = "t_axis_text_size", "Axis text size:",
                            min = 0, max = 2,
                            value = 1, step = 0.1),
                sliderInput(inputId = "t_legend_gap", "Legend gap:",
                            min = 0, max = 10,
                            value = 1, step = 0.5),
                selectInput(
                  inputId = "t_colours",
                  label = "Select a list of domain colours",
                  choices = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)],
                  multiple = TRUE,
                  selectize = TRUE,
                ),
                actionButton(inputId="t_colour_button", label = "Add colours manually?"),
                actionButton(inputId="t_button", label = "Add domains"),
                actionButton(inputId="t_off_button", label = "Remove domains"),
                actionButton(inputId="snpeff_button", label = "Add snpeff data"),
                actionButton(inputId="snpeff_offbutton", label = "Remove snpeff data"),
                plotOutput("transcript_plot"),
                downloadButton("transcript_download")
                
                )))

#g_domains <- domains %>% filter(gene == 'BaRT2v18chr3HG123260')             









server <- function(input, output){
 observeEvent(input$input,{
   observeEvent(input$shiny_input,{
   infile <- input$input
   transcripts <- read.csv(infile$datapath)
   gene_choice <- as.list(unique(transcripts$gene))
   updateSelectInput(inputId='gene', label = "Choose a gene", choices = gene_choice,
                     selected = NULL)
   
   data<- reactive({transcripts %>% filter(gene == input$gene)})
   
   shiny <- input$shiny_input
   gene_info <- read.csv(shiny$datapath)
   rownames(gene_info) <- gene_info$gene
   
   

  #data2 <- reactive({transcripts_test() %>% filter(gene == input$gene)})
  transcript_vector <- reactive({transcripts %>% filter(gene == input$gene) %>% pull(transcript) %>% unique() %>% as.vector()})
  
  #g_domains <- eventReactive(input$button, {domains %>% filter(gene == input$gene)})
  
  plotInput <- function() {
    gene.transcript.model.plot(model=data(), gene_start=gene_info[input$gene,"start"],
                               gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"],
                               xaxis=T, gap=0.2, legend_gap=input$legend_gap, axis_text_size=input$axis_text_size)
  }
  observe({
    first_start <- gene_info[input$gene,"start"]
    first_end <- gene_info[input$gene,"end"]
    updateNumericInput(session = getDefaultReactiveDomain(), inputId='browser_start', "Start", first_start, min = 0, max = NA, step = 10)
    updateNumericInput(session = getDefaultReactiveDomain(), inputId='browser_end', "End", first_end, min = 0, max = NA, step = 10)})
  output$browser_plot <- renderPlot({Transcript.Browser(gene_info=gene_info, transcripts=transcripts,
                                                        plot_start=input$browser_start,plot_end=input$browser_end)}, width= 700, height= 700, res = fig_res)
  
  
  output$gene_plot <- renderPlot({plotInput()}, width=function() input$width, height=function() input$height, res = fig_res)
  download_plot <- function(plotInput){downloadHandler(
    filename = function() {
      paste0(input$gene, ".eps")
    },
    content = function(file) {
      saveplot <- function(file){
        setEPS()
        postscript(file=file, pointsize=30, width=input$width * 0.013, height=input$height * 0.013)
        plotInput()
        dev.off()
      }
      saveplot(file)
    }
  )}
  
   output$download <- download_plot(plotInput)

        
      
      
  
  observeEvent(input$domain_input,{
    domain_input <- input$domain_input
    domains <- read.csv(domain_input$datapath)
  g_domains <- reactive({domains %>% filter(gene == input$gene)})
  #output$gene_plot <- renderPlot({gene.transcript.model.plot(model=data(), gene_start=gene_info[input$gene,"start"],gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"],xaxis=T, gap=0.2,
                                                      #legend_gap=1,axis_text_size=0.8)
                                 # })
 
  observeEvent(input$button, {
   
    g_domains <- domains %>% filter(gene == input$gene)
    plotInput_domain <- function() {gene.transcript.model.plot(model=data(), gene_start=gene_info[input$gene,"start"],gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"],xaxis=T, gap=0.2, legend_gap=input$legend_gap,axis_text_size=input$axis_text_size)
      gene.transcript.model.plot_domain(model=g_domains, 
                                        transcript_vector=transcript_vector(), gene_start=gene_info[input$gene,"start"], gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"], gap=0.2, legend_gap=input$legend_gap, legend_size = input$legend_size)
      }
    output$gene_plot <- renderPlot({plotInput_domain()
    },width=function() input$width, height=function() input$height, res = fig_res)
    
    output$download <- download_plot(plotInput_domain)
    #  output$download <- downloadHandler(
    #   filename = function() {
    #     paste0(input$gene, ".eps")
    #   },
    #   content = function(file) {
    #     saveplot <- function(file){
    #       setEPS()
    #       postscript(file=file, pointsize=30)
    #       plotInput_domain()
    #       dev.off()
    #     }
    #     saveplot(file)
    #   }
    # )
    }
  )
  
  observeEvent(input$off_button, {
    g_domains <- domains %>% filter(gene == input$gene)
    
    #output$gene_plot <- renderPlot({gene.transcript.model.plot(model=data(), gene_start=gene_info[input$gene,"start"],gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"],xaxis=T, gap=0.2, legend_gap=input$legend_gap,axis_text_size=input$axis_text_size)
    #},width=function() input$width, height=function() input$height, res = fig_res)})
    output$gene_plot <- renderPlot({plotInput()
      },width=function() input$width, height=function() input$height, res = fig_res)
    output$download <- download_plot(plotInput)
    })
    observeEvent(input$colour_button, {
      plotInput_domain_colours <- function() {gene.transcript.model.plot(model=data(), gene_start=gene_info[input$gene,"start"],gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"],xaxis=T, gap=0.2, legend_gap=input$legend_gap,axis_text_size=input$axis_text_size)
        gene.transcript.model.plot_domain(model=g_domains, 
                                          transcript_vector=transcript_vector(), gene_start=gene_info[input$gene,"start"], gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"], gap=0.2, legend_gap=input$legend_gap, legend_size = input$legend_size, manual_colours=input$colours)
      }
    g_domains <- domains %>% filter(gene == input$gene)
    output$gene_plot <- renderPlot({plotInput_domain_colours()
    },width=function() input$width, height=function() input$height, res = fig_res)
   output$download <- download_plot(plotInput_domain_colours)
   # output$gene_plot <- renderPlot({gene.transcript.model.plot(model=data(), gene_start=gene_info[input$gene,"start"],gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"],xaxis=T, gap=0.2, legend_gap=input$legend_gap,axis_text_size=input$axis_text_size)
     # gene.transcript.model.plot_domain(model=g_domains, 
                                        #transcript_vector=transcript_vector(), gene_start=gene_info[input$gene,"start"], gene_bpstop=gene_info[input$gene,"end"], orientation=gene_info[input$gene,"orientation"], gap=0.2, legend_gap=input$legend_gap,manual_colours=input$colours, legend_size = input$legend_size)
   # },width=function() input$width, height=function() input$height, res = fig_res)
    
    })
  })
  observe({
  updateSelectInput(inputId='transcript', label = "Select a transcript to view", choices = transcript_vector(),
                    selected = NULL)})
  observeEvent(input$transcript_button, {
    input_transcript <- input$transcript
    download_t_plot <- function(plotInput){downloadHandler(
      filename = function() {
        paste0(input_transcript, ".eps")
      },
      content = function(file) {
        saveplot <- function(file){
          setEPS()
          postscript(file=file, pointsize=30, width=(input$t_width * 0.013), height=(input$t_height * 0.013))
          plotInput()
          dev.off()
        }
        saveplot(file)
      }
    )}
    transcript <- transcripts %>% filter(transcript == input$transcript)
    t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
    t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])
    t_plotInput <- function(){genemodel.plot(model=transcript, start=t_start, 
                                             bpstop=t_end, orientation=gene_info[input$gene,"orientation"], xaxis=T, axis_text_size=input$t_axis_text_size, legend_gap=input$t_legend_gap)}
    output$transcript_plot <- renderPlot({t_plotInput()},
      width=function() input$t_width, height=function() input$t_height, res = fig_res)
    output$transcript_download <- download_t_plot(t_plotInput)
  observeEvent(input$snpeff_input,{
    s_input <- input$snpeff_input
    snpeff_input <- read.csv(s_input$datapath)
    observeEvent(input$snpeff_button, {
      transcript <- transcripts %>% filter(transcript == input_transcript)
      t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
      t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])  
      #
      snpeff <- snpeff_input %>% filter(Transcript == input_transcript)
      t_plotInput_snpeff <- function(){genemodel.plot(model=transcript, start=t_start, 
                                                      bpstop=t_end, orientation=gene_info[input$gene,"orientation"], xaxis=T, axis_text_size=input$t_axis_text_size,legend_gap=input$t_legend_gap)
        plot_snpeff(snpeff=snpeff)}
      output$transcript_plot <- renderPlot({t_plotInput_snpeff()}, width=function() input$t_width, height=function() input$t_height, res = fig_res)
      
      output$transcript_download <- download_t_plot(t_plotInput_snpeff)  
    })
    observeEvent(input$snpeff_offbutton, {
      transcript <- transcripts %>% filter(transcript == input_transcript)
      t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
      t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])  
      #
      
      output$transcript_plot <- renderPlot({t_plotInput()},width=function() input$t_width, height=function() input$t_height, res = fig_res)
    output$transcript_download <- download_t_plot(t_plotInput)
  })})
    
  observeEvent(input$domain_input,{
    domain_input <- input$domain_input
    domains <- read.csv(domain_input$datapath)
    t_plotInput_domain <- function(){genemodel.plot(model=transcript, start=t_start, 
                                                    bpstop=t_end, orientation=gene_info[input$gene,"orientation"], xaxis=T, axis_text_size=input$t_axis_text_size,legend_gap=input$t_legend_gap)
      t_domains <- domains %>% filter(transcript == input_transcript)
      genemodel.plot_domain(model=t_domains, start=t_start, bpstop=t_end, orientation=gene_info[input$gene,"orientation"])}
  observeEvent(input$t_button, {
    transcript <- transcripts %>% filter(transcript == input_transcript)
    t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
    t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])
    
    output$transcript_plot <- renderPlot({t_plotInput_domain()},
                                         width=function() input$t_width, height=function() input$t_height, res = fig_res)
    output$transcript_download <- download_t_plot(t_plotInput_domain)
    
    observeEvent(input$t_colour_button, {
      transcript <- transcripts %>% filter(transcript == input_transcript)
      t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
      t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])
      t_plotInput_domain_colours <- function(){genemodel.plot(model=transcript, start=t_start, 
                                                      bpstop=t_end, orientation=gene_info[input$gene,"orientation"], xaxis=T, axis_text_size=input$t_axis_text_size,legend_gap=input$t_legend_gap)
        t_domains <- domains %>% filter(transcript == input_transcript)
        genemodel.plot_domain(model=t_domains, start=t_start, bpstop=t_end, manual_colours=input$t_colours, orientation=gene_info[input$gene,"orientation"])}
      output$transcript_plot <- renderPlot({t_plotInput_domain_colours()},
                                           width=function() input$t_width, height=function() input$t_height, res = fig_res)
      output$transcript_download <- download_t_plot(t_plotInput_domain_colours)
    })})
  observeEvent(input$t_off_button, {
    transcript <- transcripts %>% filter(transcript == input_transcript)
    t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
    t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])                   
    output$transcript_plot <- renderPlot({t_plotInput()}, width=function() input$t_width, height=function() input$t_height, res = fig_res)
      
    output$transcript_download <- download_t_plot(t_plotInput)
    })
  
  observeEvent(input$snpeff_input,{
    s_input <- input$snpeff_input
    snpeff_input <- read.csv(s_input$datapath)
  observeEvent(input$snpeff_button, {
    transcript <- transcripts %>% filter(transcript == input_transcript)
    t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
    t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])  
    #
    snpeff <- snpeff_input %>% filter(Transcript == input_transcript)
    
    t_plotInput_domain_snpeff <- function(){genemodel.plot(model=transcript, start=t_start, 
                                                           bpstop=t_end, orientation=gene_info[input$gene,"orientation"], xaxis=T, axis_text_size=input$t_axis_text_size,legend_gap=input$t_legend_gap)
      t_domains <- domains %>% filter(transcript == input_transcript)
      genemodel.plot_domain(model=t_domains, start=t_start, bpstop=t_end, manual_colours=input$t_colours, orientation=gene_info[input$gene,"orientation"])
      plot_snpeff(snpeff=snpeff)}
    output$transcript_plot <- renderPlot({t_plotInput_domain_snpeff()}, width=function() input$t_width, height=function() input$t_height, res = fig_res)
    output$transcript_download <- download_t_plot(t_plotInput_domain_snpeff)
  })
  
  observeEvent(input$snpeff_offbutton, {
    transcript <- transcripts %>% filter(transcript == input_transcript)
    t_start <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[1])
    t_end <- as.numeric(stringr::str_split_fixed(transcript[1,"transcript_coordinates"], "-", 2)[2])  
    #
   
    #output$transcript_plot <- renderPlot({genemodel.plot(model=transcript, start=t_start, 
                                                         #bpstop=t_end, orientation=gene_info[input$gene,"orientation"], xaxis=T, axis_text_size=input$t_axis_text_size,legend_gap=input$t_legend_gap)
      #t_domains <- domains %>% filter(transcript == input_transcript)
      #genemodel.plot_domain(model=t_domains, start=t_start, bpstop=t_end, manual_colours=input$t_colours, orientation=gene_info[input$gene,"orientation"])
    output$transcript_plot <- renderPlot({t_plotInput_domain()}, width=function() input$t_width, height=function() input$t_height, res = fig_res)
    output$transcript_download <- download_t_plot(t_plotInput_domain)
  })
  
 })})})
 })
})
}






shinyApp(ui = ui, server = server)