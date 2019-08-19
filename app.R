library(shiny)
library(shinydashboard)
library(datasets)
library(ggplot2)
library(shinythemes)
library(tidyr)
library(reshape2)
library(data.table)
library(dplyr)

commonTheme = list(labs(color="Density",fill="Density"),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

#tcga_rna <- read.delim("~/Google Drive/epigenomics_drive/tcga/laml_tcga_pub/laml_tcga_pub/data_RNA_Seq_v2_expression_median.txt")
rnamelt <- data.frame(fread("./rnamelt.txt"))
tcga_mut <- data.frame(fread("./tcga_mut.txt"))

# Define UI for application that draws a histogram
ui <- navbarPage(title="TCGA RNA-seq Explorer", theme = shinythemes::shinytheme("sandstone"),
                 tabPanel(title = "TCGA Genes of Interest vs. Mutational Status",
                          fluidPage(
                            # Sidebar with controls to select the variable to plot against mpg
                            # and to specify whether outliers should be included
                            sidebarPanel(
                              selectInput("muts", "Separate by Mutational Status", choices = as.character(unique(tcga_mut$mutation)), multiple = T, selected = as.character(unique(tcga_mut$mutation))[c(4,50,258)]),
                              selectInput("geneselect", "Select Genes to Plot", choices = as.character(unique(rnamelt$gene)), multiple = T, selected = c("HOXA9", "CSF2RA", "KIT")),
                              tableOutput("presult")
                              #tableOutput("pvals")
                              #selectInput("fillm", "Colour genes in vitaminC region:", names(data), multiple = F, selected = "dxmr_2kbTSS"),
                            ),
                            # Show the caption and plot of the requested variable against mpg
                            mainPanel(
                              plotOutput("mutplot", 
                                         width = "100%",  height = "100%",
                                         hover = "plot_hover",
                                         click = "plot_click")
                            )
                          ))
)


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
  mutreact <- reactive({
    nmut <- tbl_df(tcga_mut) %>% select(-aa_change) %>% unique() %>% filter(mutation %in% input$muts) %>% group_by(id) %>% tally()
    tbl_df(tcga_mut) %>%
      filter(mutation %in% input$muts) %>% 
      inner_join(nmut) %>% 
      mutate(status=ifelse(n < 2, as.character(mutation), "multiple")) %>% 
      select(id, status) %>% 
      unique() %>%
      inner_join(rnamelt) %>% filter(gene %in% input$geneselect)
  })
  
#  muts <- c("IDH1", "IDH2")
#  geneselect <- c("HOXA9", "CSF2RA", "PDGFRB")
  mutreactb <- reactive({
     # nmut <-  tbl_df(tcga_mut) %>% select(-aa_change) %>% unique() %>% filter(mutation %in% muts) %>% group_by(id) %>% tally()
    #nmut <- tbl_df(tcga_mut) %>% select(-aa_change) %>% unique() %>% mutate(mutant=ifelse(mutation %in% input$muts, yes = mutation, no = "other")) %>% select(mutant, id) %>% unique()
    nmut <- tbl_df(tcga_mut) %>% select(-aa_change) %>% unique() %>% filter(mutation %in% input$muts) %>% group_by(id) %>% tally()
    mutid <- tbl_df(tcga_mut) %>%
      filter(mutation %in% input$muts) %>% 
      inner_join(nmut) %>% 
      mutate(status=ifelse(n < 2, as.character(mutation), "multiple")) %>% 
      select(id, status) %>% 
      unique() %>% data.frame()

    tbl_df(tcga_mut) %>% select(-aa_change, -mutation) %>% unique() %>% full_join(mutid, by="id") %>% 
      mutate(mutant=ifelse(is.na(status), yes = "other", no = status)) %>% select(-status) %>%
      select(mutant, id) %>%
      inner_join(rnamelt) %>% filter(gene %in% input$geneselect) %>% mutate(logtag=log(1+tagcount))
  })
  
  
  output$infod <- renderDataTable({
      scatreact() %>% tbl_df() %>% nearPoints(input$plot_click, maxpoints = 25, threshold = 10)
    })
    
  output$infor <- renderDataTable({
    mutreactb() %>% tbl_df() %>% nearPoints(input$plot_click, maxpoints = 25, threshold = 10)
  })
  
  
  
  
  output$mutplot <- renderPlot({
    mutreactb() %>%
        ggplot(aes(x=mutant, y=logtag))+
        geom_boxplot(position="dodge", width=0.7, fill="honeydew3", alpha=0.3, outlier.colour = "white")+
        geom_point(data=mutreactb(), aes(x=mutant, y=logtag, colour=mutant), position="jitter", size=0.75)+
        ylab("Expression log(1+tagcount)")+
        xlab("")+
      theme_bw()+
        facet_wrap(~gene, scales="free")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      guides(color=F)+
      scale_colour_brewer(palette = "Set2", name="Mutant Status")
    }, width = 900, height = 450)
  
  output$mutable <- renderDataTable({
    # renderPrint({
    # if(is.null(input$plot_hover)) return("hover over a point")
    # With base graphics, need to tell it what the x and y variables are.
    mutreactb() %>% group_by(gene, mutant) %>% summarize(meanlog1plus=mean(log(1+tagcount)), sdlog1plus=sd(log(1+tagcount)))
    # nearPoints() also works with hover and dblclick events
    #})
  })
  
  output$caption <- renderText({ "# of patients with mutations in:" })
  
  output$presult <- renderTable({
    print <- mutreactb() %>% select(id,mutant) %>% unique() %>% group_by(mutant) %>% tally() %>% data.frame()
    colnames(print) <- c("Mutated Gene", "# of patients")
    print
  })
  
  
  output$stats <- renderTable({
  mutreact() %>%
    select(-id) %>%
    group_by(status, gene) %>% 
    summarise(value = list(tagcount)) %>% 
    spread(gene, value) %>% 
    group_by(status) %>% 
    mutate(p_value = t.test(unlist(input$ttesta), unlist(input$ttestb))$p.value,
           t_value = t.test(unlist(input$ttesta), unlist(input$ttestb))$statistic)
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)



