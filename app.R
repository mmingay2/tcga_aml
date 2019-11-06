library(shiny)
library(shinydashboard)
library(datasets)
library(ggplot2)
library(shinythemes)
library(tidyr)
library(reshape2)
library(data.table)
library(dplyr)

#############################################################################################################################################################
#load datasets
# contains gene expression data
rnamelt <- data.frame(fread("./rnamelt.txt"))
# contains patient id <-> mutation status
tcga_mut <- data.frame(fread("./tcga_mut.txt"))

#############################################################################################################################################################

# Define UI for application that will displays options, plots and interact with our users
ui <- 
  ##Use fluid page to create a page and hold our sidebarPanel and mainPanel
  fluidPage(#give the page a title
                title = "TCGA Genes of Interest vs. Mutational Status",
                #lets apply a theme using the 'shinythemes' package we loaded'
                theme = shinythemes::shinytheme("sandstone"),
                #hold our input selection boxes within a "sidebarPanel" within our FluidPage
                        sidebarPanel(
                              #Define a select box for mutation status
                                selectInput(
                                      #our inputId will be how we access the value selected.
                                      #it is accessed using `input$muts` or more generally `input$inputId``
                                      inputId = "muts", label = "Separate by Mutational Status", #text to display above the box
                                      #this parameter is where we define the choices available          
                                      choices = as.character(unique(tcga_mut$mutation)), 
                                      #can the user select multiple options?
                                      multiple = TRUE, 
                                      #what will be selected by default?
                                      selected = as.character(unique(tcga_mut$mutation))[c(4,50,258)]),
                                #Same as above but assigning selected genes to `input$geneselect`
                                selectInput(inputId = "geneselect", label = "Select Genes to Plot", 
                                            choices = as.character(unique(rnamelt$gene)), 
                                            multiple = T, selected = c("HOXA9", "CSF2RA", "KIT")),
                                #This function will show a nicely table defined in the server function as output$table
                                tableOutput("table")
                      ),#End of SidebarPanel
              
                      # Create a mainPanel to display the dynamic plot
                        mainPanel(
                                #this function is needed to display plots created in the server function
                                plotOutput(
                                  # display the plot called `output$mutplot` defined in the server function
                                  outputId = "mutplot", width = "100%", height = "100%")
                        ) #end of mainPanel
              )#end of FluidPage and UI element

#############################################################################################################################################################

# Define server logic required to apply user input to modify plots as `server`

server <- shinyServer(function(input, output) {
  
  # this reactive function outputs our counts (from rnamelt) filtered for genes select by user (input$geneselect)
  # ({ ... }) is the syntax used to define reactive functions. 
  # this data here can be accessed using `genereact()`
  genereact <- reactive({
    filter(rnamelt, gene %in% input$geneselect) %>% 
      # create a new column called `logtag` for log(1+count)
      mutate(logtag=log(1+tagcount))
  })

  # identify patient id that have mutations selected by user; also see if they have >1 mutation (multiple)
  mutid <-  reactive({
    filter(tcga_mut, mutation %in% input$muts) %>% 
    inner_join(nmut()) %>% 
    mutate(status=ifelse(n < 2, as.character(mutation), "multiple")) %>% 
    select(id, status) %>% 
    unique() %>% data.frame()
    })
  
  # create a reactive function that detects changes in `input$muts` and `input$geneselect` and updates accordingly
  mutreact <- reactive({
      select(tcga_mut, -aa_change, -mutation) %>% unique() %>%
      full_join(mutid(), by="id") %>%
      # group patients by `status`...if they don't have one of selected mutations they are "other"
      mutate(mutant=ifelse(is.na(status), yes = "other", no = status)) %>% 
      # select mutant, id columns
      select(mutant, id) %>%
      # join with gene expression values
      inner_join(genereact())
  })
  
  # get the number of patients with the selected mutations (from input$muts)
  nmut <- reactive({
    filter(tcga_mut, mutation %in% input$muts) %>% group_by(id) %>% tally()
  })

  # create the plot using the reactive data `mutreact()`
  output$mutplot <- renderPlot({
    mutreact() %>%
        ggplot(aes(x=mutant, y=logtag))+
        geom_boxplot(position="dodge", width=0.7, fill="honeydew3", alpha=0.3, outlier.colour = "white")+
        geom_point(data=mutreact(), aes(x=mutant, y=logtag, colour=mutant), position="jitter", size=0.75)+
        ylab("Expression log(1+tagcount)")+
        xlab("")+
      theme_bw()+
        facet_wrap(~gene, scales="free")+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      guides(color=F)+
      scale_colour_brewer(palette = "Set2", name="Mutant Status")
    }, width = 900, height = 450)
  
  # create a table for # of patients w/ different mutations
  output$table <- renderTable({
    print <- mutreact() %>% select(id,mutant) %>% unique() %>% 
      group_by(mutant) %>% tally() %>% data.frame()
    colnames(print) <- c("Mutated Gene", "# of patients")
    print})
})

#############################################################################################################################################################

# Run the application 
shinyApp(ui = ui, server = server)



