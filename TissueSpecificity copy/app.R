#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(ggplot2)
library(ggpubr)
library(scales)
library(viridis)
library(plyr)
library(dplyr)
library(tidyr)
library(cowplot)
library(rstatix)
library(tidyverse)
library(RColorBrewer)
library(data.table)
library(stringr)
library(ps)

#read in file with Age
#Tissue
#Expression
#isTestisPeak
#isOvaryPeak
#isS2Peak
#nTF
#Tf activity in Tissue
#hide errors
#Put citation at bottom, also for flyatlas

#then, for lists of input genes, convert symbol to fbgn, calculate: propexpressed, propmax, tau
#make ggplots with option to facet by Age, istestisPeak, isOvaryPeak, isS2Peak, isS2+ovary-testis, S2+testis, ovary+testis, just ovary, just testis, just s2




flydata<-read.table("201030_flydata_gene_tissue_age_exp_tau_maxtissue.txt")
flydata$Tissue<-gsub(flydata$Tissue, pattern = "Spermatheca", replacement=" Spermatheca")
flydata2<-flydata[!is.na(flydata$Age),]
genestsv<-read.table("genes.tsv")
names(genestsv)<-c("ID", "Gene")


# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(tags$style(".shiny-output-error{color: white;}")),
    # Application title
    titlePanel("Drosophila tissue specificity explorer "),

    # Sidebar with a box for gene input and checkbox for faceting plots by gene age
    sidebarLayout(
        sidebarPanel(
          textAreaInput("geneInput", "Enter a list of FBgn ids",value = "FBgn0000008 FBgn0000014 FBgn0000015 FBgn0000018 FBgn0000022 FBgn0000024"),
          radioButtons("radio", "split by", choices=list("Tissue", "Gene Age", "Neither"), selected="Neither"),
          checkboxGroupInput("checkGroup", "Tissues to display:", choices=unique(flydata$Tissue), selected=unique(flydata$Tissue)),
          downloadButton(outputId = 'downloadData', label = "Download source data")
    ),

        # Show a plot of the generated distribution
    mainPanel(
      p("This app is a companion to Witt, Svetec, Benjamin and Zhao 2020. It uses RNA-seq data from Flyatlas2: Leader, D. P., Krause, S. A. , Pandit, A., Davies, S. A. and Dow, J. A.T. (2018) FlyAtlas 2: a new version of the Drosophila melanogaster expression atlas with RNA-Seq, miRNA-Seq and sex-specific data. Nucleic Acids Research 46 D809–D815. "),
        tags$a(href="http://flyatlas.gla.ac.uk/FlyAtlas2/index.html?page=help", "Flyatlas2"),
        tags$br(),
        tags$a(href="https://zhaolab.rockefeller.edu/", "Zhao lab"),  
        plotOutput("coolplot2", width="100%"),
        plotOutput("coolplot3", width="100%", height="600px"),
        plotOutput("coolplot4", width="100%",height="800px")
        )
    )
)

server <- function(input, output) {
  output$downloadData <- downloadHandler(contentType = "text/csv",
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(flydata2,file)
    }
  )
  

  
  state <- reactiveValues()
 # If name is FBgn, convert to gene symbol
  observe({
   state$x<-input$geneInput
    if( state$x %in% genestsv$ID){
      state$y<- subset(genestsv, ID==input$geneInput )$Gene
    }else{
      state$y<-input$geneInput
      #state$y is the gene symbol for plotting
    }
    
  
  })
  #output$downloadReport <- downloadHandler(
   # filename = 'graphs.pdf',
  #  content = function(file) {
  #   pdf(file)
   #   print()
    #    dev.off()
    #})
  

  output$reactive_choices<-reactive({
    str_extract_all(input$geneInput, genestsv$ID)
    
    
  })
  
  
  

  #plot gene expression, faceted by gene age if box is checked
 observe(   
  output$coolplot2<-renderPlot(
    
      if (input$radio== "Neither"){
              ggplot(
               subset(
                  unique(flydata[,c("gene", "Tissue","logFPKM")]), gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup),
                aes(x=Tissue, y=logFPKM))+
                geom_violin(aes(color=Tissue, fill=Tissue))+geom_boxplot(aes(color=Tissue), fill="white" , width=0.3)+
                  theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+scale_fill_viridis_d()+scale_color_viridis_d()
     }
     else if (input$radio== "Gene Age")
       {
         ggplot(
           subset(
             unique(flydata2[,c("gene", "Tissue","logFPKM", "Age")]), gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup),
           aes(x=Tissue, y=logFPKM))+
           geom_violin(aes(color=Tissue, fill=Tissue))+geom_boxplot(aes(color=Tissue), fill="white" , width=0.3)+
         theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+
          # geom_jitter(alpha=input$sliderAlpha, size=0.5)+
           theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+facet_wrap(~Age)+scale_fill_viridis_d()+scale_color_viridis_d()
     }  
     else {
       ggplot(
         subset( unique(flydata2[,c("gene", "Tissue","logFPKM", "Age")]), gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup),
         aes(x=Age, y=logFPKM))+
         geom_violin(aes(color=Age, fill=Age))+geom_boxplot(aes(color=Age), fill="white", width=0.3 )+
         #geom_jitter(alpha=input$sliderAlpha, size=0.5)+
         theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+facet_wrap(~Tissue)+scale_fill_viridis_d()+scale_color_viridis_d()
       
       }#+scale_color_gradientn(limits=c(input$Featureplotslider[1], input$Featureplotslider[2]),oob = scales::squish,colors=c(input$minColor, input$maxColor))+ggtitle(state$y)+theme(text=element_text(size=15),plot.title=element_text(size=40, face="italic"))})
  )
 )
    #compare tau vs reference distribution
    output$coolplot3<-renderPlot(
       if (input$radio=="Gene Age")
      {
        tmp1<-data.frame(subset(
          unique(flydata2[,c("gene", "tau", "Age")]), gene %in% str_extract_all(input$geneInput, genestsv$ID)))
        tmp1$group="User-selected genes"
        
        tmp2<-data.frame(subset(
          unique(flydata2[,c("gene", "tau","Age" )]), !gene %in% str_extract_all(input$geneInput, genestsv$ID)))
        tmp2$group="other genes"
        tmp3<-rbind(tmp1, tmp2)
        ggplot(tmp3,
               aes(x=group, y=tau))+geom_violin(aes(color=group, fill=group))+geom_boxplot(aes(color=group), fill="white", width=0.2 )+#geom_jitter(alpha=input$sliderAlpha, size=0.5)+
          theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+stat_compare_means(label.y=1.05)+facet_wrap(~Age)+scale_fill_viridis_d()+scale_color_viridis_d()
      } else {
        
          
          tmp1<-data.frame(subset(
            unique(flydata[,c("gene", "tau")]), gene %in% str_extract_all(input$geneInput, genestsv$ID)))
          tmp1$group="User-selected genes"
          
          tmp2<-data.frame(subset(
            unique(flydata[,c("gene", "tau" )]), !gene %in% str_extract_all(input$geneInput, genestsv$ID)))
          tmp2$group="other genes"
          tmp3<-rbind(tmp1, tmp2)
          ggplot(tmp3,
                 aes(x=group, y=tau))+geom_point()+geom_violin(aes(color=group, fill=group))+geom_boxplot(aes(color=group), fill="white", width=0.2 )+#geom_jitter(alpha=input$sliderAlpha, size=0.5)+
            theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+stat_compare_means(label.y=1.05)+scale_fill_viridis_d()+scale_color_viridis_d()
    })
    
      
    #plot TF activity
    output$coolplot4<-renderPlot(
      if (input$radio=="Gene Age")
      {
        tmp1<-data.frame(subset(
          unique(flydata2[,c("gene", "connectivity", "Age", "Tissue")]), gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup))
        tmp1$group="User-selected genes"
        
        tmp2<-data.frame(subset(
          unique(flydata2[,c("gene", "connectivity","Age", "Tissue" )]), !gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup))
        tmp2$group="other genes"
        tmp3<-rbind(tmp1, tmp2)
        ggplot(tmp3,
               aes(x=group, y=connectivity))+ylab("TF expression of genes in tissue")+ylab("Upstream TF expression of gene in tissue")+geom_violin(aes(color=group, fill=group, group=Age))+geom_boxplot(aes(color=group), fill="white", width=0.2 )+#geom_jitter(alpha=input$sliderAlpha, size=0.5)+
          theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+stat_compare_means(label.y=30)+ylim(0,35)+facet_wrap(Age~Tissue)+scale_fill_viridis_d()+scale_color_viridis_d()
      } else  {
        tmp1<-data.frame(subset(
          unique(flydata2[,c("gene", "connectivity", "Tissue")]), gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup))
        tmp1$group="User-selected genes"
        
        tmp2<-data.frame(subset(
          unique(flydata2[,c("gene", "connectivity","Tissue" )]), !gene %in% str_extract_all(input$geneInput, genestsv$ID) & Tissue %in% input$checkGroup))
        tmp2$group="other genes"
        tmp3<-rbind(tmp1, tmp2)
        ggplot(tmp3,
               aes(x=group, y=connectivity))+geom_violin(aes(color=group, fill=group))+ylab("Upstream TF expression of gene in tissue")+geom_boxplot(aes(color=group), fill="white", width=0.2 )+#geom_jitter(alpha=input$sliderAlpha, size=0.5)+
          theme_classic()+theme(legend.position="none", axis.text=element_text(angle=90))+stat_compare_means(label.y=30)+ylim(0,35)+facet_wrap(~Tissue)+scale_fill_viridis_d()+scale_color_viridis_d()
    })
  #  observe(
  #  { if 
  #    (state$y %in% genestsv$Gene)
  #    {
  #      output$coolplot2 <- renderPlot(
  #        {ggplot(
  #          subset(
  #            flydata, gene %in% geneInput2),
  #          aes(x=Tissue, y=logFPKM))+
  #          geom_point()+
  #          geom_violin(aes(color=Tissue))+
  #            theme_classic()#+scale_color_gradientn(limits=c(input$Featureplotslider[1], input$Featureplotslider[2]),oob = scales::squish,colors=c(input$minColor, input$maxColor))+ggtitle(state$y)+theme(text=element_text(size=15),plot.title=element_text(size=40, face="italic"))})
  #        }
  #      )
  #    }
  #  }
  #)
   

}
# Run the application 
shinyApp(ui = ui, server = server)
