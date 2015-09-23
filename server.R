library(shiny)
library(sqldf)



# Define server logic required to summarize and view the 
# selected dataset

# Define server logic required to summarize and view the 
# selected dataset
shinyServer(function(input, output) {

  datain<-reactive({
    input$topic
           })

  
  output$summary<-renderTable({

    if (datain()=='ALL'){
      subset<-x4
    }
    else{
    subset<-subX[subX$ID==datain(),]
    }
    
    subset2<-subset[order(subset$pvalue),]                         
    subset2[c(1:5),]
})
  output$histogram<-renderPlot({
    if (datain()=='ALL'){
      subset<-x4
    }
    else{
      subset<-subX[subX$ID==datain(),]
    }
    
    subset2<-subset[order(-subset$pvalue),]
    hist(-log(subset2$pvalue), breaks=100, xlab="-log10(p-value)", main="Mouse phenotypes to human disease test results",col="blue")
  
})

  
})
