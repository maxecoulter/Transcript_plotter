#Code for trancsript_plotter.

#Original transcript code developed by J Grey Munroe (https://cran.r-project.org/web/packages/genemodel/index.html)

# Max Coulter added functions to display interproscan domains and multiple transcripts simultaneously.



genemodel.plot<-function(model, start, bpstop, orientation, xaxis=TRUE,axis_text_size=1,legend_gap=1)
{
  ylim_bottom <- -0.3 - legend_gap
  lwidth <- 0.4
  par(mar=c(1,1,3,1), cex=0.2, lwd=lwidth)#plot coordinates?
  model<-cbind(model[,4], as.data.frame(stringr::str_split_fixed(model$coordinates, "-", 2)))#Splits coordinates into two columns
  colnames(model)<-c("feature", "start", "bpstop")
  model$start<-as.numeric(as.character(model$start));model$bpstop<-as.numeric(as.character(model$bpstop))
  
  length<-bpstop-start
  
  if (orientation=="reverse")
  {
    model$newstart<-bpstop-model$bpstop+1
    model$newstop<-bpstop-model$start+1
    model<-model[which(model$feature!="exon"),]
    model<-model[which(model$feature!="ORF"),]
    model<-model[order(model$newstart),]
    
    plot(1, type="l",axes=F,ann=FALSE, xlim=c(start-.03*length, bpstop), ylim=c(ylim_bottom, .5))
    for (i in 2:nrow(model))
    {
      type<-model$feature[i]
      
      if (type=="coding_region")
      {
        rect(model$newstart[i], 0, model$newstop[i], .2, col = "steelblue3", border="dodgerblue4", lwd=lwidth)
        
      }
      
      if (type=="intron")
      {
        
        mid<-mean(c(model$newstart[i], model$newstop[i]))
        segments(x0=model$newstart[i],y0=.1,x1=mid,y1=.2, lwd=lwidth, col="dodgerblue4")
        segments(x0=model$newstop[i],y0=.1,x1=mid,y1=.2, lwd=lwidth, col="dodgerblue4")
      }
      
      if (type=="3' utr")
      {
        rect(model$newstart[i], 0, model$newstop[i], .2, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        
      }
      if (type=="5' utr")
      {
        rect(model$newstart[i], 0, model$newstop[i], .2, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        
      }
    }
    x<-c(model$newstop[1], model$newstop[1], model$newstart[1], start-.02*length, model$newstart[1])
    y<-c(0,.2,.2,.1,0)
    endtype<-model$feature[1]
    if (endtype=="coding_region") polygon(x,y, border = "dodgerblue4" , col ="steelblue3" , lwd=lwidth)
    else polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=lwidth)
  }
  
  if (orientation=="forward")
  {
    model$newstart<-start+model$start-1
    model$newstop<-start+model$bpstop-1
    model<-model[which(model$feature!="exon"),]
    model<-model[which(model$feature!="ORF"),]
    model<-model[order(model$newstop, decreasing = T),]
    
    plot(1, type="l",yaxt='n',xaxt='n', ann=FALSE, xlim=c(start, bpstop+.03*length), ylim=c(ylim_bottom, .22))
    for (i in 2:nrow(model))
    {
      type<-model$feature[i]
      
      if (type=="coding_region")
      {
        rect(model$newstart[i], 0, model$newstop[i], .2, col = "steelblue3", border="dodgerblue4", lwd=lwidth)
        
      }
      
      if (type=="intron")
      {
        
        mid<-mean(c(model$newstart[i], model$newstop[i]))
        segments(x0=model$newstart[i],y0=.1,x1=mid,y1=.2, lwd=lwidth, col="dodgerblue4")
        segments(x0=model$newstop[i],y0=.1,x1=mid,y1=.2, lwd=lwidth, col="dodgerblue4")
      }
      
      if (type=="3' utr")
      {
        rect(model$newstart[i], 0, model$newstop[i], .2, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        
      }
      if (type=="5' utr")
      {
        rect(model$newstart[i], 0, model$newstop[i], .2, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        
      }
    }
    x<-c(model$newstart[1], model$newstart[1], model$newstop[1], bpstop+.02*length, model$newstop[1])
    y<-c(0,.2,.2,.1,0)
    endtype<-model$feature[1]
    if (endtype=="coding_region") {polygon(x,y, border = "dodgerblue4" , col ="steelblue3" , lwd=lwidth) }
    else {polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=lwidth)}
  }
  if (xaxis==T)
  {
    Axis(side=3, labels=T, cex.axis=axis_text_size, lwd=lwidth)
  }
}

genemodel.plot_domain<-function(model, start, bpstop, orientation,legend_size=2,manual_colours = c()){
  #Model is list of domains with coordinates
  lwidth <- 0.4
  par(mar=c(1,1,3,1), cex=0.2, lwd=lwidth)#plot coordinates?
  colours = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]#Get all 433 R colours
  model<-cbind(model[,4], as.data.frame(stringr::str_split_fixed(model$coordinates, "-", 2)))#Splits coordinates into two columns
  colnames(model)<-c("feature", "start", "bpstop")
  model$size <- as.numeric(as.character(model$bpstop)) - as.numeric(as.character(model$start))
  model$start<-as.numeric(as.character(model$start));model$bpstop<-as.numeric(as.character(model$bpstop))
  all_types <- as.vector(unique(model$feature))
  length<-bpstop-start
  if (length(all_types) == 0){
    stop("No predicted domains!")
  }
  type_colours <- c()
  for (c in 1:length(all_types)){
    type_colours <- c(type_colours, colours[sample(1:length(colours), 1)])
  }
  colours_types <- data.frame("domain"=all_types, "colour"=type_colours)
  #Manual option for domain colours:
  if (length(manual_colours) != 0){
    if (length(manual_colours) != length(colours_types$colour)){
      stop("Incorrect number of colours specified!")
    }
    colours_types$colour <- manual_colours
  }
  rownames(colours_types) <- colours_types$domain
  if (orientation=="reverse"){
    model$newstart<-bpstop-model$bpstop+1
    model$newstop<-bpstop-model$start+1
    model<-model[order(model$size, decreasing = T),]
    #plot(1, type="l",yaxt='n',ann=FALSE, xlim=c(start, bpstop+.03*length), ylim=c(-1, .5))
    #colours_types <- data.frame("domain"=c(), "colour"=c())
    for (i in 1:nrow(model)){
      type<- as.vector(model$feature)[i]
      colour <- as.character(colours_types[type, "colour"])
      #rect(model$newstart[i], 0, model$newstop[i], .2, col = "white", border="dodgerblue4", lwd=1)
      rect(model$newstart[i], 0, model$newstop[i], .2, col = colour, border="dodgerblue4", lwd=lwidth)#
      print(colour)
    }
    
    
  }
  if (orientation=="forward"){
    model$newstart<-start+model$start-1
    model$newstop<-start+model$bpstop-1
    model<-model[order(model$size, decreasing = T),]
    #plot(1, type="l",yaxt='n',ann=FALSE, xlim=c(start, bpstop+.03*length), ylim=c(-1, .5))
    for (i in 1:nrow(model)){
      type<- as.vector(model$feature)[i]
      colour <- as.character(colours_types[type, "colour"])
      rect(model$newstart[i], 0, model$newstop[i], .2, col = colour, border="dodgerblue4", lwd=lwidth)
      print(colour)
      
    }
    
  }
  colours_types <- rbind(colours_types, data.frame("domain"=c("CDS","UTR"), "colour"=c("steelblue3", "lightsteelblue1")))
  legend("bottomright", inset=.02, title="Domains",
         legend=as.vector(colours_types$domain), fill=as.vector(colours_types$colour), horiz=FALSE, cex=legend_size)
  
  colours_types
}


gene.transcript.model.plot<-function(model, gene_start, gene_bpstop, orientation, xaxis=TRUE, gap=0, legend_gap=1,axis_text_size=0.8)
{#Plot all the transcripts from a gene, gap = additional space between transcripts
  #Input: gene, transcript, transcript_coordinates, type, coordinates
  lwidth <- 0.4
  par(mar=c(1,1,3,1), cex=0.2, lwd=lwidth)#plot settings
  transcripts <- unique(as.vector(model$transcript))
  ylim_bottom <- ((-0.3 - gap) * length(transcripts)) - legend_gap
  length<-gene_bpstop-gene_start
  plot(1, type="l",yaxt='n',xaxt='n', ann=FALSE, xlim=c(gene_start, gene_bpstop+.03*length), ylim=c(ylim_bottom, .22), lwd=lwidth)
  plot_bottom <- 0
  plot_top <- 0.2
  
  for (t in 1:length(transcripts)){
    model_transcript <- model %>% filter(transcript == transcripts[t])
    transcript_coordinates <- as.vector(stringr::str_split_fixed(unique(model_transcript$transcript_coordinates), "-", 2))
    start<- as.numeric(transcript_coordinates[1])
    bpstop <- as.numeric(transcript_coordinates[2])
    model_transcript<-cbind(model_transcript[,4], as.data.frame(stringr::str_split_fixed(model_transcript$coordinates, "-", 2)))#Splits coordinates into two columns
    colnames(model_transcript)<-c("feature", "start", "bpstop")
    model_transcript$start<-as.numeric(as.character(model_transcript$start));model_transcript$bpstop<-as.numeric(as.character(model_transcript$bpstop))
    length<-bpstop-start
    
    if (orientation=="forward")
    {
      model_transcript$newstart<-start+model_transcript$start-1
      model_transcript$newstop<-start+model_transcript$bpstop
      if ("coding_region" %in% unique(model_transcript$feature))
      {
        model_transcript<-model_transcript[which(model_transcript$feature!="exon"),]
      }
      model_transcript<-model_transcript[which(model_transcript$feature!="ORF"),]
      model_transcript<-model_transcript[order(model_transcript$newstop, decreasing = T),]
      print(model_transcript$feature)
      if (nrow(model_transcript) > 1)
      {
      for (i in 2:nrow(model_transcript))
      {
        type<-model_transcript$feature[i]
        print(type)
        if (type=="exon" & !("coding_region" %in% unique(model_transcript$feature)))#No CDS transcript
        {
          rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        }
        
        if (type=="coding_region")
        {
          rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "steelblue3", border="dodgerblue4", lwd=lwidth)
          
        }
        
        if (type=="intron")
        {
          
          mid<-mean(c(model_transcript$newstart[i], model_transcript$newstop[i]))
          segments(x0=model_transcript$newstart[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
          segments(x0=model_transcript$newstop[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
        }
        
        if (type=="3' utr")
        {
          rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
          
        }
        if (type=="5' utr")
        {
          rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
          
        }
        
      }
      }
      x<-c(model_transcript$newstart[1], model_transcript$newstart[1], model_transcript$newstop[1], bpstop+.02*length, model_transcript$newstop[1])
      y <- c(plot_bottom, plot_top, plot_top, plot_bottom+0.1, plot_bottom)
      #y<-c(0,.2,.2,.1,0)
      endtype<-model_transcript$feature[1]
      if (endtype=="coding_region") {polygon(x,y, border = "dodgerblue4" , col ="steelblue3" , lwd=lwidth) }
      else {polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=lwidth)}
    }
    
  
  if (orientation=="reverse")
  {
    model_transcript$newstart<-bpstop-model_transcript$bpstop+1
    model_transcript$newstop<-bpstop-model_transcript$start+2
    if ("coding_region" %in% unique(model_transcript$feature))
    {
      model_transcript<-model_transcript[which(model_transcript$feature!="exon"),]
    }
   
    model_transcript<-model_transcript[which(model_transcript$feature!="ORF"),]
    model_transcript<-model_transcript[order(model_transcript$newstart),]
    if (nrow(model_transcript) > 1)
    {
    for (i in 2:nrow(model_transcript))
    {
      type<-model_transcript$feature[i]
      if (type=="exon" & !("coding_region" %in% unique(model_transcript$feature)))#No CDS transcript
      {
        rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
      }
      
      if (type=="coding_region")
      {
        rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "steelblue3", border="dodgerblue4", lwd=lwidth)
        
      }
      
      if (type=="intron")
      {
        
        mid<-mean(c(model_transcript$newstart[i], model_transcript$newstop[i]))
        segments(x0=model_transcript$newstart[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
        segments(x0=model_transcript$newstop[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
      }
      
      if (type=="3' utr")
      {
        rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        
      }
      if (type=="5' utr")
      {
        rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
        
      }
      
    }
    }
    x<-c(model_transcript$newstop[1], model_transcript$newstop[1], model_transcript$newstart[1], start-.02*length, model_transcript$newstart[1])
    y <- c(plot_bottom, plot_top, plot_top, plot_bottom+0.1, plot_bottom)
    endtype<-model_transcript$feature[1]
    if (endtype=="coding_region") polygon(x,y, border = "dodgerblue4" , col ="steelblue3" , lwd=lwidth)
    else polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=lwidth)
  }
    if (xaxis==T)
    {
      Axis(side=3, labels=T, cex.axis=axis_text_size, lwd=lwidth)
    }
    plot_bottom <- plot_bottom - 0.3 - gap
    plot_top <- plot_top - 0.3 - gap
    xaxis=FALSE
    
  }
  
}


gene.transcript.model.plot_domain<-function(model, transcript_vector, gene_start, gene_bpstop, orientation,gap=0.2,legend_gap=1,legend_size=1.3, manual_colours = c()){
  print(length(transcript_vector))
  #Model is list of domains with coordinates
  lwidth <- 0.4
  par(mar=c(1,1,3,1), cex=0.2, lwd=lwidth)
  transcripts <- unique(as.vector(model$transcript))
  print(transcripts)
  ylim_bottom <- ((-0.3 - gap) * length(transcripts)) - legend_gap
  #plot(1, type="l",yaxt='n',ann=FALSE, xlim=c(gene_start, gene_bpstop+.03*length), ylim=c(ylim_bottom, .22))
  plot_bottom <- 0
  plot_top <- 0.2
  
  colours = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]#Get all 433 R colours
  all_types <- as.vector(unique(model$type))
  if (length(all_types) == 0){
    #stop("No predicted domains!")
    return()
  }
  type_colours <- c()
  for (c in 1:length(all_types)){
    type_colours <- c(type_colours, colours[sample(1:length(colours), 1)])
  }
  #print(type_colours)
  colours_types <- data.frame("domain"=all_types, "colour"=type_colours)
  #Manual option for domain colours:
  if (length(manual_colours) != 0){
    if (length(manual_colours) != length(colours_types$colour)){
      stop("Incorrect number of colours specified!")
    }
    colours_types$colour <- manual_colours
  }
  rownames(colours_types) <- colours_types$domain
  for (t in 1:length(transcript_vector)){
    if (!(transcript_vector[t] %in% transcripts)){
      plot_bottom <- plot_bottom - 0.3 - gap
      plot_top <- plot_top - 0.3 - gap
      print("No domains!")
      next#if transcript does not have any domains
    }
    model_transcript <- model %>% filter(transcript == transcript_vector[t])
    transcript_coordinates <- as.vector(stringr::str_split_fixed(unique(model_transcript$transcript_coordinates), "-", 2))
    start<- as.numeric(transcript_coordinates[1])
    bpstop <- as.numeric(transcript_coordinates[2])
    model_transcript<-cbind(model_transcript[,4], as.data.frame(stringr::str_split_fixed(model_transcript$coordinates, "-", 2)))#Splits coordinates into two columns
    colnames(model_transcript)<-c("feature", "start", "bpstop")
    model_transcript$start<-as.numeric(as.character(model_transcript$start));model_transcript$bpstop<-as.numeric(as.character(model_transcript$bpstop))
    length<-bpstop-start
    model_transcript$size <- as.numeric(as.character(model_transcript$bpstop)) - as.numeric(as.character(model_transcript$start))
    #print(transcript_vector[t])
    
  
  if (orientation=="reverse"){
    model_transcript$newstart<-bpstop-model_transcript$bpstop+1
    model_transcript$newstop<-bpstop-model_transcript$start+2
    model_transcript<-model_transcript[order(model_transcript$size, decreasing = T),]
    
    #plot(1, type="l",yaxt='n',ann=FALSE, xlim=c(start, bpstop+.03*length), ylim=c(-1, .5))
    for (i in 1:nrow(model_transcript)){
      type<- as.vector(model_transcript$feature)[i]
      colour <- as.character(colours_types[type, "colour"])
      #rect(model_transcript$newstart[i], 0, model_transcript$newstop[i], .2, col = "white", border="dodgerblue4", lwd=1)
      rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = colour, border="dodgerblue4", lwd=lwidth)#
      
    }
  }
  if (orientation=="forward"){
    model_transcript$newstart<-start+model_transcript$start-1
    model_transcript$newstop<-start+model_transcript$bpstop
    model_transcript<-model_transcript[order(model_transcript$size, decreasing = T),]
    #plot(1, type="l",yaxt='n',ann=FALSE, xlim=c(start, bpstop+.03*length), ylim=c(-1, .5))
    for (i in 1:nrow(model_transcript)){
      type<- as.vector(model_transcript$feature)[i]
      colour <- as.character(colours_types[type, "colour"])
      rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = colour, border="dodgerblue4", lwd=lwidth)
      
      
    }
    
  }
    #Adjust coordinates:
    plot_bottom <- plot_bottom - 0.3 - gap
    plot_top <- plot_top - 0.3 - gap
  }
  colours_types <- rbind(colours_types, data.frame("domain"=c("CDS","UTR"), "colour"=c("steelblue3", "lightsteelblue1")))
  legend("bottomright", inset=.02, title="Domains",
         legend=as.vector(colours_types$domain), fill=as.vector(colours_types$colour), horiz=FALSE, cex=legend_size, box.lwd=lwidth, pt.lwd = rep(0.5,nrow(colours_types)))
  
  colours_types
}

mutation.plot<-function(start, stop, text="", drop=-0.15, col="red", haplotypes=NULL,text_offset=10,lwidth= 0.4)
{
  rect(start, .2, stop, drop+.01*length(haplotypes), col=col, border=col,lwd=lwidth)
  text(stop + text_offset, drop, text, cex=0.7, col=col, pos=4, offset=0.1+.1*length(haplotypes))
  for (i in 1:length(haplotypes)) points(stop, drop-0.005+(i-1)*0.1, col=haplotypes[i], pch=20, cex=2)
}







Transcript.model.browser.plot<-function(model, gene_info, max_transcripts, gene_start, gene_bpstop, xaxis=TRUE, gap=0, legend_gap=0,axis_text_size=1.8)
{#Plot all the transcripts from a gene, gap = additional space between transcripts
  #Input: gene, transcript, transcript_coordinates, type, coordinates
  lwidth <- 0.4
  par(mar=c(1,1,3,1), cex=0.2, lwd=lwidth)#plot settings
  transcripts <- unique(as.vector(model$transcript))
  ylim_bottom <- ((-0.3 - gap) * max_transcripts) - legend_gap
  length<-gene_bpstop-gene_start
  plot(1, type="l",yaxt='n',xaxt='n', ann=FALSE, xlim=c(gene_start, gene_bpstop+.03*length), ylim=c(ylim_bottom, .22), lwd=lwidth)
  plot_bottom <- 0
  plot_top <- 0.2
  previous_gene <- ""
  for (t in 1:length(transcripts)){
    
    
    gene <- stringr::str_split_fixed(transcripts[t],"\\.",2)[1]
    
    orientation <- gene_info[gene,"orientation"]
    
    model_transcript <- model %>% filter(transcript == transcripts[t])
    transcript_coordinates <- as.vector(stringr::str_split_fixed(unique(model_transcript$transcript_coordinates), "-", 2))
    start<- as.numeric(transcript_coordinates[1])
    bpstop <- as.numeric(transcript_coordinates[2])
    model_transcript<-cbind(model_transcript[,4], as.data.frame(stringr::str_split_fixed(model_transcript$coordinates, "-", 2)))#Splits coordinates into two columns
    colnames(model_transcript)<-c("feature", "start", "bpstop")
    model_transcript$start<-as.numeric(as.character(model_transcript$start));model_transcript$bpstop<-as.numeric(as.character(model_transcript$bpstop))
    length<-bpstop-start
    
    if (orientation=="forward")
    {
      model_transcript$newstart<-start+model_transcript$start-1
      model_transcript$newstop<-start+model_transcript$bpstop
      if ("coding_region" %in% unique(model_transcript$feature))
      {
        model_transcript<-model_transcript[which(model_transcript$feature!="exon"),]
      }
      model_transcript<-model_transcript[which(model_transcript$feature!="ORF"),]
      model_transcript<-model_transcript[order(model_transcript$newstop, decreasing = T),]
      print(model_transcript$feature)
      if (nrow(model_transcript) > 1)
      {
        for (i in 2:nrow(model_transcript))
        {
          type<-model_transcript$feature[i]
          #print(type)
          if (type=="exon" & !("coding_region" %in% unique(model_transcript$feature)))#No CDS transcript
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
          }
          
          if (type=="coding_region")
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "steelblue3", border="dodgerblue4", lwd=lwidth)
            
          }
          
          if (type=="intron")
          {
            
            mid<-mean(c(model_transcript$newstart[i], model_transcript$newstop[i]))
            segments(x0=model_transcript$newstart[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
            segments(x0=model_transcript$newstop[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
          }
          
          if (type=="3' utr")
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
            
          }
          if (type=="5' utr")
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
            
          }
          
        }
      }
      x<-c(model_transcript$newstart[1], model_transcript$newstart[1], model_transcript$newstop[1], bpstop+.02*length, model_transcript$newstop[1])
      y <- c(plot_bottom, plot_top, plot_top, plot_bottom+0.1, plot_bottom)
      #y<-c(0,.2,.2,.1,0)
      endtype<-model_transcript$feature[1]
      if (endtype=="coding_region") {polygon(x,y, border = "dodgerblue4" , col ="steelblue3" , lwd=lwidth) }
      else {polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=lwidth)}
    }
    
    
    if (orientation=="reverse")
    {
      model_transcript$newstart<-bpstop-model_transcript$bpstop+1
      model_transcript$newstop<-bpstop-model_transcript$start+2
      if ("coding_region" %in% unique(model_transcript$feature))
      {
        model_transcript<-model_transcript[which(model_transcript$feature!="exon"),]
      }
      
      model_transcript<-model_transcript[which(model_transcript$feature!="ORF"),]
      model_transcript<-model_transcript[order(model_transcript$newstart),]
      if (nrow(model_transcript) > 1)
      {
        for (i in 2:nrow(model_transcript))
        {
          type<-model_transcript$feature[i]
          if (type=="exon" & !("coding_region" %in% unique(model_transcript$feature)))#No CDS transcript
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
          }
          
          if (type=="coding_region")
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "steelblue3", border="dodgerblue4", lwd=lwidth)
            
          }
          
          if (type=="intron")
          {
            
            mid<-mean(c(model_transcript$newstart[i], model_transcript$newstop[i]))
            segments(x0=model_transcript$newstart[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
            segments(x0=model_transcript$newstop[i],y0=plot_bottom + 0.1,x1=mid,y1=plot_top, lwd=lwidth, col="dodgerblue4")
          }
          
          if (type=="3' utr")
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
            
          }
          if (type=="5' utr")
          {
            rect(model_transcript$newstart[i], plot_bottom, model_transcript$newstop[i], plot_top, col = "lightsteelblue1", border="dodgerblue4", lwd=lwidth)
            
          }
          
        }
      }
      x<-c(model_transcript$newstop[1], model_transcript$newstop[1], model_transcript$newstart[1], start-.02*length, model_transcript$newstart[1])
      y <- c(plot_bottom, plot_top, plot_top, plot_bottom+0.1, plot_bottom)
      endtype<-model_transcript$feature[1]
      if (endtype=="coding_region") polygon(x,y, border = "dodgerblue4" , col ="steelblue3" , lwd=lwidth)
      else polygon(x,y, border = "dodgerblue4" , col ="lightsteelblue1" , lwd=lwidth)
    }
    if (xaxis==T)
    {
      Axis(side=3, labels=T, cex.axis=axis_text_size, lwd=lwidth)
    }
    next_transcript <- transcripts[t+1]
    #print("Next transcript:")
    print(next_transcript)
    next_gene <- stringr::str_split_fixed(next_transcript,"\\.",2)[1]
    if (gene == next_gene){
      plot_bottom <- plot_bottom - 0.3 - gap
      plot_top <- plot_top - 0.3 - gap
      xaxis=FALSE
    }else{
      plot_bottom <- 0
      plot_top <- 0.2
    }
    
    
  }
  
}


Transcript.Browser <- function(gene_info,
                               transcripts,
                               plot_start,
                               plot_end)
{
  filtered_genes <- gene_info %>% filter(end > plot_start) %>% filter(start < plot_end)
  
  transcripts_filtered <- transcripts %>% filter(gene %in% filtered_genes$gene)
  max_transcripts<- max(unlist(lapply(transcripts_filtered$transcript,t_no)))
  
  Transcript.model.browser.plot(model=transcripts_filtered, gene_info=filtered_genes,max_transcripts= max_transcripts, gene_start=plot_start, gene_bpstop=plot_end, xaxis=T, gap=0.2)
  
}

# Transcript.Browser <- function(){
#   
#   model
# }


