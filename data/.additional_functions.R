
# Additional Scripts -------------------------------------------------------------------


progress <- function(max=3) {
  
  phrases <-c(rep(c("Tom is Cool!", 
                    "Spatial Transcriptomics is fun!", 
                    "Single cell RNA seq is very expensive!", 
                    "IMC stands for imaging massive cytometry!", 
                    "Public data is free!",
                    " :^) ", 
                    " :^( ", 
                    " :^D ", 
                    " >:^) ",
                    
                    "Nothing is true, and everything is possible", 
                    "The mitochondria is the powerhouse of the cell"), each=3),
              "Did you know that this is a joke?",
              "Did you know that this progress bar is fake?"
  )
  
  cat(bold(green("\n")))
  
  a <- ceiling(runif(1, min=0, max=length(phrases)))
  line<-paste(white("\r|--------|\n"),yellow(bold("Fun Fact:",phrases[a])))
  cat(line)
  phrases <- phrases[-a]
  Sys.sleep(ceiling(runif(1, min = 1, max=max)))
  flush.console()
  cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
  
  a <- ceiling(runif(1, min=0, max=length(phrases)-1))
  line<-paste(white("\r|##------|\n"),yellow(bold("Fun Fact:",phrases[a])))
  cat(line)
  phrases <- phrases[-a]
  Sys.sleep(ceiling(runif(1, min = 1, max=max)))
  flush.console()
  cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
  
  
  a <- ceiling(runif(1, min=0, max=length(phrases)-2))
  line<-paste(white("\r|####----|\n"),yellow(bold("Fun Fact:",phrases[a])))
  cat(line)
  phrases <- phrases[-a]
  Sys.sleep(ceiling(runif(1, min = 1, max=max)))
  flush.console()
  cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
  
  a <- ceiling(runif(1, min=0, max=length(phrases)-3))
  line<-paste(white("\r|######--|\n"),yellow(bold("Fun Fact:",phrases[a])))
  cat(line)
  phrases <- phrases[-a]
  Sys.sleep(ceiling(runif(1, min = 1, max=max)))
  flush.console()
  cat("\r", paste(rep(" ", nchar(line)), collapse=""), sep="")
  
  cat(white("\r|########|\n\n"),bold(yellow("Done!\n")))
}

cycle <- function(seconds=5) {
  cat(yellow("\nLoading:\n"))
  x=0
  while(x<seconds){
    Sys.sleep(0.1);x=x+0.1
    cat(white("\r   |"))
    Sys.sleep(0.1);x=x+0.1
    cat(white("\r   /"))
    Sys.sleep(0.1);x=x+0.1
    cat(white("\r   \u2500"))
    Sys.sleep(0.1);x=x+0.1
    cat(white("\r   \\"))
    Sys.sleep(0.1);x=x+0.1
  }
  cat(white("\n\n"))
  
}
