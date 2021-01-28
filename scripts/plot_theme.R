theme_thesis <- function(){ 
    font <- "Noto Sans Condensed"   #assign font family up front
    
    theme_classic() %+replace%    #replace elements we want to change
    
    theme(
      #text elements
      axis.title = element_text(
                   family = font,
                   size = 12),
      axis.text = element_text(
                   family = font,
                   size = 11),
      axis.text.x = element_text(
                    margin=margin(5, b = 10))
    )
}