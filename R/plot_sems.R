library(DiagrammeR)
library(piecewiseSEM)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)

########
#this function modifies a coefs_df from sem to plot the sem later
coefsdf_to_grViz <- function(coefs.df, mode = c('single','with_CI','pseudoposterior')){
  ##grViz block
  #replace dots in variable names
  coefs.df$Response<-gsub(coefs.df$Response,pattern='\\.',replacement='_')
  coefs.df$Predictor<-gsub(coefs.df$Predictor,pattern='\\.',replacement='_')
  #this controls the width of the arrows
  coefs.df$Estimate <- as.numeric(coefs.df$Estimate)
  coefs.df$edge_width<-abs(coefs.df$Estimate)*20
  coefs.df$edge_number<-coefs.df$Estimate
  if('Std.Error' %in% colnames(coefs.df)){
    coefs.df$edge_error<-coefs.df$Std.Error
  }
  coefs.df$edge_significance<-coefs.df$P.Value
  coefs.df$edge_colour<-NA
  #assign colours to positive and negative effects
  coefs.df[coefs.df$edge_number>0,'edge_colour']<-'tomato'
  coefs.df[coefs.df$edge_number<0,'edge_colour']<-'steelblue'
  #change colour of number
  #coefs.df$edge_font<-'black'
  coefs.df$edge_font<-NA
  coefs.df[coefs.df$edge_number>0,'edge_font']<-'tomato'
  coefs.df[coefs.df$edge_number<0,'edge_font']<-'steelblue'
  #use grey for arrows and font of non-significant edges
  coefs.df[coefs.df$edge_significance>0.05,'edge_font']<-'grey'
  coefs.df[coefs.df$edge_significance>0.05,'edge_colour']<-'grey'
  
  #create vectors to be passed on to
  edge_width<-coefs.df$edge_width
  edge_number<-round(coefs.df$edge_number,3)
  if('Std.Error' %in% colnames(coefs.df)){
    edge_error<-round(coefs.df$edge_error,3)
  }
  edge_significance<-coefs.df$P.Value
  response<-coefs.df$Response
  predictor<-coefs.df$Predictor
  #changing variable names here
  response<-gsub(response,pattern="mean_",replacement="")
  response<-gsub(response,pattern="elevation_ETOPO_land",replacement="present_elevation")
  response<-gsub(response,pattern="mean_present_T",replacement="present_T")
  response<-gsub(response,pattern="present_minus_past",replacement="change")
  response<-gsub(response,pattern="temperature",replacement="T")
  response<-gsub(response,pattern="elevation_gain",replacement="gain_elevation")
  response<-gsub(response,pattern="elevation_loss",replacement="loss_elevation")
  
  predictor<-gsub(predictor,pattern="mean_",replacement="")
  predictor<-gsub(predictor,pattern="elevation_ETOPO_land",replacement="present_elevation")
  predictor<-gsub(predictor,pattern="mean_present_T",replacement="present_T")
  predictor<-gsub(predictor,pattern="present_minus_past",replacement="change")
  predictor<-gsub(predictor,pattern="temperature",replacement="T")
  predictor<-gsub(predictor,pattern="elevation_gain",replacement="gain_elevation")
  predictor<-gsub(predictor,pattern="elevation_loss",replacement="loss_elevation")
  
  edge_colour<-coefs.df$edge_colour
  edge_font<-coefs.df$edge_font
  var_names<-unique(c(response,predictor))
  if(mode == 'single'){
    edge_label <- edge_number
  }else if (mode =='with_CI'){
    edge_label <- paste0(edge_number, 
                         ' (', 
                         round(edge_number-(1.96*edge_error),3),
                         ' - ',
                         round(edge_number+(1.96*edge_error),3), 
                         ')')
  }else if (mode == 'pseudoposterior'){
    edge_label <- edge_number 
    for (i in c(1:length(coefs.df$CIlow))){
      if(coefs.df$CIup[i] != coefs.df$CIlow[i]){
        edge_label[i] = paste0(edge_label[i],
                               ' (', 
                               round(coefs.df$CIlow[i],3),
                               ' - ',
                               round(coefs.df$CIup[i],3), 
                               ')')
      }
    }
  }
  
  grViz_input <- list(edge_width,
                      edge_number,
                      edge_error,
                      response,
                      predictor,
                      edge_colour,
                      edge_font,
                      var_names,
                      edge_label)
  names(grViz_input) <- c('edge_width',
                          'edge_number',
                          'edge_error',
                          'response',
                          'predictor',
                          'edge_colour',
                          'edge_font',
                          'var_names',
                          'edge_label')
  
  list2env(grViz_input, envir = .GlobalEnv)
  grViz <- ("
                    digraph boxes_and_circles {
                      #   a 'graph' statement
                      
                      graph [overlap = true, fontsize = 10]
                      node [shape = box,
                      fontname = Helvetica]
                      
                      subgraph {
                      rank = same; '@@7-5'; '@@7-4'
                      }
                      subgraph {
                      rank = same;'@@7-2';'@@7-3'
                      }
                      # several 'node' statements
                      '@@7-1'; '@@7-2'; '@@7-3'; '@@7-4'; '@@7-5'
                      # several 'edge' statements
                      '@@4-1'->'@@3-1'[penwidth='@@1-1',label='@@8-1',color='@@5-1',fontcolor='@@6-1']
                      '@@4-2'->'@@3-2'[penwidth='@@1-2',label='@@8-2',color='@@5-2',fontcolor='@@6-2']
                      '@@4-3'->'@@3-3'[penwidth='@@1-3',label='@@8-3',color='@@5-3',fontcolor='@@6-3']
                      '@@4-4'->'@@3-4'[penwidth='@@1-4',label='@@8-4',color='@@5-4',fontcolor='@@6-4']
                      '@@4-5'->'@@3-5'[penwidth='@@1-5',label='@@8-5',color='@@5-5',fontcolor='@@6-5']
                      '@@4-6'->'@@3-6'[penwidth='@@1-6',label='@@8-6',color='@@5-6',fontcolor='@@6-6']
                      '@@4-7'->'@@3-7'[penwidth='@@1-7',label='@@8-7',color='@@5-7',fontcolor='@@6-7']
                      '@@4-8'->'@@3-8'[penwidth='@@1-8',label='@@8-8',color='@@5-8',fontcolor='@@6-8']
                      '@@4-9'->'@@3-9'[penwidth='@@1-9',label='@@8-9',color='@@5-9',fontcolor='@@6-9']
                       }
                       [1]: edge_width
                       [2]: edge_number
                       [3]: response
                       [4]: predictor
                       [5]: edge_colour
                       [6]: edge_font
                       [7]: var_names
                       [8]: edge_label
                       ")
  
  return(grViz)
  
}
########