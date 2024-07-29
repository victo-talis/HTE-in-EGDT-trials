
# Original code modified from treeshap package
# Much of this code has been modified slightly further to accept output from shapr

plot_feature_importance<-function(shap,
                                  desc_sorting = TRUE,
                                  max_vars = ncol(shaps),
                                  title = "",
                                  subtitle = NULL,
                                  custom_var_names = NULL,
                                  tsize=12,psize=2,ljust=c(0.5,0),
                                  BREAKS=c(0,0.5,1.0,1.5)) {
  ljustx<-ljust[1]
  ljusty<-ljust[2]
  
  shaps <- shap
  
  # argument check
  # if (!is.treeshap(treeshap)) {
  #   stop("treeshap parameter has to be correct object of class treeshap. Produce it using treeshap function.")
  # }
  
  if (!is.logical(desc_sorting)) {
    stop("desc_sorting is not logical.")
  }
  
  if (!is.numeric(max_vars)) {
    stop("max_vars is not numeric.")
  }
  
  if (max_vars > ncol(shaps)) {
    warning("max_vars exceeded number of explained variables. All variables will be shown.")
    max_vars <- ncol(shaps)
  }
  
  mean <- colMeans(abs(shaps))
  
  if(!is.null(custom_var_names)){
    df <- data.frame(variable = factor(custom_var_names), 
                     importance = as.vector(mean))
  } else {
    df <- data.frame(variable = factor(names(mean)), 
                     importance = as.vector(mean))
  }
  df$variable <- reorder(df$variable, df$importance * ifelse(desc_sorting, 1, -1))
  df <- df[order(df$importance, decreasing = TRUE)[1:max_vars], ]

  p <- ggplot(df, aes(x = variable, y = importance, color=importance)) +
    geom_bar(stat = "identity", fill = "#000000") +
    scale_color_gradientn(colours = c("white"),
                          oob=squish,
                          limits=c(-3,3),
                          labels=c("","","","","","",""),
                          guide=guide_colorbar(barwidth = 15,
                                               barheight = 0.3,
                                               ticks=FALSE,
                                               title.position="top"))
  
  
  p + coord_flip() +
    labs(title = title, subtitle = subtitle, colour="") +
    #theme_drwhy_vertical() +
    ylab("Variable importance") + xlab("") +
    labs(title = title, subtitle = subtitle) +
    scale_y_continuous(labels = number_format(accuracy=0.1),
                       breaks= BREAKS) +
    #theme(legend.position = "none") +
    theme(legend.position = "top",
          legend.direction="horizontal",
          legend.title.align=ljustx,
          legend.justification = c(ljustx,ljusty),
          legend.margin=ggplot2::margin(-5,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10),
          legend.title=element_text(size=tsize),#,family="Arial"),
          axis.text=element_text(size=tsize),#,family="Arial"),
          legend.text=element_text(size=tsize),#, family="Arial"),
          axis.title=element_text(size=tsize),#,family="Arial"),
          plot.margin=ggplot2::margin(0,5,0,0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    )
}



plot_contribution<-function(shap,
                            x_test,
                            #mean_pred,
                            obs = 1,
                            max_vars = 5,
                            min_max = NA,
                            digits = 3,
                            manual_yaxis_vals=NULL,
                            explain_deviation = FALSE,
                            title = "SHAP Break-Down",
                            subtitle = "",
                            axtxtsz,
                            axttlsz,
                            pltttlsz) {
  #browser()
  pred_none<-as.numeric(shap[obs,1])
  shap <- shap[obs,-1]
  pred <- rowSums(shap)+pred_none
  x <- x_test[obs, ]

  # argument check
  # if (!is.treeshap(treeshap)) {
  #   stop("treeshap parameter has to be correct object of class treeshap. Produce it using treeshap function.")
  # }
  
  if (max_vars > ncol(shap)) {
    warning("max_vars exceeds number of variables. All variables will be shown.")
    max_vars <- ncol(shap)
  }
  if (nrow(shap) != 1) {
    warning("Only 1 observation can be plotted. Plotting 1st one.")
    shap <- shap[1, ]
  }
  
  # setting intercept
  #browser()
  # mean_prediction <- mean(predict.model_unified(treeshap$unified_model, treeshap$unified_model$data))
  # if (explain_deviation) {
  #   mean_prediction <- 0
  # }
  mean_prediction<-mean(pred)
  
  #browser()
  df <- data.frame(variable = colnames(shap), 
                   contribution = as.numeric(shap))
  #browser()
  # setting variable names to showing their value
  df$variable <- paste0(df$variable, " = ", as.character(round(x,digits)))
  
  #browser()
  # selecting max_vars most important variables
  is_important <- order(abs(df$contribution), decreasing = TRUE)[1:max_vars]
  other_variables_contribution_sum <- sum(df$contribution[-is_important])
  df <- df[is_important, ]
  df$position <- 2:(max_vars + 1)
  if (max_vars < ncol(shap)) {
    df <- rbind(df, data.frame(variable = "+ All others",
                               contribution = other_variables_contribution_sum,
                               position = max(df$position) + 1))
  }
  #browser()
  # adding "prediction" bar
  df <- rbind(df, data.frame(variable = ifelse(explain_deviation, "prediction deviation", "prediction"),
                             contribution = pred,
                             position = max(df$position) + 1))
  
  df$sign <- ifelse(df$contribution >= 0, "1", "0")
  
  # adding "intercept" bar
  df <- rbind(df, data.frame(variable = "Intercept",
                             contribution = pred_none,
                             position = 1,
                             sign = 2))
  
  # ordering
  df <- df[order(df$position), ]
  
  # adding columns needed by plot
  df$cumulative <- cumsum(df$contribution)
  df$prev <- df$cumulative - df$contribution
  df$text <- as.character(round(df$contribution, digits))
  df$text[df$contribution > 0] <- paste0("+", df$text[df$contribution > 0])
  
  # intercept bar corrections:
  df$prev[1] <- df$contribution[1]
  df$text[1] <- as.character(round(df$contribution[1], digits))
  
  # prediction bar corrections:
  df$prev[nrow(df)] <- df$contribution[1]
  df$cumulative[nrow(df)] <- df$cumulative[max_vars + 2]
  if (!explain_deviation) { #  assuring it doesn't differ from prediction because of some numeric errors
    df$cumulative[nrow(df)] <- pred
  }
  df$sign[nrow(df)] <- 2
  df$text[nrow(df)] <- as.character(round(df$contribution[nrow(df)], digits))
  
  # removing intercept bar if requested by explain_deviation argument
  if (explain_deviation) {
    df <- df[-1, ]
  }
  
  # reversing postions to sort bars decreasing
  df$position <- rev(df$position)
  
  # VBT: Manual entry of Y-axis values
  if(!is.null(manual_yaxis_vals)){
    df$variable<-manual_yaxis_vals
  }
  
  #df$sign<-factor(df$sign,ordered=TRUE,levels=0:2)
  
  if(all(as.numeric(df$sign[-c(1,nrow(df))])==0)){
    colors<-c("#F8766D", "#619CFF")
  } else if (all(as.numeric(df$sign[-c(1,nrow(df))])==1)){
    colors<-c("#00BA38", "#619CFF")
  } else {
    colors<-c("#F8766D", "#00BA38", "#619CFF")
  }
  
  
  #browser()
  # base plot
  p <- ggplot(df, aes(x = position + 0.5,
                      y = pmax(cumulative, prev),
                      xmin = position + 0.15, 
                      xmax = position + 0.85,
                      ymin = cumulative, 
                      ymax = prev,
                      fill = sign,
                      label = text))
    
  #browser()
  # add rectangles and hline
  p <- p +
    geom_errorbarh(data = df[-c(nrow(df), if (explain_deviation) nrow(df) - 1), ],
                   aes(xmax = position - 0.85,
                       xmin = position + 0.85,
                       y = cumulative), height = 0,
                   color = "#371ea3") +
    geom_rect(alpha = 0.9) +
    if (!explain_deviation) (geom_hline(data = df[df$variable == "intercept", ],
                                        aes(yintercept = contribution),
                                        lty = 3, alpha = 0.5, color = "#371ea3"))
  
  #browser()
  # add adnotations
  drange <- diff(range(df$cumulative))
  p <- p + geom_text(aes(y = pmax(cumulative,  cumulative - contribution)),
                     vjust = 0.5,
                     nudge_y = drange * 0.05,
                     hjust = 0,
                     color = "#371ea3") 
  
  # set limits for contributions
  if (any(is.na(min_max))) {
    x_limits <- scale_y_continuous(expand = c(0.05, 0.15), name = "SHAP value", labels = scales::comma)
  } else {
    x_limits <- scale_y_continuous(name = "SHAP value", limits = min_max, labels = scales::comma)
  }
  #browser()
  p <- p + x_limits +
    scale_x_continuous(labels = df$variable, breaks = df$position + 0.5, name = "") # +
    #scale_fill_manual(values = colors_breakdown_drwhy())
  
  #browser()
  # add theme
  p + coord_flip() +
    ylab("SHAP value") + #theme_drwhy_vertical() +
    theme(legend.position = "none",
          text=element_text(size=axtxtsz, family="Arial"),
          axis.text=element_text(size=axtxtsz, family="Arial"),
          axis.title=element_text(size=axttlsz, family="Arial"),
          plot.title=element_text(size=pltttlsz, face="bold",family="Arial",hjust = -0.2),
          plot.margin=ggplot2::margin(4,5,0,0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_fill_manual(values=colors)
}




plot_feature_beeswarm2<-function(shap, 
                                 x_test,
                                 title = "Feature Importance", 
                                 subtitle = NULL,
                                 custom_var_names=NULL,
                                 width=1,
                                 wh_bin=NULL,
                                 xlim=c(-1,1),
                                 tsize=10,
                                 psize=1,
                                 ljust,
                                 max_vars){

  xlim_min<-xlim[1]
  xlim_max<-xlim[2]
  ljustx<-ljust[1]
  ljusty<-ljust[2]
  
  # shaps <- treeshap$shaps[,xvars]
  # feat<-treeshap$observations[,xvars]
  # if (!is.treeshap(treeshap)) {
  #   stop("treeshap parameter has to be correct object of class treeshap. Produce it using treeshap function.")
  # }
  shaps<-as.data.frame(shap)
  feat<-as.data.frame(x_test)
  xvars<-names(feat)
  if(is.null(custom_var_names)){
    plot_names<-names(xvars)
  } else {
    plot_names<-custom_var_names
  }

  mean <- colMeans(abs(shaps))
  ord<-order(mean,decreasing=TRUE)
  if(is.null(max_vars)){
    ord2<-ord
  } else {
    ord2<-ord[1:max_vars]
  }
  
  wh_bin2<-wh_bin[which(wh_bin %in% ord2)] #indices of bin vars in set to be plotted.
  xvars2<-xvars[ord2]
  shaps2<-shaps[,ord2]
  feat2<-feat[,ord2]
  mean2<-mean[ord2]
  plot_names2<-plot_names[ord2]

  #try to trim the extremes to get a better color gradient
  
  feat3<-feat2
  for(iii in 1:ncol(feat3)){
    if(iii %in% wh_bin2){
      feat3[feat3[,iii]<0.5,iii]<-(-10)
      feat3[feat3[,iii]>=0.5,iii]<-10
    } else {
      feat3[,iii]<-(feat3[,iii]-mean(feat3[,iii]))/sd(feat3[,iii])
    }
  }
  
  df<-data.frame(variable=paste0(plot_names2[1]),
                 shap=shaps2[,1],
                 Feature=feat3[,1])
  
  for(j in 2:ncol(shaps2)){
    df0<-data.frame(variable=paste0(plot_names2[j]),
                    shap=shaps2[,j],
                    Feature=feat3[,j])
    df<-rbind(df,df0)
  }
  
  varnames<-rev(plot_names2)
  
  # Add the sum of the other variables
  if(!is.null(max_vars)){
    if(max_vars != ncol(shaps)){
      shaps0<-shaps[,-ord2]
      shaps_other<-apply(shaps0,1,sum)
      df0<-data.frame(variable="Other variables",
                      shap=shaps_other,
                      Feature=10)
      df<-rbind(df,df0)
      
      varnames<-rev(c(plot_names2,"Other variables"))
    }
  }
  
  #browser()
  
  df$variable<-factor(df$variable,
                      levels=varnames)

  p <- ggplot(df, aes(x = variable, y = shap, color=Feature)) + geom_quasirandom(size=psize)
  
  p + scale_color_gradientn(colours = c("red","red","green","darkblue","darkblue"),
                            values=rescale(x=c(-10,-3,-2,2,3,10)),
                            oob=squish,
                            limits=c(-3,3),
                            labels=c("Low/Absent","","","","","","High/Present"),
                            guide=guide_colorbar(barwidth = 15, 
                                                 barheight = 0.3,
                                                 ticks=FALSE,
                                                 title.position="top")) +
    #scale_color_brewer() +
    coord_flip() + #theme_drwhy_vertical() + 
    ylab("SHAP value") + 
    xlab("") + 
    labs(title = title, subtitle = subtitle, colour="Covariate Value") + 
    #scale_y_discrete(limits=c("")) + 
    theme(legend.position = "top",
          legend.direction="horizontal",
          legend.title.align=ljustx,
          legend.justification = c(ljustx,ljusty),
          legend.margin=ggplot2::margin(-5,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10),
          legend.title=element_text(size=tsize),#,family="Arial"),
          axis.text=element_text(size=tsize),#,family="Arial"),
          legend.text=element_text(size=tsize),#, family="Arial"),
          axis.title=element_text(size=tsize),#,family="Arial"),
          plot.margin=ggplot2::margin(0,5,0,0),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    ) +
    ylim(xlim_min,xlim_max)
}
