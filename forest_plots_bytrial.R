
library(forestplot)
library(dplyr)

machine<-"C:/Users/VIT13/OneDrive - University of Pittsburgh/"
src.dir<-"CATES in Critical Care/EGDT manuscript 1/Summer 2023/For github/"
qmat_c2a<-readRDS(paste0(machine,src.dir,"quantile_mat_c2a.RDS"))
qmat_a2c<-readRDS(paste0(machine,src.dir,"quantile_mat_a2c.RDS"))


# Validation of ProCESS model

base_data <- tibble::tibble(mean  = round(qmat_c2a$arr,3)*100,
                            lower = round(qmat_c2a$arr.hi,3)*100,
                            upper = round(qmat_c2a$arr.lo,3)*100,
                            quintile = paste0("Q",1:5),
                            y_egdt = as.character(qmat_c2a$x1),
                            n_egdt = as.character(qmat_c2a$n1),
                            y_uc = as.character(qmat_c2a$x0),
                            n_uc = as.character(qmat_c2a$n0),
                            mean_iarr  = as.character(round(qmat_c2a$mean_iarr,3)*100),
                            arr_paste = paste0(round(qmat_c2a$arr,3)*100, " (",
                                               round(qmat_c2a$arr.hi,3)*100,",",
                                               round(qmat_c2a$arr.lo,3)*100,")"),
                            rr_paste = paste0(round(qmat_c2a$rr,2), " (",
                                               round(qmat_c2a$rr.lo,2),",",
                                               round(qmat_c2a$rr.hi,2),")"))

#dev.new(height=3,width=12,noRStudioGD = TRUE)

# pdf(file="C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/Lancet RM submission files/fig1.pdf",
#     height=3,
#     width=12,
#     onefile=FALSE)

tiff(file="C:/Users/VIT13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/CCM review/fig1_c2a_tiff.tif",
    height=3,
    width=12,
    units="in",
    res=1000)

base_data |>
  forestplot(labeltext = c(quintile, y_egdt, n_egdt, y_uc, n_uc, mean_iarr, arr_paste,rr_paste),
             clip = c(-20,20),
             xlog = FALSE,
             graphwidth="auto",
             colgap=unit(3,"mm"),
             txt_gp = fpTxtGp(ticks=gpar(cex=0.9),
                              xlab=gpar(cex=0.9)),
             lwd.xaxis=1.5,
             lwd.ci=2,
             lwd.zero=2,
             lty.zero=2,
             boxsize=0.3,
             cex=2,
             vertices=TRUE,
             graph.pos=7,
             xlab="% Absolute risk difference",
             xticks=c(-20,-15,-10,-5,0,5,10,15,20)) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_header(quintile = c("", "Quintile"),
                y_egdt = c("Deaths", "EGDT"),
                n_egdt = c("Total N", "EGDT"),
                y_uc = c("Deaths", "Usual Care"),
                n_uc = c("Total N", "Usual Care"),
                mean_iarr = c("Mean iARD","prediction"),
                arr_paste = c("% ARD", "(95% CI)"),
                rr_paste = c("Relative Risk", "(95% CI)")) |>
  # fp_append_row(mean  = 0.531,
  #               lower = 0.386,
  #               upper = 0.731,
  #               study = "Summary",
  #               OR = "0.53",
  #               is.summary = TRUE) |> 
  fp_set_zebra_style("#EFEFEF")

dev.off()




# Validation of Arise model


tiff(file="C:/Users/VIT13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/CCM review/fig1_a2c_tiff.tif",
     height=3,
     width=12,
     units="in",
     res=1000)

base_data <- tibble::tibble(mean  = round(qmat_a2c$arr,3)*100,
                            lower = round(qmat_a2c$arr.hi,3)*100,
                            upper = round(qmat_a2c$arr.lo,3)*100,
                            quintile = paste0("Q",1:5),
                            y_egdt = as.character(qmat_a2c$x1),
                            n_egdt = as.character(qmat_a2c$n1),
                            y_uc = as.character(qmat_a2c$x0),
                            n_uc = as.character(qmat_a2c$n0),
                            mean_iarr  = as.character(round(qmat_a2c$mean_iarr,3)*100),
                            arr_paste = paste0(round(qmat_a2c$arr,3)*100, " (",
                                               round(qmat_a2c$arr.hi,3)*100,",",
                                               round(qmat_a2c$arr.lo,3)*100,")"),
                            rr_paste = paste0(round(qmat_a2c$rr,2), " (",
                                              round(qmat_a2c$rr.lo,2),",",
                                              round(qmat_a2c$rr.hi,2),")"))

base_data |>
  forestplot(labeltext = c(quintile, y_egdt, n_egdt, y_uc, n_uc, mean_iarr, arr_paste,rr_paste),
             clip = c(-25,25),
             xlog = FALSE,
             graphwidth="auto",
             colgap=unit(3,"mm"),
             txt_gp = fpTxtGp(ticks=gpar(cex=0.9),
                              xlab=gpar(cex=0.9)),
             lwd.xaxis=1.5,
             lwd.ci=2,
             lwd.zero=2,
             lty.zero=2,
             boxsize=0.3,
             cex=2,
             vertices=TRUE,
             graph.pos=7,
             xlab="% Absolute risk difference",
             xticks=c(-25,-20,-15,-10,-5,0,5,10,15,20,25)) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |>
  fp_add_header(quintile = c("", "Quintile"),
                y_egdt = c("Deaths", "EGDT"),
                n_egdt = c("Total N", "EGDT"),
                y_uc = c("Deaths", "Usual Care"),
                n_uc = c("Total N", "Usual Care"),
                mean_iarr = c("Mean iARD","prediction"),
                arr_paste = c("% ARD", "(95% CI)"),
                rr_paste = c("Relative Risk", "(95% CI)")) |>
  # fp_append_row(mean  = 0.531,
  #               lower = 0.386,
  #               upper = 0.731,
  #               study = "Summary",
  #               OR = "0.53",
  #               is.summary = TRUE) |> 
  fp_set_zebra_style("#EFEFEF")

dev.off()

