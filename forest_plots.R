
library(forestplot)
library(dplyr)


#machine<-"D:/Box Sync/Box Sync/Box Sync/"
machine<-"C:/Users/vit13/OneDrive - University of Pittsburgh/"
qmat_e2f<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/quantile_mat_e2f_v2.RDS"))
qmat_f2e<-readRDS(paste0(machine,"_BoxMigration/Faraaz PROWESS/datasets/EGDT Trials/quantile_mat_f2e_v2.RDS"))


# Validation of trial E model (trial A in the paper; Figure 1)

base_data <- tibble::tibble(mean  = round(qmat_e2f$arr,3)*100,
                            lower = round(qmat_e2f$arr.lo,3)*100,
                            upper = round(qmat_e2f$arr.hi,3)*100,
                            quintile = paste0("Q",1:5),
                            y_egdt = as.character(qmat_e2f$x1),
                            n_egdt = as.character(qmat_e2f$n1),
                            y_uc = as.character(qmat_e2f$x0),
                            n_uc = as.character(qmat_e2f$n0),
                            mean_iarr  = as.character(round(qmat_e2f$mean_iarr,3)*100),
                            arr_paste = paste0(round(qmat_e2f$arr,3)*100, " (",
                                               round(qmat_e2f$arr.lo,3)*100,",",
                                               round(qmat_e2f$arr.hi,3)*100,")"),
                            rr_paste = paste0(round(qmat_e2f$rr,2), " (",
                                               round(qmat_e2f$rr.lo,2),",",
                                               round(qmat_e2f$rr.hi,2),")"))

#dev.new(height=3,width=12,noRStudioGD = TRUE)

# pdf(file="C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/Summer 2023/Lancet RM submission files/fig1.pdf",
#     height=3,
#     width=12,
#     onefile=FALSE)

tiff(file="C:/Users/vit13/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/CCM review/fig1_tiff.tif",
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




# Validation of trial F model (trial B in the paper; eFigure)

base_data <- tibble::tibble(mean  = round(qmat_f2e$arr,3)*100,
                            lower = round(qmat_f2e$arr.lo,3)*100,
                            upper = round(qmat_f2e$arr.hi,3)*100,
                            quintile = paste0("Q",1:5),
                            y_egdt = as.character(qmat_f2e$x1),
                            n_egdt = as.character(qmat_f2e$n1),
                            y_uc = as.character(qmat_f2e$x0),
                            n_uc = as.character(qmat_f2e$n0),
                            mean_iarr  = as.character(round(qmat_f2e$mean_iarr,3)*100),
                            arr_paste = paste0(round(qmat_f2e$arr,3)*100, " (",
                                               round(qmat_f2e$arr.lo,3)*100,",",
                                               round(qmat_f2e$arr.hi,3)*100,")"),
                            rr_paste = paste0(round(qmat_f2e$rr,2), " (",
                                              round(qmat_f2e$rr.lo,2),",",
                                              round(qmat_f2e$rr.hi,2),")"))

tiff(file="C:/Users/vbtal/OneDrive - University of Pittsburgh/CATES in Critical Care/EGDT manuscript 1/CCM review/fig1_f2e_tiff.tif",
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

