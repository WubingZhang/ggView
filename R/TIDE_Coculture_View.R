#' Visualize the selection level of gene in TIDE coculture screens.
#'
#' @docType methods
#' @name TIDE_Coculture_View
#' @rdname TIDE_Coculture_View
#' @param gene A vector of gene names.
#' @param PanDir The path to Pan's screen.
#' @param MangusoDir1 The path to Manguso's screen.
#' @examples
#' cocultureView("Traf3", width = 6, height = 3)
#' @export
TIDE_Coculture_View <- function(gene = "Traf3",
                                PanDir = "~/Jobs/Archive/CRISPRscreens/TIDE_coculture/Pan_2018_normcount.txt",
                                MangusoDir1 = "~/Jobs/Archive/CRISPRscreens/TIDE_coculture/Manguso_2017/Manguso_lfc.txt",
                                IgGDir = "~/Jobs/Archive/CRISPRscreens/TIDE_coculture/Kearney_2018/test/MC38_IgGvsT30.sgrna_summary.txt",
                                PD1Dir = "~/Jobs/Archive/CRISPRscreens/TIDE_coculture/Kearney_2018/test/MC38_PD1vsT30.sgrna_summary.txt",
                                NK10Dir = "~/Jobs/Archive/CRISPRscreens/TIDE_coculture/Kearney_2018/test/NK_Hit10vsTend.sgrna_summary.txt",
                                NK20Dir = "~/Jobs/Archive/CRISPRscreens/TIDE_coculture/Kearney_2018/test/NK_Hit20vsTend.sgrna_summary.txt"){
  # (1) B16F10 cells were harvested for genomic DNA isolation prior to the screen;
  # (2) B16F10 cells + control OT-I T cells at a 1:1 ratio (control condition, no Ova peptide);
  # (3) B16F10 cells were co-cultured with Pmel-1 T cells at a 1:1 ratio (experimental condition, Ova peptide).
  Pan = read.table(PanDir, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  genes = unique(gsub("@.*", "", rownames(Pan)))
  Pan$OT1.lfc = Pan$`OT-1_1st;OVA-+`-Pan$`OT-1_1st;OVA-`
  Pan$Pmel.lfc = Pan$`Pmel-1_1st;IFNg+`-Pan$`Pmel-1_1st;OT1_Ctrl_IFNg+`

  ## Pmel1: sgRNA rank figures
  gg1 = data.frame(sgrna = rownames(Pan), Gene = gsub("@.*", "", rownames(Pan)),
                  LFC = Pan$Pmel.lfc)
  gg1 = gg1[gg1$Gene%in%gene, ]
  gg1$Gene = "Pan2018_Pmel1"

  ## OT1: sgRNA rank figures
  gg2 = data.frame(sgrna = rownames(Pan), Gene = gsub("@.*", "", rownames(Pan)),
                  LFC = Pan$OT1.lfc)
  gg2 = gg2[gg2$Gene%in%gene, ]; gg2$Gene = "Pan2018_OT1"
  # gg$LFC[gg$LFC>xlim[2]] = xlim[2]
  # gg$LFC[gg$LFC< xlim[1]] = xlim[1]
  # p = sgRankView(gg, gene = gene, top = 0, bottom = 0)
  # p = p + labs(title = "B16F10-OT1 system")
  # p = p + xlim(xlim[1], xlim[2])
  # ggsave(paste0(outdir, "sgRankView_OT1_", gene[1], "_", date, ".png"), p,
  #        width = 4, height = 3, dpi = 200)

  # ## Pmel1: Gene rank figures
  # library(data.table)
  # tmp = Pan
  # tmp$Gene = gsub("@.*", "", rownames(tmp))
  # tmp = setDT(tmp)
  # Pmel = tmp[, .(LFC = median(Pmel.lfc)), by = Gene]
  # rankdata = Pmel$LFC; names(rankdata) = Pmel$Gene
  # p = RankView(rankdata, genelist = gene, top = 6, bottom = 0)
  # p = p + labs(x = "Log2(Fold enrichment)", title = "B16F10-Pmel1 system")
  # ggsave(paste0(outdir, "RankView_Pmel1_", gene[1], "_", date, ".png"), p,
  #        width = 4, height = 3, dpi = 200)
  #
  # ## OT1: Gene rank figures
  # OT1 = tmp[, .(LFC = median(OT1.lfc)), by = Gene]
  # rankdata = OT1$LFC; names(rankdata) = OT1$Gene
  # p = RankView(rankdata, genelist = gene, top = 6, bottom = 0)
  # p = p + labs(x = "Log2(Fold enrichment)", title = "B16F10-Pmel1 system")
  # ggsave(paste0(outdir, "RankView_OT1_", gene[1], "_", date, ".png"), p,
  #        width = 4, height = 3, dpi = 200)

  # The library was divided into four sub-pools, each containing one sgRNA per gene and 100 non-targeting control sgRNAs.
  # For each sub-pool, B16 cells were implanted into 10 Tcra−/− mice, 10 wild-type mice treated with GVAX,
  # and 10 wild-type mice treated with GVAX and PD-1 blockade.
  Manguso_1 = read.table(MangusoDir1, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  gg3 = Manguso_1[Manguso_1$Symbol%in%gene, ]
  gg3 = gg3[,-2]
  colnames(gg3) = c("sgrna", "Manguso2017_GVAX", "Manguso2017_GVAX+PD-1")
  gg3 = data.table::melt(gg3, id.vars = "sgrna", variable.name="Gene", value.name = "LFC")


  ## Kearney 2018
  MC38_IgG = ReadsgRRA(IgGDir)
  MC38_PD1 = ReadsgRRA(PD1Dir)
  NK_10 = ReadsgRRA(NK10Dir)
  NK_20 = ReadsgRRA(NK20Dir)
  MC38_IgG = MC38_IgG[MC38_IgG$Gene%in%gene, 1:3]
  MC38_IgG$Gene = "Kearney2018_T_IgG"
  MC38_PD1 = MC38_PD1[MC38_PD1$Gene%in%gene, 1:3]
  MC38_PD1$Gene = "Kearney2018_T_PD1"
  NK_10 = NK_10[NK_10$Gene%in%gene, 1:3]
  NK_10$Gene = "Kearney2018_NK_10"
  NK_20 = NK_20[NK_20$Gene%in%gene, 1:3]
  NK_20$Gene = "Kearney2018_NK_20"
  gg = rbind.data.frame(NK_10, NK_20, gg1, gg2, gg3, MC38_IgG, MC38_PD1)
  p = sgRankView(gg, gene = intersect(c("Pan2018_OT1", "Pan2018_Pmel1", "Manguso2017_GVAX",
                                        "Manguso2017_GVAX+PD-1", "Kearney2018_T_IgG", "Kearney2018_T_PD1",
                                        "Kearney2018_NK_10", "Kearney2018_NK_20"), gg$Gene),
                 top = 0, bottom = 0)
  p = p + labs(title = paste0("Ref_Screen - ", gene[1]))
  p = p + xlim(-max(abs(gg$LFC)), max(abs(gg$LFC)))
  p
}
