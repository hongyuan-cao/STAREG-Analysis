library(SPARK)
library(STAREG)
library(ggplot2)

source('./funcs/SimuFunc.R')
source("./funcs/ROC_funcs.R")

tau1     = 0.5
tau2     = 0.8
ipt      = "Pattern I"
n.gene   = 10000
sig1.str = "Strong"
sig2.str = "Moderate"
xi00 = 0.9
xi01 = 0.025
xi10 = 0.025
xi11 = 0.05

h = sample(0:3, n.gene, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
states1 = rep(0, n.gene)
states1[which(h==2|h==3)] = 1
states2 = rep(0, n.gene)
states2[which(h==1|h==3)] = 1
data.obj <- SimuData(tau1 = tau1,
                     tau2 = tau2,
                     ipt = ipt,
                     sig1.str = sig1.str,
                     sig2.str = sig2.str,
                     truth1   = states1,
                     truth2   = states2)
closeAllConnections()

p1 <- data.obj$pvals1
p2 <- data.obj$pvals2
states1 <- data.obj$truth1
states2 <- data.obj$truth2

# BH
padj1.bh <- p.adjust(p1, method = 'BH')
padj2.bh <- p.adjust(p2, method = 'BH')

# MaxP
maxp <- apply(cbind(p1, p2), 1, max)
padj.maxp <- p.adjust(maxp, method = "BH")

# STAREG
res.rep <- STAREG(p1, p2, trace = FALSE)
padj.rep <- res.rep$fdr.rep
f1.rep = res.rep$f1
f2.rep = res.rep$f2
xi00.hat = res.rep$xi00
xi01.hat = res.rep$xi01
xi10.hat = res.rep$xi10
xi11.hat = res.rep$xi11

##--------------------------------------------------------------------
rroc.bh <- revisedROC2(states1*states2, padj1.bh, padj2.bh)
rroc.maxp <- revisedROC(states1*states2, padj.maxp)
rroc.rep <- revisedROC(states1*states2, padj.rep)

rroc.bh <- as.data.frame(rroc.bh)
rroc.maxp <- as.data.frame(rroc.maxp)
rroc.rep <- as.data.frame(rroc.rep)

rroc.bh$method <- "BH"
rroc.maxp$method <- "MaxP"
rroc.rep$method <- "STAREG"

res <- rbind(rroc.maxp, rroc.bh, rroc.rep)

ggplot(res, aes(fdr,tpr, group = method, color = method, shape = method)) +
  scale_colour_manual(name="",
                      values = c("MaxP"="#6496D2", "BH"="#F4B183", "STAREG"="#8FBC8F")) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_line(size = 1.5) + xlab("Empirical FDR") + ylab("Power") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, colour = "#44546A", size = 1, linetype = 2) +
  theme(legend.position = "none", # legend.position = c(0.01,1.03)
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.6, linetype = "solid"))
