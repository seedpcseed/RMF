selbal.aux2 <- function (x, y, th.imp = 0, covar = NULL, logit.acc = "AUC", 
                         logt = T, maxV = 1e+10, zero.rep = "bayes") 
{
  classy <- class(y)
  f.class <- "gaussian"
  numy <- y
  if (classy == "factor") {
    ylev <- levels(y)
    numy <- as.numeric(y) - 1
    f.class <- "binomial"
  }
  suppressMessages(library(zCompositions))
  var.nam <- rem.nam <- colnames(x)
  if (!is.null(covar)) {
    dat <- data.frame(cbind(numy, covar))
  }else {
    dat <- data.frame(numy)
  }
  if (logt == F) {
    logCounts <- x
  }else {
    logCounts <- log(x + 1)
  }
  first.bal <- function(logCounts, Y, covar = NULL) {
    n <- ncol(logCounts)
    nam <- colnames(logCounts)
    if (classy == "factor") {
      M <- matrix(0, nrow = n, ncol = n)
    }
    else {
      M <- matrix(1e+99, nrow = n, ncol = n)
    }
    row.names(M) <- colnames(M) <- nam
    if (classy == "factor") {
      suppressWarnings(suppressMessages(library("CMA")))
      suppressMessages(detach("package:CMA", unload = TRUE))
      suppressMessages(library(pROC))
      for (i in 2:n) {
        for (j in 1:(i - 1)) {
          TAB <- data.frame(cbind(Y, logCounts[, i] - 
                                    logCounts[, j], covar))
          FIT <- glm(Y ~ ., data = TAB, family = f.class)
          ifelse(FIT$coefficients[2] > 0, M[i, j] <- logit.cor(FIT, 
                                                               y = y, covar = covar, logit.acc), M[j, i] <- logit.cor(FIT, 
                                                                                                                      y = y, covar = covar, logit.acc))
        }
      }
      r <- which(M == max(M), arr.ind = T)
    }
    else {
      for (i in 2:n) {
        for (j in 1:(i - 1)) {
          TAB <- data.frame(cbind(Y, logCounts[, i] - 
                                    logCounts[, j], covar))
          FIT <- glm(Y ~ ., data = TAB, family = f.class)
          ifelse(FIT$coefficients[2] > 0, M[i, j] <- mean(FIT$residuals^2), 
                 M[j, i] <- mean(FIT$residuals^2))
        }
      }
      r <- which(M == min(M), arr.ind = T)
    }
    return(r)
  }
  
  balance.adding.cor <- function(x, LogCounts, POS, NEG, numy, 
                                 covar = NULL) {
    S1.pos <- rowM(LogCounts[, c(POS, x)])
    s1 <- length(POS) + 1
    S2.pos <- rowM(LogCounts[, NEG])
    s2 <- length(NEG)
    BAL <- sqrt((s1 * s2)/(s1 + s2)) * (S1.pos - S2.pos)
    D.pos <- data.frame(cbind(numy, BAL, covar))
    FIT.pos <- glm(numy ~ ., data = D.pos, family = f.class)
    if (classy == "numeric") {
      C.pos <- mean(FIT.pos$residuals^2)
    }
    else {
      C.pos <- logit.cor(FIT.pos, y, covar = covar, logit.acc)
    }
    S1.neg <- rowM(LogCounts[, POS])
    s1 <- length(POS)
    S2.neg <- rowM(LogCounts[, c(NEG, x)])
    s2 <- length(NEG) + 1
    BAL <- sqrt((s1 * s2)/(s1 + s2)) * (S1.neg - S2.neg)
    D.neg <- data.frame(cbind(numy, BAL, covar))
    FIT.neg <- glm(numy ~ ., data = D.neg, family = f.class)
    if (classy == "numeric") {
      C.neg <- mean(FIT.neg$residuals^2)
    }
    else {
      C.neg <- logit.cor(FIT.neg, y, covar = covar, logit.acc)
    }
    COR <- c(C.pos, C.neg)
    return(COR)
  }
  A1 <- first.bal(logCounts, Y = numy, covar = covar)
  POS <- colnames(x)[A1[1, 1]]
  NEG <- colnames(x)[A1[1, 2]]
  INC.VAR <- c(POS, NEG)
  rem.nam <- setdiff(rem.nam, INC.VAR)
  S1 <- logCounts[, POS]
  S2 <- logCounts[, NEG]
  B <- sqrt(1/2) * (S1 - S2)
  Tab.var <- data.frame(Taxa = c(POS, NEG), Group = c("NUM", 
                                                      "DEN"))
  Tab.var[, 1] <- as.character(Tab.var[, 1])
  dat.ini <- cbind(dat, B)
  FIT.initial <- glm(numy ~ ., data = dat.ini, family = f.class)
  suppressWarnings(suppressMessages(library("CMA")))
  suppressMessages(detach("package:CMA", unload = TRUE))
  suppressMessages(library(pROC))
  if (classy == "numeric") {
    ACC.Bal <- mean(FIT.initial$residuals^2)
  }
  else {
    ACC.Bal <- logit.cor(FIT.initial, y, covar = covar, logit.acc)
  }
  ACC.ref <- ACC.Bal
  ACC.set <- ACC.ref
  nV <- 2
  if (classy == "numeric") {
    while (ACC.set <= ACC.ref && length(rem.nam) != 0 && 
           nV < maxV) {
      ACC.ref <- ACC.set
      add2bal.ACC <- matrix(, nrow = 0, ncol = 2)
      suppressWarnings(suppressMessages(library("CMA")))
      suppressMessages(detach("package:CMA", unload = TRUE))
      suppressMessages(library(pROC))
      add2bal.ACC <- t(sapply(rem.nam, function(x) balance.adding.cor(x, 
                                                                      LogCounts = logCounts, POS, NEG, numy = numy, 
                                                                      covar = covar)))
      row.names(add2bal.ACC) <- rem.nam
      ACC.opt <- which(add2bal.ACC == min(add2bal.ACC), 
                       arr.ind = T)[1, ]
      ACC.set <- min(add2bal.ACC)
      if ((ACC.set - ACC.ref) < th.imp) {
        INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
        nV <- nV + 1
        if (ACC.opt[2] == 1) {
          POS <- c(POS, rem.nam[ACC.opt[1]])
          Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], 
                                      "NUM"))
        }
        else if (ACC.opt[2] == 2) {
          NEG <- c(NEG, rem.nam[ACC.opt[1]])
          Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], 
                                      "DEN"))
        }
        else {
          ACC.set <- 0
        }
        rem.nam <- rem.nam[-ACC.opt[1]]
      }
    }
  }
  else {
    while (ACC.set >= ACC.ref && length(rem.nam) != 0 && 
           nV < maxV) {
      ACC.ref <- ACC.set
      add2bal.ACC <- matrix(, nrow = 0, ncol = 2)
      suppressWarnings(suppressMessages(library("CMA")))
      suppressMessages(detach("package:CMA", unload = TRUE))
      suppressMessages(library(pROC))
      add2bal.ACC <- t(sapply(rem.nam, function(x) balance.adding.cor(x, 
                                                                      LogCounts = logCounts, POS, NEG, numy = numy, 
                                                                      covar = covar)))
      row.names(add2bal.ACC) <- rem.nam
      ACC.opt <- which(add2bal.ACC == max(add2bal.ACC), 
                       arr.ind = T)[1, ]
      ACC.set <- max(add2bal.ACC)
      if ((ACC.set - ACC.ref) > th.imp) {
        INC.VAR <- c(INC.VAR, rem.nam[ACC.opt[1]])
        nV <- nV + 1
        if (ACC.opt[2] == 1) {
          POS <- c(POS, rem.nam[ACC.opt[1]])
          Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], 
                                      "NUM"))
        }
        else if (ACC.opt[2] == 2) {
          NEG <- c(NEG, rem.nam[ACC.opt[1]])
          Tab.var <- rbind(Tab.var, c(rem.nam[ACC.opt[1]], 
                                      "DEN"))
        }
        else {
          ACC.set <- 0
        }
      }
      rem.nam <- rem.nam[-ACC.opt[1]]
    }
  }
  k1 <- length(POS)
  k2 <- length(NEG)
  FINAL.BAL <- sqrt((k1 * k2)/(k1 + k2)) * (rowM(logCounts[, 
                                                           POS]) - rowM(logCounts[, NEG]))
  return(Tab.var)
}

selbal.cv2<- function (x, y, n.fold = 5, n.iter = 10, seed = 31415, covar = NULL, 
          col = c("steelblue1", "tomato1"), col2 = c("darkgreen", "steelblue4", 
                                                     "tan1"), logit.acc = "AUC", maxV = 20, zero.rep = "bayes", 
          opt.cri = "1se", user_numVar = NULL) 
{
  suppressMessages(library(plyr))
  classy <- class(y)
  f.class <- "gaussian"
  numy <- y
  if (classy == "factor") {
    ylev <- levels(y)
    numy <- as.numeric(y) - 1
    f.class <- "binomial"
  }
  x.nam <- colnames(x)
  ROB.TAB <- matrix(0, nrow = ncol(x), ncol = 3)
  colnames(ROB.TAB) <- c("Prop. Included", "Prop_Numerator", 
                         "Prop_Denominator")
  row.names(ROB.TAB) <- x.nam
  BAL.resume <- matrix(0, nrow = n.fold * n.iter, ncol = ncol(x))
  colnames(BAL.resume) <- colnames(x)
  if (!is.null(covar)) {
    dat <- data.frame(cbind(numy, covar))
  }else {
    dat <- data.frame(numy)
  }
  cat(paste("\n\n###############################################################", 
            "\n STARTING selbal.cv FUNCTION", "\n###############################################################"))
  cat(paste("\n\n#-------------------------------------------------------------#", 
            "\n# ZERO REPLACEMENT . . .\n\n"))
  logc <- log(x + 1)
  cat(paste("\n, . . . FINISHED.", "\n#-------------------------------------------------------------#"))
  suppressMessages(library(CMA))
  set.seed(seed)
  CV.groups <- GenerateLearningsets(y = y, fold = n.fold, niter = n.iter, 
                                    method = "CV", strat = ifelse(class(y) != "factor", F, 
                                                                  T))@learnmatrix
  suppressMessages(detach("package:CMA", unload = TRUE))
  suppressMessages(library(pROC))
  cv.MSE <- function(k) {
    CV.group <- CV.groups[((k - 1) * n.fold + 1):(k * n.fold), 
    ]
    Bal.List <- list()
    ACC.Mat <- matrix(0, nrow = maxV - 1, ncol = nrow(CV.group))
    for (i in 1:nrow(CV.group)) {
      train.idx <- CV.group[i, ]
      x.train <- logc[train.idx, ]
      x.test <- logc[-train.idx, ]
      y.train <- y[train.idx]
      y.test <- y[-train.idx]
      covar.train <- covar[train.idx, ]
      covar.test <- covar[-train.idx, ]
      BAL <- selbal.aux(x.train, y.train, th.imp = 0, covar = covar.train, 
                        logit.acc, logt = F, maxV = maxV)
      Bal.List[[i]] <- BAL
      PRED.y <- matrix(0, nrow = length(y.test), ncol = nrow(BAL) - 
                         1)
      for (l in 2:min(maxV, nrow(BAL))) {
        df <- data.frame(Y = y.train, B = bal.value(BAL[1:l, 
        ], x.train))
        df.test <- data.frame(Y = y.test, B = bal.value(BAL[1:l, 
        ], x.test))
        if (!is.null(covar)) {
          df <- cbind(df, cov = covar.train)
          df.test <- cbind(df.test, cov = covar.test)
        }
        FIT <- glm(Y ~ ., data = df, family = f.class)
        PRED.y[, l - 1] <- predict(FIT, df.test, type = "response")
      }
      ACC.eval <- function(y, pred, classy, logit.acc = NULL) {
        if (classy == "numeric") {
          ACC <- apply(pred, 2, function(x) mean((y - 
                                                    x)^2))
        }
        else {
          if (logit.acc == "AUC") {
            library(pROC)
            ACC <- apply(pred, 2, function(x) auc(y, 
                                                  x, quiet = TRUE))
          }
          else if (logit.acc == "Rsq") {
            ACC <- apply(pred, 2, function(x) cor(y, 
                                                  x)^2)
          }
          else if (logit.acc == "Tjur") {
            ACC <- apply(pred, 2, function(x) mean(x[y == 
                                                       1]) - mean(x[y == 0]))
          }
          else if (logit.acc == "Dev") {
            ACC <- apply(pred, 2, function(x) 1 - (deviance(glm(y ~ 
                                                                  x, data = df, family = binomial()))/glm(y ~ 
                                                                                                            1, family = binomial())[[10]]))
          }
        }
        return(ACC)
      }
      R <- ACC.eval(numy[-train.idx], PRED.y, classy = classy, 
                    logit.acc)
      ACC.Mat[, i] <- c(R, rep(R[length(R)], maxV - length(R) - 
                                 1))
    }
    return(list(Bal.List, ACC.Mat))
  }
  cat(paste("\n\n#-------------------------------------------------------------#", 
            "\n# Starting the cross - validation procedure . . ."))
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  no_cores <- detectCores() - 2
  registerDoParallel(no_cores)
  comb <- function(x, ...) {
    lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), 
                                                      function(y) y[[i]])))
  }
  INTEREST <- foreach(h = 1:n.iter, .export = c("logit.cor", 
                                                "rowM", "selbal.aux", "bal.value", "logit.acc", "cmultRepl", 
                                                "cmultRepl2"), .combine = "comb", .multicombine = TRUE, 
                      .init = list(list(), list())) %dopar% {
                        cv.MSE(h)
                      }
  stopImplicitCluster()
  cat(paste("\n . . . finished.", "\n#-------------------------------------------------------------#", 
            "\n###############################################################"))
  Balances <- unlist(INTEREST[[1]], recursive = F)
  ACC.Matrix <- do.call(cbind, INTEREST[[2]])
  ACC.mean <- apply(ACC.Matrix, 1, mean)
  ACC.se <- apply(ACC.Matrix, 1, function(x) sd(x)/sqrt(length(x)))
  if (classy == "numeric") {
    m <- which.min(ACC.mean)
    if (length(which((ACC.mean < (ACC.mean[m] + ACC.se[m])) == 
                     T)) > 0) {
      mv <- min(which((ACC.mean < (ACC.mean[m] + ACC.se[m])) == 
                        T))
    }
    else {
      mv <- m
    }
    if (opt.cri == "1se") {
      opt.M <- mv + 1
    }else {
      opt.M <- m + 1
    }
  }else {
    m <- which.max(ACC.mean)
    if (length(which((ACC.mean > (ACC.mean[m] - ACC.se[m])) == 
                     T)) > 0) {
      mv <- min(which((ACC.mean > (ACC.mean[m] - ACC.se[m])) == 
                        T))
    }
    else {
      mv <- m
    }
    if (opt.cri == "1se") {
      opt.M <- mv + 1
    }
    else {
      opt.M <- m + 1
    }
  }
  cat(paste("\n\n The optimal number of variables is:", opt.M, 
            "\n\n"))
  if (!is.null(user_numVar)) {
    opt.M <- user_numVar
  }
  suppressMessages(BAL <- selbal.aux2(x, y, th.imp = 0, covar = covar, 
                                     logit.acc, logt = T, maxV = opt.M))
  NUM <- BAL[BAL[, 2] == "NUM", "Taxa"]
  DEN <- BAL[BAL[, 2] == "DEN", "Taxa"]
  k1 <- length(NUM)
  k2 <- length(DEN)
  FINAL.BAL <- sqrt((k1 * k2)/(k1 + k2)) * (rowM(logc[, NUM]) - 
                                              rowM(logc[, DEN]))
  U <- data.frame(dat, FINAL.BAL)
  colnames(U)[ncol(U)] <- "V1"
  FIT.final <- glm(numy ~ ., data = U, family = f.class)
  if (classy == "numeric") {
    df.boxplot <- data.frame(mean = ACC.mean, se = ACC.se, 
                             n = 2:maxV)
    library(ggplot2)
    MSE.Boxplot <- ggplot(df.boxplot, aes(x = n, y = mean)) + 
      geom_errorbar(aes(ymin = mean - se, ymax = mean + 
                          se), width = 0.2, col = "gray") + geom_vline(xintercept = opt.M, 
                                                                       linetype = "dotted", col = "blue") + geom_point(color = "red", 
                                                                                                                       lwd = 2) + theme_bw() + xlab("Number of variables") + 
      ylab("Mean-Squared Error") + scale_x_continuous(breaks = seq(2, 
                                                                   maxV, 1)) + theme(strip.text.x = element_text(size = 12, 
                                                                                                                 angle = 0, face = "bold", colour = "white"), strip.text.y = element_text(size = 12, 
                                                                                                                                                                                          face = "bold"), strip.background = element_rect(colour = "black", 
                                                                                                                                                                                                                                          fill = "black"), plot.title = element_text(size = 20, 
                                                                                                                                                                                                                                                                                     vjust = 2.25, hjust = 0.5, face = "bold"), panel.grid.major = element_blank(), 
                                                                                     panel.grid.minor = element_blank())
  }
  else {
    df.boxplot <- data.frame(mean = ACC.mean, se = ACC.se, 
                             n = 2:maxV)
    ylabelName = "Accuracy (AUC)"
    if (logit.acc == "Dev") {
      ylabelName = "Explained Deviance"
    }
    library(ggplot2)
    MSE.Boxplot <- ggplot(df.boxplot, aes(x = n, y = mean)) + 
      geom_errorbar(aes(ymin = mean - se, ymax = mean + 
                          se), width = 0.2, col = "gray") + geom_vline(xintercept = opt.M, 
                                                                       linetype = "dotted", col = "blue") + geom_point(color = "red", 
                                                                                                                       lwd = 2) + theme_bw() + xlab("Number of variables") + 
      ylab(ylabelName) + scale_x_continuous(breaks = seq(2, 
                                                         maxV, 1)) + theme(strip.text.x = element_text(size = 12, 
                                                                                                       angle = 0, face = "bold", colour = "white"), strip.text.y = element_text(size = 12, 
                                                                                                                                                                                face = "bold"), strip.background = element_rect(colour = "black", 
                                                                                                                                                                                                                                fill = "black"), plot.title = element_text(size = 20, 
                                                                                                                                                                                                                                                                           vjust = 2.25, hjust = 0.5, face = "bold"), panel.grid.major = element_blank(), 
                                                                           panel.grid.minor = element_blank())
  }
  Sub.Balances <- lapply(Balances, function(x) x[1:min(nrow(x), 
                                                       opt.M), ])
  BR <- matrix(0, nrow = n.fold * n.iter, ncol = length(x.nam))
  colnames(BR) <- x.nam
  for (i in 1:length(Sub.Balances)) {
    BR[i, colnames(BR) %in% Sub.Balances[[i]][Sub.Balances[[i]][, 
                                                                "Group"] == "NUM", "Taxa"]] <- 1
    BR[i, colnames(BR) %in% Sub.Balances[[i]][Sub.Balances[[i]][, 
                                                                "Group"] == "DEN", "Taxa"]] <- 2
  }
  ROB.TAB[, 1] <- apply(BR != 0, 2, function(x) 100 * mean(x))
  ROB.TAB[, 2] <- apply(BR == 1, 2, function(x) 100 * mean(x))
  ROB.TAB[, 3] <- apply(BR == 2, 2, function(x) 100 * mean(x))
  fil1 <- which(ROB.TAB[, 1] != 0)
  ord.ROB.TAB <- ROB.TAB[fil1, ]
  sel.ord <- order(ord.ROB.TAB[, 1], decreasing = F)
  ord.ROB.TAB <- ord.ROB.TAB[sel.ord, ]
  BAL.SEL.TAB <- data.frame(name = row.names(ord.ROB.TAB), 
                            sel = ord.ROB.TAB[, 1])
  BAL.SEL.TAB$name <- factor(BAL.SEL.TAB$name, levels = BAL.SEL.TAB$name)
  COLOR.BAL <- rep(col[2], nrow(BAL.SEL.TAB))
  vDEN <- row.names(ord.ROB.TAB)[which(ord.ROB.TAB[, "Prop_Denominator"] != 
                                         0)]
  COLOR.BAL[row.names(BAL.SEL.TAB) %in% vDEN] <- col[1]
  BAL.SEL.TAB$COLOR.BAL <- factor(COLOR.BAL, levels = col, 
                                  labels = col)
  suppressMessages(library(ggplot2))
  IMP.plot <- ggplot(BAL.SEL.TAB, aes(x = factor(name), y = sel)) + 
    geom_bar(stat = "identity", aes(fill = COLOR.BAL), size = 1) + 
    guides(size = FALSE) + scale_fill_manual(name = "Group of . . .", 
                                             values = c(col[1], col[2]), breaks = c(col[1], col[2]), 
                                             labels = c("DEN", "NUM")) + scale_color_manual(name = "Variables \n appearing in . . .", 
                                                                                            values = c(col2, "white"), breaks = c(col2, "white"), 
                                                                                            labels = c("Both", "Global", "CV", "NONE"), drop = F, 
                                                                                            guide = guide_legend(override.aes = list(fill = "gray90"))) + 
    ylab("% of times included in a balance") + xlab("") + 
    theme_bw() + coord_flip() + ggtitle("Cross validation in balance selection") + 
    theme(strip.text.x = element_text(size = 12, angle = 0, 
                                      face = "bold", colour = "white"), strip.text.y = element_text(size = 12, 
                                                                                                    face = "bold"), strip.background = element_rect(colour = "black", 
                                                                                                                                                    fill = "black"), plot.title = element_text(size = 20, 
                                                                                                                                                                                               vjust = 2.25, hjust = 0.5, face = "bold"), legend.title = element_text(face = "bold"), 
          legend.text = element_text(face = "bold"))
  BAL.str <- apply(BR, 1, function(x) paste(x, collapse = ""))
  BAL.tab <- prop.table(table(BAL.str))
  nam.str <- names(BAL.tab)
  nam.A <- t(sapply(nam.str, FUN = function(x) unlist(strsplit(x, 
                                                               ""))))
  INC <- apply(nam.A, 1, function(x) x.nam[x != 0])
  INC.NUM <- alply(nam.A, 1, function(x) x.nam[x == 1])
  INC.DEN <- alply(nam.A, 1, function(x) x.nam[x == 2])
  UNIQUE.VAR <- unique(c(as.vector(unlist(INC)), NUM, DEN))
  RESUME.BAL <- as.data.frame(matrix(0, nrow = length(UNIQUE.VAR), 
                                     ncol = length(BAL.tab)))
  row.names(RESUME.BAL) <- UNIQUE.VAR
  RESUME.BAL[sapply(INC.NUM, function(x) UNIQUE.VAR %in% x)] <- "NUM"
  RESUME.BAL[sapply(INC.DEN, function(x) UNIQUE.VAR %in% x)] <- "DEN"
  RESUME.BAL <- rbind(RESUME.BAL, FREQ = as.numeric(BAL.tab))
  RESUME.BAL <- cbind(RESUME.BAL, 0, 0)
  RESUME.BAL[row.names(RESUME.BAL) %in% NUM, ncol(RESUME.BAL)] <- "NUM"
  RESUME.BAL[row.names(RESUME.BAL) %in% DEN, ncol(RESUME.BAL)] <- "DEN"
  RESUME.BAL[-nrow(RESUME.BAL), ncol(RESUME.BAL) - 1] <- ROB.TAB[row.names(RESUME.BAL)[-nrow(RESUME.BAL)], 
                                                                 1]
  RESUME.BAL <- RESUME.BAL[, c(ncol(RESUME.BAL), ncol(RESUME.BAL) - 
                                 1, order(RESUME.BAL[nrow(RESUME.BAL), -c(ncol(RESUME.BAL), 
                                                                          ncol(RESUME.BAL) - 1)], decreasing = T))]
  RESUME.BAL[nrow(RESUME.BAL), 1:2] <- "-"
  dat <- RESUME.BAL[, c(1, 2:(min(5, ncol(RESUME.BAL))))]
  W <- which(apply(dat[, -2] == 0, 1, mean) == 1)
  if (length(W) != 0) {
    dat <- dat[-as.numeric(W), ]
  }
  dat <- dat[, c(2, 1, 3:ncol(dat))]
  colnames(dat)[1:2] <- c("%", "Global")
  dat <- dat[c(order(as.numeric(dat[-nrow(dat), 1]), decreasing = T), 
               nrow(dat)), ]
  ifelse(classy %in% c("numeric", "integer"), y.plot <- y, 
         y.plot <- 1:length(y))
  PLOT.Global <- plot.bal(NUM, DEN, logc, y, covar, col = col, 
                          logit.acc)
  cat(paste("\n\n###############################################################", 
            "\n . . . FINISHED.", "\n###############################################################"))
  L <- list(accuracy.nvar = MSE.Boxplot, var.barplot = IMP.plot, 
            global.plot = PLOT.Global$Global.plot, global.plot2 = PLOT.Global$Global.plot2, 
            ROC.plot = PLOT.Global$ROC.plot, cv.tab = dat, cv.accuracy = ACC.Matrix[(opt.M - 
                                                                                       1), ], global.balance = BAL, glm = FIT.final, opt.nvar = opt.M)
  return(L)
}


