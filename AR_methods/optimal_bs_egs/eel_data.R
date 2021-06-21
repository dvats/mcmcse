load("optimal_bs_egs/truth_eel_1e3")
library(dismo)
library(MCMCpack)

data(Anguilla_train)

ang_sub <- data.frame(Anguilla_train$Angaus, Anguilla_train$SegSumT, Anguilla_train$DSDist, Anguilla_train$USNative,
	Anguilla_train$Method, Anguilla_train$DSMaxSlope, Anguilla_train$USSlope)

names(ang_sub) <- c("Angaus", "SegSumT", "DSDist", "USNative", "Method", "DSMaxSlope", "USSlope")

dat <- ang_sub
truth <- colMeans(mean.reps)