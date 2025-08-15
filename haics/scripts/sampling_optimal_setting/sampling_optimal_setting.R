# Script that calculates and save the scaling values and HSLs

#function from analysis
fnReadFile = function(strFilePath)
{
  ##check arguments---------------------------------------#
  if(length(strFilePath)!=1 || !is.character(strFilePath) || !file.access(strFilePath, mode = 4)==0) stop("\n ERROR: couldn't find ",strFilePath)
  
  ##main--------------------------------------------------#
  dfInput = read.table(file=strFilePath, header = F, sep = "", colClasses=c("character", "character"), comment.char="#")
  lsInput = setNames(split(dfInput[, 2], seq(nrow(dfInput))), nm = dfInput[, 1])
  
  ##output------------------------------------------------#
  return(lsInput)
}

cat("\n\tReading arguments...")

#internally or externally provided?
args = commandArgs(TRUE)

ID_BATCH = as.character(args[1])
NUM_OF_RUNS = as.numeric(args[2])
ITER_SAMPLING = as.numeric(args[3])
ID_NEW = as.character(args[4])

cat("\n\tChecking arguments...")

dirData = paste("data/", ID_BATCH, sep = "")
if( !dir.exists(dirData) ) stop("\n ERROR: couldn't find \t", dirData)

fileHsl = paste(dirData, "/hsl-", ID_BATCH, ".txt", sep = "")
if( !file.exists(fileHsl) ) stop("\n ERROR: couldn't find", fileHsl)

fileHslFiltAR = paste(dirData, "/hsl_filtAR-", ID_BATCH, ".txt", sep = "")
if( !file.exists(fileHslFiltAR) ) stop("\n ERROR: couldn't find", fileHslFiltAR)

fileHslGupta = paste(dirData, "/hsl_Gupta-", ID_BATCH, ".txt", sep = "")
if( !file.exists(fileHslGupta) ) stop("\n ERROR: couldn't find", fileHslGupta)

fileFittingFactor = paste(dirData, "/fitting_factor-", ID_BATCH, ".txt", sep = "")
if( !file.exists(fileFittingFactor) ) stop("\n ERROR: couldn't find", fileFittingFactor)

fileInputfile = paste(dirData, "/inputfile_", ID_BATCH, ".txt", sep = "")
if( !file.exists(fileInputfile) ) stop("\n ERROR: couldn't find", fileInputfile)

#read the hsl files
dfHsl = read.table(fileHsl)
dfHsl_filtAR = read.table(fileHslFiltAR)
dfHsl_Gupta = read.table(fileHslGupta)

#read the fitting factor file
dfFittingFactor = read.table(fileFittingFactor)

#read the inputfile, most of the entries will be the same
lsInput = fnReadFile(strFilePath = fileInputfile)

dfInput = read.table(fileInputfile)

#directory of the new simulation
dirSimulation = paste("../../output/", ID_NEW, sep = "")

#creating directory
dirInput = paste(dirSimulation, "/input/", sep = "")
if(!dir.exists(dirSimulation)){
  invisible(dir.create(dirSimulation))
  invisible(dir.create(dirInput))
}

#common parameters setting
model = lsInput$model
data = lsInput$data
method = "GHMC"
seed = lsInput$seed
iter_sampling = ITER_SAMPLING
iter_burn_in = "0"
integrator = "3sAIA-HMC"
t_L = "1"
t_stepsize = "3"
thinning = lsInput$thinning
t_Phi = "1"
Phi = "0.5"
iter_tune = lsInput$iter_tune
AR_target = lsInput$AR_target
delta_AR_target = lsInput$delta_AR_target
num_basis_spline = as.numeric(lsInput$num_basis_spline)
spline_pol_degree = lsInput$spline_pol_degree
tfin_spl_bas = lsInput$tfin_spl_bas
num_out_cvodes = lsInput$num_out_cvodes
S0_prior = lsInput$S0_prior
S0_param_1 = lsInput$S0_param_1
S0_param_2 = lsInput$S0_param_2
E0_prior = lsInput$E0_prior
E0_param_1 = lsInput$E0_param_1
E0_param_2 = lsInput$E0_param_2
I0_prior = lsInput$I0_prior
I0_param_1 = lsInput$I0_param_1
I0_param_2 = lsInput$I0_param_2
alpha_prior = lsInput$alpha_prior
alpha_param_1 = lsInput$alpha_param_1
alpha_param_2 = lsInput$alpha_param_2
gamma_prior = lsInput$gamma_prior
gamma_param_1 = lsInput$gamma_param_1
gamma_param_2 = lsInput$gamma_param_2
phi_inv_prior = lsInput$phi_inv_prior
phi_inv_param_1 = lsInput$phi_inv_param_1
phi_inv_param_2 = lsInput$phi_inv_param_2
tau_prior = lsInput$tau_prior
tau_param_1 = lsInput$tau_param_1
tau_param_2 = lsInput$tau_param_2
num_comp = lsInput$num_comp
num_comp_E = lsInput$num_comp_E
do_under_report = lsInput$do_under_report
T0_under = lsInput$T0_under
U0_under = lsInput$U0_under
T1_under = lsInput$T1_under
U1_under = lsInput$U1_under
do_grad_desc = lsInput$do_grad_desc
learning_rate = lsInput$learning_rate
max_iter = lsInput$max_iter
gammaFixed = lsInput$gammaFixed
gamma_fixed_value = lsInput$gamma_fixed_value
gammaBounded = lsInput$gammaBounded
gammaUpper = lsInput$gammaUpper
gammaLower = lsInput$gammaLower
beta_prior = lsInput$beta_prior
beta_param_1 = lsInput$beta_param_1
beta_param_2 = lsInput$beta_param_2

for (i in c(1:NUM_OF_RUNS)) {
  
  #Write the initialpoint and splinebasis files
  fileFinalPoint = paste(dirData, "/finalpoint_", i, ".txt", sep = "")
  
  new_points = suppressWarnings(read.table(fileFinalPoint))
  
  D = as.numeric(length(new_points))
  
  L = round(2 * D / 3 - 1)
  stepsize = dfHsl$V2[i] * 3
  scaling_value = dfFittingFactor$V2[i]
  stepsize_delta = (dfHsl_Gupta$V2[i] - dfHsl_filtAR$V2[i]) / 2 * 3

  init_points = new_points[1,1:(D - num_basis_spline)]
  spl_basis = new_points[1, (D - num_basis_spline + 1):D]
  
  fileInitialpointsNew = paste(dirInput, "initialpoint_", ID_NEW, "_", i, ".txt", sep ="")
  if(file.exists(fileInitialpointsNew)){
    invisible(file.remove(fileInitialpointsNew))
  }
  write.table(x = init_points, file = fileInitialpointsNew, append = T, row.names = F, col.names = F, quote = F)
  
  fileSplinebasisNew = paste(dirInput, "splinebasis_", ID_NEW, "_", i, ".txt", sep ="")
  if(file.exists(fileSplinebasisNew)){
    invisible(file.remove(fileSplinebasisNew))
  }
  write.table(x = spl_basis, file = fileSplinebasisNew, append = T, row.names = F, col.names = F, quote = F)

  inputfileTable <- rbind(c("model", model))
  inputfileTable <- rbind(inputfileTable, c("data", data))
  inputfileTable <- rbind(inputfileTable, c("method", method))
  inputfileTable <- rbind(inputfileTable, c("seed", seed))
  inputfileTable <- rbind(inputfileTable, c("iter_sampling", iter_sampling))
  inputfileTable <- rbind(inputfileTable, c("iter_burn_in", iter_burn_in))
  inputfileTable <- rbind(inputfileTable, c("integrator", integrator))
  inputfileTable <- rbind(inputfileTable, c("t_L", t_L))
  inputfileTable <- rbind(inputfileTable, c("L", L))
  inputfileTable <- rbind(inputfileTable, c("t_stepsize", t_stepsize))
  inputfileTable <- rbind(inputfileTable, c("stepsize", stepsize))
  inputfileTable <- rbind(inputfileTable, c("thinning", thinning))
  inputfileTable <- rbind(inputfileTable, c("t_Phi", t_Phi))
  inputfileTable <- rbind(inputfileTable, c("Phi", Phi))
  inputfileTable <- rbind(inputfileTable, c("scaling_value", scaling_value))
  inputfileTable <- rbind(inputfileTable, c("stepsize_delta", stepsize_delta))
  inputfileTable <- rbind(inputfileTable, c("iter_tune", iter_tune))
  inputfileTable <- rbind(inputfileTable, c("AR_target", AR_target))
  inputfileTable <- rbind(inputfileTable, c("delta_AR_target", delta_AR_target))
  inputfileTable <- rbind(inputfileTable, c("num_basis_spline", num_basis_spline))
  inputfileTable <- rbind(inputfileTable, c("spline_pol_degree", spline_pol_degree))
  inputfileTable <- rbind(inputfileTable, c("tfin_spl_bas", tfin_spl_bas))
  inputfileTable <- rbind(inputfileTable, c("num_out_cvodes", num_out_cvodes))
  inputfileTable <- rbind(inputfileTable, c("S0_prior", S0_prior))
  inputfileTable <- rbind(inputfileTable, c("S0_param_1", S0_param_2))
  inputfileTable <- rbind(inputfileTable, c("S0_param_2", S0_param_2))
  inputfileTable <- rbind(inputfileTable, c("E0_prior", E0_prior))
  inputfileTable <- rbind(inputfileTable, c("E0_param_1", E0_param_2))
  inputfileTable <- rbind(inputfileTable, c("E0_param_2", E0_param_2))
  inputfileTable <- rbind(inputfileTable, c("I0_prior", I0_prior))
  inputfileTable <- rbind(inputfileTable, c("I0_param_1", I0_param_2))
  inputfileTable <- rbind(inputfileTable, c("I0_param_2", I0_param_2))
  inputfileTable <- rbind(inputfileTable, c("alpha_prior", alpha_prior))
  inputfileTable <- rbind(inputfileTable, c("alpha_param_1", alpha_param_1))
  inputfileTable <- rbind(inputfileTable, c("alpha_param_2", alpha_param_2))
  inputfileTable <- rbind(inputfileTable, c("gamma_prior", gamma_prior))
  inputfileTable <- rbind(inputfileTable, c("gamma_param_1", gamma_param_1))
  inputfileTable <- rbind(inputfileTable, c("gamma_param_2", gamma_param_2))
  inputfileTable <- rbind(inputfileTable, c("phi_inv_prior", phi_inv_prior))
  inputfileTable <- rbind(inputfileTable, c("phi_inv_param_1", phi_inv_param_1))
  inputfileTable <- rbind(inputfileTable, c("phi_inv_param_2", phi_inv_param_2))
  inputfileTable <- rbind(inputfileTable, c("tau_prior", tau_prior))
  inputfileTable <- rbind(inputfileTable, c("tau_param_1", tau_param_1))
  inputfileTable <- rbind(inputfileTable, c("tau_param_2", tau_param_2))
  inputfileTable <- rbind(inputfileTable, c("num_comp", num_comp))
  inputfileTable <- rbind(inputfileTable, c("num_comp_E", num_comp_E))
  inputfileTable <- rbind(inputfileTable, c("do_under_report", do_under_report))
  inputfileTable <- rbind(inputfileTable, c("T0_under", T0_under))
  inputfileTable <- rbind(inputfileTable, c("U0_under", U0_under))
  inputfileTable <- rbind(inputfileTable, c("T1_under", T1_under))
  inputfileTable <- rbind(inputfileTable, c("U1_under", U1_under))
  inputfileTable <- rbind(inputfileTable, c("do_grad_desc", do_grad_desc))
  inputfileTable <- rbind(inputfileTable, c("learning_rate", learning_rate))
  inputfileTable <- rbind(inputfileTable, c("max_iter", max_iter))
  inputfileTable <- rbind(inputfileTable, c("gammaFixed", gammaFixed))
  inputfileTable <- rbind(inputfileTable, c("gamma_fixed_value", gamma_fixed_value))
  inputfileTable <- rbind(inputfileTable, c("gammaBounded", gammaBounded))
  inputfileTable <- rbind(inputfileTable, c("gammaUpper", gammaUpper))
  inputfileTable <- rbind(inputfileTable, c("gammaLower", gammaLower))
  inputfileTable <- rbind(inputfileTable, c("beta_prior", beta_prior))
  inputfileTable <- rbind(inputfileTable, c("beta_param_1", beta_param_1))
  inputfileTable <- rbind(inputfileTable, c("beta_param_2", beta_param_2))

  file_inputfile = paste(dirInput, "inputfile_", ID_NEW, "_", i, ".txt", sep = "")
  if(file.exists(file_inputfile)){
    invisible(file.remove(file_inputfile))
  }
  
  write.table(inputfileTable, file_inputfile, row.names = FALSE, col.names = FALSE, quote = F)
}

#THE END#######################################################################
cat("\n\tEnd of sampling_optimal_setting.R\n\n")