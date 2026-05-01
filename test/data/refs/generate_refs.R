# Generates reference results from R's `survival` package for each dataset
# pulled by `download_data.jl` from
# https://github.com/vincentarelbundock/Rdatasets. R operates on the same
# `<name>_data.csv` files the Julia tests load, so the two sides see identical
# data.
#
# Run from the package root:
#   julia test/data/refs/download_data.jl
#   Rscript test/data/refs/generate_refs.R

suppressPackageStartupMessages(library(survival))

ref_dir <- "test/data/refs"

read_data <- function(name) {
    read.csv(file.path(ref_dir, paste0(name, "_data.csv")),
             stringsAsFactors = TRUE)
}

write_km_na <- function(name, time, status) {
    # Survival.jl rejects `t == 0` observations (`EventTable` throws), so the
    # caller is expected to filter them out; mirror that here for parity. None
    # of the current upstream datasets actually contain t = 0, but the filter
    # keeps the reference output stable if a future refresh introduces them.
    keep <- time > 0
    time <- time[keep]
    status <- status[keep]
    fit <- survfit(Surv(time, status) ~ 1, conf.type = "log-log")

    n <- fit$n.risk
    d <- fit$n.event
    aalen_var <- cumsum(ifelse(n > 0, d * (n - d) / n^3, 0))

    df <- data.frame(
        time     = fit$time,
        n_risk   = fit$n.risk,
        n_event  = fit$n.event,
        n_censor = fit$n.censor,
        surv     = fit$surv,
        std_err  = fit$std.err,
        lower    = fit$lower,
        upper    = fit$upper,
        cumhaz   = fit$cumhaz,
        chaz_se  = sqrt(aalen_var)
    )
    write.csv(df, file.path(ref_dir, paste0(name, "_km.csv")), row.names = FALSE)
    invisible(NULL)
}

write_cox <- function(name, fit) {
    s <- summary(fit)
    ci <- confint(fit)
    df <- data.frame(
        variable = rownames(s$coefficients),
        coef     = s$coefficients[, "coef"],
        se       = s$coefficients[, "se(coef)"],
        z        = s$coefficients[, "z"],
        pvalue   = s$coefficients[, "Pr(>|z|)"],
        lo95     = ci[, 1],
        hi95     = ci[, 2]
    )
    write.csv(df, file.path(ref_dir, paste0(name, "_cox.csv")), row.names = FALSE)
    meta <- data.frame(
        loglik      = fit$loglik[2],
        null_loglik = fit$loglik[1],
        n           = fit$n,
        nevent      = fit$nevent
    )
    write.csv(meta, file.path(ref_dir, paste0(name, "_cox_meta.csv")), row.names = FALSE)
    invisible(NULL)
}

# --- lung (== survival::cancer) ---
lung <- read_data("lung")
write_km_na("lung", lung$time, lung$status == 2)
lung_cc <- na.omit(lung[, c("time", "status", "age", "sex", "ph_ecog", "ph_karno", "wt_loss")])
fit_lung <- coxph(Surv(time, status == 2) ~ age + sex + ph_ecog + ph_karno + wt_loss,
                  data = lung_cc, ties = "efron")
write_cox("lung", fit_lung)

# --- leukemia ---
leukemia <- read_data("leukemia")
write_km_na("leukemia", leukemia$time, leukemia$status == 1)
fit_leuk <- coxph(Surv(time, status) ~ x, data = leukemia, ties = "efron")
write_cox("leukemia", fit_leuk)

# --- mgus ---
mgus <- read_data("mgus")
write_km_na("mgus", mgus$futime, mgus$death == 1)
mgus_cc <- na.omit(mgus[, c("futime", "death", "age", "sex", "alb", "creat", "hgb", "mspike")])
fit_mgus <- coxph(Surv(futime, death) ~ age + sex + alb + creat + hgb + mspike,
                  data = mgus_cc, ties = "efron")
write_cox("mgus", fit_mgus)

# --- nwtco ---
nwtco <- read_data("nwtco")
write_km_na("nwtco", nwtco$edrel, nwtco$rel == 1)
fit_nwtco <- coxph(
    Surv(edrel, rel) ~ factor(histol) + factor(stage) + age + factor(study),
    data = nwtco, ties = "efron"
)
write_cox("nwtco", fit_nwtco)

# --- ovarian ---
ovarian <- read_data("ovarian")
write_km_na("ovarian", ovarian$futime, ovarian$fustat == 1)
fit_ovarian <- coxph(Surv(futime, fustat) ~ age + ecog_ps + rx + resid_ds,
                     data = ovarian, ties = "efron")
write_cox("ovarian", fit_ovarian)

# --- pbc ---
pbc <- read_data("pbc")
write_km_na("pbc", pbc$time, pbc$status == 2)
pbc_cc <- na.omit(pbc[, c("time", "status", "age", "edema", "bili", "albumin",
                          "protime", "stage")])
fit_pbc <- coxph(
    Surv(time, status == 2) ~ age + edema + log(bili) + log(albumin) +
        log(protime) + stage,
    data = pbc_cc, ties = "efron"
)
write_cox("pbc", fit_pbc)

# --- stanford2 ---
stanford2 <- read_data("stanford2")
write_km_na("stanford2", stanford2$time, stanford2$status == 1)
stan_cc <- na.omit(stanford2[, c("time", "status", "age", "t5")])
fit_stan <- coxph(Surv(time, status) ~ age + t5, data = stan_cc, ties = "efron")
write_cox("stanford2", fit_stan)

# --- veteran ---
veteran <- read_data("veteran")
veteran$celltype <- factor(veteran$celltype,
                           levels = c("squamous", "smallcell", "adeno", "large"))
write_km_na("veteran", veteran$time, veteran$status == 1)
fit_vet <- coxph(
    Surv(time, status) ~ trt + celltype + karno + diagtime + age + prior,
    data = veteran, ties = "efron"
)
write_cox("veteran", fit_vet)

# --- kidney (treat recurrent events as independent) ---
kidney <- read_data("kidney")
kidney$disease <- factor(kidney$disease, levels = c("Other", "GN", "AN", "PKD"))
write_km_na("kidney", kidney$time, kidney$status == 1)
fit_kidney <- coxph(Surv(time, status) ~ age + sex + disease,
                    data = kidney, ties = "efron")
write_cox("kidney", fit_kidney)

# --- colon (death event subset) ---
colon <- read_data("colon")
colon$rx <- factor(colon$rx, levels = c("Obs", "Lev", "Lev+5FU"))
colon_death <- subset(colon, etype == 2)
write_km_na("colon", colon_death$time, colon_death$status == 1)
colon_cc <- na.omit(colon_death[, c("time", "status", "rx", "sex", "age",
                                    "obstruct", "perfor", "adhere", "nodes",
                                    "differ", "extent", "surg")])
fit_colon <- coxph(
    Surv(time, status) ~ rx + sex + age + obstruct + perfor + adhere + nodes +
        differ + extent + surg,
    data = colon_cc, ties = "efron"
)
write_cox("colon", fit_colon)

cat("Wrote reference files to ", ref_dir, "\n", sep = "")
