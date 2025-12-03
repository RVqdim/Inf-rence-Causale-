N  <- 250        # Nb individus
P  <- 28000      # Nb SNPs
K  <- 80         # Nb metabolites
P_causal <- 30   # Nb SNPs causaux



# Fréquences alléliques réalistes (entre 0.05 et 0.5)
allele_freqs <- runif(P, 0.05, 0.5)

simulate_snp <- function(p, N){
  # Génotype 0/1/2 avec HWE
  rbinom(N, 2, p)
}

G <- sapply(allele_freqs, simulate_snp, N = N)
G <- as.matrix(G)   # matrice SNP (N × P)



causal_idx <- sample(1:P, P_causal)
G_causal <- G[, causal_idx]



# Effets SNP → Métabolite
# - la plupart des SNPs n'affectent que quelques métabolites
# - effets modulables
beta <- matrix(0, nrow = P_causal, ncol = K)

for(i in 1:P_causal){
  affected <- sample(1:K, sample(1:5, 1))  # chaque SNP touche 1 à 5 métabolites
  beta[i, affected] <- rnorm(length(affected), 0.2, 0.1)
}


# Simuler des "voies métaboliques" corrélées
n_latent <- 5
latent <- matrix(rnorm(N * n_latent), N, n_latent)
loadings <- matrix(rnorm(K * n_latent, 0, .3), n_latent, K)

latent_effect <- latent %*% loadings   # (N × K)



# Contribution génétique causale
genetic_effect <- G_causal %*% beta

# Bruit individuel
noise <- matrix(rnorm(N * K, 0, 1), N, K)

# Métabolites simulés finaux
M <- genetic_effect #+ noise + latent_effect


colnames(G) <- paste0("SNP_", 1:P)
colnames(M) <- paste0("Met_", 1:K)


# Matrice estimée des effets causaux
beta_est <- matrix(0, nrow = ncol(G_causal), ncol = ncol(M))

for(j in 1:ncol(M)) {
  # On ajuste un modèle pour chaque métabolite
  fit <- lm(M[, j] ~ G_causal)
  # On récupère les coefficients sans l’intercept
  beta_est[, j] <- coef(fit)[-1]
}

colnames(beta_est) <- colnames(M)
rownames(beta_est) <- paste0("SNP_causal_", 1:ncol(G_causal))


# Corrélation globale
overall_cor <- cor(as.vector(beta_est), as.vector(beta))

overall_cor




# Choisir combien de composantes garder
r <- 80   # à modifier au besoin (ex: tester 1, 5, 10, 20)

# 1) ACP sur les métabolites
pca <- prcomp(M, center = TRUE, scale. = FALSE)

# 2) Reconstruction à partir des r premières composantes
scores_r <- pca$x[, 1:r, drop = FALSE]             # N x r
loadings_r <- pca$rotation[, 1:r, drop = FALSE]    # K x r
M_proj <- scores_r %*% t(loadings_r)               # N x K

# 3) Réestimation des betas avec M_proj
beta_est_proj <- matrix(0, nrow = ncol(G_causal), ncol = ncol(M_proj))

for(j in 1:ncol(M_proj)) {
  fit <- lm(M_proj[, j] ~ G_causal)
  beta_est_proj[, j] <- coef(fit)[-1]
}

# 4) Mesure de préservation causale
overall_cor_pca <- cor(as.vector(beta_est_proj), as.vector(beta))

overall_cor_pca




